#!/usr/bin/env python
"""Reproduce the GSE266933 Coq8a spatial-transcriptomics analysis.

The analysis uses GEO raw count matrices and the authors' deposited spatial
coordinates, injury-zone calls, and cell2location outputs. It does not
recluster beads or retrain cell2location.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.lines import Line2D
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import fisher_exact, rankdata
from statsmodels.stats.multitest import multipletests


SEED = 266933
RADII_UM = (25, 50, 75, 100, 150)
CELL_TYPES = (
    "B cells",
    "Dendritic cells",
    "Endothelial cells",
    "Erythrocytes",
    "FAPs",
    "Fusing Myocytes",
    "Monocytes/Macrophages",
    "MuSCs",
    "Myonuclei",
    "NK cells",
    "Neutrophils",
    "Pericytes and Smooth muscle cells",
    "Schwann and Neural/Glial cells",
    "T cells",
    "Tenocytes",
)
GENES = ("Coq8a", "Myh3", "Myh8", "Myog", "Mymk", "Mymx", "Ncam1", "Pax7", "Myod1")
SAMPLES = {
    "Young": {
        "file": "GSM8257020_Young_TA_ST_5dpi.h5ad",
        "age_months": 5,
        "raw_beads": 77642,
        "count20_beads": 32156,
        "final_beads": 32104,
        "retained_genes": 11927,
        "color": "#3E78B2",
    },
    "Geriatric": {
        "file": "GSM8257021_Geriatric_TA_ST_5dpi.h5ad",
        "age_months": 26,
        "raw_beads": 78514,
        "count20_beads": 22974,
        "final_beads": 22920,
        "retained_genes": 10148,
        "color": "#D64B45",
    },
}
CELL_COLORS = {
    "B cells": "#8E6C8A",
    "Dendritic cells": "#7F8C8D",
    "Endothelial cells": "#2A9D8F",
    "Erythrocytes": "#9B2226",
    "FAPs": "#E9C46A",
    "Fusing Myocytes": "#F28E2B",
    "Monocytes/Macrophages": "#8C564B",
    "MuSCs": "#3E78B2",
    "Myonuclei": "#B9B9B9",
    "NK cells": "#9467BD",
    "Neutrophils": "#BCBD22",
    "Pericytes and Smooth muscle cells": "#17BECF",
    "Schwann and Neural/Glial cells": "#59A14F",
    "T cells": "#E15759",
    "Tenocytes": "#76B7B2",
}


def dense1(x) -> np.ndarray:
    return x.toarray().ravel() if sparse.issparse(x) else np.asarray(x).ravel()


def finite(x) -> float | None:
    value = float(x)
    return value if np.isfinite(value) else None


def load_sample(data_dir: Path, metadata: pd.DataFrame, sample: str):
    cfg = SAMPLES[sample]
    raw_path = data_dir / cfg["file"]
    counts = ad.read_h5ad(raw_path)
    if counts.n_obs != cfg["raw_beads"]:
        raise ValueError(f"{sample}: unexpected raw bead count")

    raw_totals = dense1(counts.X.sum(axis=1))
    if int(np.sum(raw_totals >= 20)) != cfg["count20_beads"]:
        raise ValueError(f"{sample}: count >=20 checkpoint failed")

    frame = metadata.loc[metadata["sample"].eq(sample)].copy().set_index("barcode", drop=False)
    if len(frame) != cfg["final_beads"]:
        raise ValueError(f"{sample}: unexpected final metadata bead count")
    missing = frame.index.difference(counts.obs_names)
    if len(missing):
        raise ValueError(f"{sample}: metadata barcodes missing from GEO counts")

    final_raw = counts[frame.index].copy()
    counts_before_gene_filter = dense1(final_raw.X.sum(axis=1))
    detected_per_gene = dense1((final_raw.X > 0).sum(axis=0))
    keep_gene = detected_per_gene >= 10

    mt_mask = np.asarray([str(g).lower().startswith("mt-") for g in final_raw.var_names])
    mt_counts = dense1(final_raw.X[:, mt_mask].sum(axis=1))
    total_counts = dense1(final_raw.X.sum(axis=1))
    mt_pct = np.divide(
        100 * mt_counts,
        total_counts,
        out=np.zeros_like(total_counts, dtype=float),
        where=total_counts > 0,
    )

    normalized = final_raw[:, keep_gene].copy()
    if normalized.n_vars != cfg["retained_genes"]:
        raise ValueError(f"{sample}: gene >=10-bead checkpoint failed")
    sc.pp.normalize_total(normalized, target_sum=1e4)
    sc.pp.log1p(normalized)

    frame["raw_counts"] = counts_before_gene_filter
    frame["mt_pct_recomputed"] = mt_pct
    for gene in GENES:
        if gene in normalized.var_names:
            values = dense1(normalized[:, gene].X)
            frame[f"{gene}_log1p10k"] = values
            frame[f"{gene}_detected"] = values > 0
        else:
            frame[f"{gene}_log1p10k"] = np.nan
            frame[f"{gene}_detected"] = False

    qc = {
        "sample": sample,
        "age_months": cfg["age_months"],
        "raw_beads": cfg["raw_beads"],
        "beads_counts_ge_20": cfg["count20_beads"],
        "final_beads": len(frame),
        "spatial_neighbor_removed": cfg["count20_beads"] - len(frame),
        "retained_genes": normalized.n_vars,
        "median_raw_counts": np.median(total_counts),
        "median_mt_pct": np.median(mt_pct),
        "beads_mt_gt20_pct": np.mean(mt_pct > 20) * 100,
        "coq8a_positive": int(frame["Coq8a_detected"].sum()),
        "coq8a_detection_pct": frame["Coq8a_detected"].mean() * 100,
    }
    return frame.reset_index(drop=True), qc


def celltype_counts(frame: pd.DataFrame) -> list[dict]:
    rows = []
    for region, mask in (
        ("Full section", np.ones(len(frame), dtype=bool)),
        ("Injury", frame["injury_classification"].eq("injury").to_numpy()),
        ("Outside", frame["injury_classification"].eq("outside").to_numpy()),
    ):
        total = int(mask.sum())
        counts = frame.loc[mask, "max_pred_celltype"].value_counts()
        for cell_type in CELL_TYPES:
            n = int(counts.get(cell_type, 0))
            rows.append(
                {
                    "sample": frame["sample"].iloc[0],
                    "region": region,
                    "cell_type": cell_type,
                    "beads": n,
                    "percent": 100 * n / total,
                }
            )
    return rows


def coq_region_summary(frame: pd.DataFrame) -> list[dict]:
    injury = frame["injury_classification"].eq("injury").to_numpy()
    coq = frame["Coq8a_detected"].to_numpy()
    table = np.array(
        [
            [np.sum(coq & injury), np.sum(~coq & injury)],
            [np.sum(coq & ~injury), np.sum(~coq & ~injury)],
        ]
    )
    odds, p_value = fisher_exact(table, alternative="two-sided")
    rows = []
    for region, mask in (("Injury", injury), ("Outside", ~injury)):
        rows.append(
            {
                "sample": frame["sample"].iloc[0],
                "region": region,
                "beads": int(mask.sum()),
                "coq8a_positive": int(np.sum(coq & mask)),
                "coq8a_detection_pct": 100 * coq[mask].mean(),
                "mean_coq8a_log1p10k": frame.loc[mask, "Coq8a_log1p10k"].mean(),
                "injury_vs_outside_odds_ratio": odds,
                "injury_vs_outside_p": p_value,
            }
        )
    return rows


def coq_celltype_association(frame: pd.DataFrame) -> list[dict]:
    injury = frame.loc[frame["injury_classification"].eq("injury")].copy()
    coq = injury["Coq8a_detected"].to_numpy()
    rows = []
    for cell_type in CELL_TYPES:
        label = injury["max_pred_celltype"].eq(cell_type).to_numpy()
        table = np.array(
            [
                [np.sum(coq & label), np.sum(coq & ~label)],
                [np.sum(~coq & label), np.sum(~coq & ~label)],
            ]
        )
        odds, p_value = fisher_exact(table, alternative="two-sided")
        among_positive = table[0, 0] / table[0].sum()
        baseline = label.mean()
        rows.append(
            {
                "sample": frame["sample"].iloc[0],
                "cell_type": cell_type,
                "coq8a_positive_in_label": int(table[0, 0]),
                "coq8a_positive_total": int(table[0].sum()),
                "label_beads": int(label.sum()),
                "injury_beads": len(injury),
                "fraction_label_among_coq8a_positive": among_positive,
                "baseline_label_fraction": baseline,
                "relative_enrichment": among_positive / baseline,
                "odds_ratio": odds,
                "p_fisher": p_value,
            }
        )
    result = pd.DataFrame(rows)
    result["p_adj_bh_within_section_15_labels"] = multipletests(result["p_fisher"], method="fdr_bh")[1]
    return result.to_dict("records")


def fusing_marker_audit(frame: pd.DataFrame) -> list[dict]:
    injury = frame.loc[frame["injury_classification"].eq("injury")].copy()
    label = injury["max_pred_celltype"].eq("Fusing Myocytes")
    rows = []
    for gene in ("Myh3", "Myh8", "Myog", "Mymk", "Mymx", "Ncam1"):
        inside = injury.loc[label, f"{gene}_detected"].mean()
        outside = injury.loc[~label, f"{gene}_detected"].mean()
        rows.append(
            {
                "sample": frame["sample"].iloc[0],
                "gene": gene,
                "fusing_beads": int(label.sum()),
                "detection_fusing": inside,
                "detection_other_injury": outside,
                "detection_fold": inside / outside,
            }
        )
    return rows


def direct_musc_summary(frame: pd.DataFrame) -> dict:
    injury = frame["injury_classification"].eq("injury")
    musc = injury & frame["max_pred_celltype"].eq("MuSCs")
    positive = musc & frame["Coq8a_detected"]
    fusing = injury & frame["max_pred_celltype"].eq("Fusing Myocytes")
    distances = cKDTree(frame.loc[fusing, ["x_um", "y_um"]]).query(
        frame.loc[musc, ["x_um", "y_um"]], k=1
    )[0]
    return {
        "sample": frame["sample"].iloc[0],
        "injury_musc_beads": int(musc.sum()),
        "coq8a_positive_musc_beads": int(positive.sum()),
        "coq8a_detection_pct_musc": 100 * frame.loc[musc, "Coq8a_detected"].mean(),
        "median_distance_to_fusing_um": np.median(distances),
    }


def local_coq_enrichment(
    frame: pd.DataFrame,
    radius_um: int,
    n_permutations: int,
    rng: np.random.Generator,
    extra_filter: np.ndarray | None = None,
) -> dict:
    scope = frame["injury_classification"].eq("injury").to_numpy()
    if extra_filter is not None:
        scope &= extra_filter
    data = frame.loc[scope].copy()
    xy = data[["x_um", "y_um"]].to_numpy()
    tree = cKDTree(xy)
    detected = data["Coq8a_detected"].to_numpy(dtype=float)
    fusing = data["max_pred_celltype"].eq("Fusing Myocytes").to_numpy()
    local_density = np.full(len(data), np.nan)
    for i, neighbors in enumerate(tree.query_ball_point(xy, r=radius_um)):
        neighbors = [j for j in neighbors if j != i]
        if neighbors:
            local_density[i] = detected[neighbors].mean()
    valid = np.flatnonzero(np.isfinite(local_density))
    focal = fusing & np.isfinite(local_density)
    observed = local_density[focal].mean()
    k = int(focal.sum())
    null = np.empty(n_permutations)
    for i in range(n_permutations):
        null[i] = local_density[rng.choice(valid, size=k, replace=False)].mean()
    null_mean = null.mean()
    p_value = (np.sum(np.abs(null - null_mean) >= abs(observed - null_mean)) + 1) / (n_permutations + 1)
    return {
        "sample": frame["sample"].iloc[0],
        "radius_um": radius_um,
        "injury_beads": int(scope.sum()),
        "fusing_foci": k,
        "coq8a_positive_injury": int(detected.sum()),
        "observed_local_detection": observed,
        "random_foci_local_detection": null_mean,
        "observed_to_random_ratio": observed / null_mean,
        "p_permutation": p_value,
        "permutations": n_permutations,
    }


def neighborhood_composition(
    frame: pd.DataFrame,
    radius_um: int,
    n_permutations: int,
    rng: np.random.Generator,
) -> list[dict]:
    data = frame.loc[frame["injury_classification"].eq("injury")].copy()
    xy = data[["x_um", "y_um"]].to_numpy()
    labels = pd.Categorical(data["max_pred_celltype"], categories=CELL_TYPES).codes
    tree = cKDTree(xy)
    composition = np.zeros((len(data), len(CELL_TYPES)), dtype=np.float32)
    valid = np.zeros(len(data), dtype=bool)
    for i, neighbors in enumerate(tree.query_ball_point(xy, r=radius_um)):
        neighbors = [j for j in neighbors if j != i]
        if neighbors:
            composition[i] = np.bincount(labels[neighbors], minlength=len(CELL_TYPES)) / len(neighbors)
            valid[i] = True
    focal = data["max_pred_celltype"].eq("Fusing Myocytes").to_numpy() & valid
    observed = composition[focal].mean(axis=0)
    eligible = np.flatnonzero(valid)
    null = np.empty((n_permutations, len(CELL_TYPES)), dtype=np.float32)
    for i in range(n_permutations):
        null[i] = composition[rng.choice(eligible, size=int(focal.sum()), replace=False)].mean(axis=0)
    null_mean = null.mean(axis=0)
    rows = []
    for j, cell_type in enumerate(CELL_TYPES):
        p_value = (np.sum(np.abs(null[:, j] - null_mean[j]) >= abs(observed[j] - null_mean[j])) + 1) / (n_permutations + 1)
        rows.append(
            {
                "sample": frame["sample"].iloc[0],
                "radius_um": radius_um,
                "neighbor_cell_type": cell_type,
                "observed_fraction": observed[j],
                "random_foci_fraction": null_mean[j],
                "observed_to_random_ratio": observed[j] / null_mean[j] if null_mean[j] > 0 else np.nan,
                "p_permutation": p_value,
                "permutations": n_permutations,
            }
        )
    return rows


def spearman_permutation(x: np.ndarray, y: np.ndarray, n_permutations: int, rng: np.random.Generator):
    xr = rankdata(x).astype(float)
    yr = rankdata(y).astype(float)
    xr -= xr.mean()
    yr -= yr.mean()
    denom = np.sqrt(np.sum(xr**2) * np.sum(yr**2))
    if denom == 0:
        return np.nan, np.nan
    observed = np.dot(xr, yr) / denom
    exceed = 0
    for _ in range(n_permutations):
        value = np.dot(xr, rng.permutation(yr)) / denom
        exceed += abs(value) >= abs(observed)
    return observed, (exceed + 1) / (n_permutations + 1)


def spatial_block_association(
    frame: pd.DataFrame,
    block_um: int,
    n_permutations: int,
    rng: np.random.Generator,
) -> dict:
    data = frame.loc[frame["injury_classification"].eq("injury")].copy()
    data["block_x"] = np.floor((data["x_um"] - data["x_um"].min()) / block_um).astype(int)
    data["block_y"] = np.floor((data["y_um"] - data["y_um"].min()) / block_um).astype(int)
    data["fusing_fraction"] = data["max_pred_celltype"].eq("Fusing Myocytes").astype(float)
    data["coq8a_detection"] = data["Coq8a_detected"].astype(float)
    data["coq8a_expression"] = data["Coq8a_log1p10k"].astype(float)
    blocks = (
        data.groupby(["block_x", "block_y"], observed=True)
        .agg(
            beads=("coq8a_detection", "size"),
            fusing_fraction=("fusing_fraction", "mean"),
            coq8a_detection=("coq8a_detection", "mean"),
            coq8a_expression=("coq8a_expression", "mean"),
        )
        .query("beads >= 10")
        .reset_index(drop=True)
    )
    rho_detection, p_detection = spearman_permutation(
        blocks["fusing_fraction"].to_numpy(),
        blocks["coq8a_detection"].to_numpy(),
        n_permutations,
        rng,
    )
    rho_expression, p_expression = spearman_permutation(
        blocks["fusing_fraction"].to_numpy(),
        blocks["coq8a_expression"].to_numpy(),
        n_permutations,
        rng,
    )
    return {
        "sample": frame["sample"].iloc[0],
        "block_um": block_um,
        "minimum_beads_per_block": 10,
        "blocks": len(blocks),
        "rho_detection": rho_detection,
        "p_detection_permutation": p_detection,
        "rho_expression": rho_expression,
        "p_expression_permutation": p_expression,
        "permutations": n_permutations,
    }


def add_scale_bar(ax, length_um: int = 500, side: str = "left"):
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    if side == "right":
        left = max(x0, x1) - 0.06 * abs(x1 - x0) - length_um
    else:
        left = min(x0, x1) + 0.06 * abs(x1 - x0)
    bottom = max(y0, y1) - 0.07 * abs(y1 - y0)
    ax.plot([left, left + length_um], [bottom, bottom], color="black", lw=2.2, solid_capstyle="butt")
    ax.text(left + length_um / 2, bottom - 0.025 * abs(y1 - y0), f"{length_um} µm", ha="center", va="top", fontsize=10)


def plot_qc(qc: pd.DataFrame, output: Path):
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.8), gridspec_kw={"width_ratios": [1.25, 1]})
    stages = ["Raw beads", "Counts ≥20", "Spatial QC pass"]
    x = np.arange(len(stages))
    width = 0.34
    for i, sample in enumerate(("Young", "Geriatric")):
        row = qc.loc[qc["sample"].eq(sample)].iloc[0]
        values = [row["raw_beads"], row["beads_counts_ge_20"], row["final_beads"]]
        axes[0].bar(x + (i - 0.5) * width, values, width, color=SAMPLES[sample]["color"], label=sample)
    axes[0].set_xticks(x, stages)
    axes[0].set_ylabel("Beads")
    axes[0].set_title("Published bead-level QC")
    axes[0].legend(frameon=False)

    metrics = ["median_raw_counts", "median_mt_pct", "coq8a_detection_pct"]
    labels = ["Median counts", "Median mtRNA (%)", "Coq8a+ beads (%)"]
    xpos = np.arange(len(metrics))
    for offset, sample in ((-0.10, "Young"), (0.10, "Geriatric")):
        row = qc.loc[qc["sample"].eq(sample)].iloc[0]
        values = [row[m] for m in metrics]
        axes[1].scatter(xpos + offset, values, s=70, color=SAMPLES[sample]["color"], label=sample)
        for xi, value in zip(xpos + offset, values):
            axes[1].text(xi, value, f" {value:.1f}", va="center", fontsize=9)
    axes[1].set_xticks(xpos, labels, rotation=18, ha="right")
    axes[1].set_title("Final spatial beads")
    axes[1].set_ylim(bottom=0)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=10)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_celltype_map(frame: pd.DataFrame, sample: str, output: Path):
    fig, ax = plt.subplots(figsize=(7.8, 7.0))
    for cell_type in CELL_TYPES:
        subset = frame["max_pred_celltype"].eq(cell_type)
        size = 3.0 if cell_type != "Myonuclei" else 1.4
        ax.scatter(frame.loc[subset, "x_um"], frame.loc[subset, "y_um"], s=size, color=CELL_COLORS[cell_type], rasterized=True)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.axis("off")
    add_scale_bar(ax, side="right")
    ax.set_title(f"{sample} TA, 5 dpi: author cell2location labels", fontsize=17, pad=12)
    handles = [Line2D([0], [0], marker="o", linestyle="", color=CELL_COLORS[c], markersize=6, label=c) for c in CELL_TYPES]
    ax.legend(handles=handles, frameon=False, fontsize=9.2, bbox_to_anchor=(1.01, 0.5), loc="center left", handletextpad=0.4)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_coq_map(frame: pd.DataFrame, sample: str, output: Path):
    fig, ax = plt.subplots(figsize=(7.4, 7.0))
    injury = frame["injury_classification"].eq("injury")
    outside = ~injury
    fusing = injury & frame["max_pred_celltype"].eq("Fusing Myocytes")
    musc = injury & frame["max_pred_celltype"].eq("MuSCs")
    coq = frame["Coq8a_detected"]
    ax.scatter(frame.loc[outside, "x_um"], frame.loc[outside, "y_um"], s=1.0, color="#E8E8E8", rasterized=True)
    ax.scatter(frame.loc[injury, "x_um"], frame.loc[injury, "y_um"], s=1.3, color="#C9C9C9", rasterized=True)
    ax.scatter(frame.loc[fusing, "x_um"], frame.loc[fusing, "y_um"], s=4.0, color="#F28E2B", rasterized=True)
    ax.scatter(frame.loc[musc, "x_um"], frame.loc[musc, "y_um"], s=17, facecolors="none", edgecolors="#2F6DB0", linewidths=0.65)
    ax.scatter(frame.loc[coq, "x_um"], frame.loc[coq, "y_um"], s=12, color="#7B2CBF", edgecolors="white", linewidths=0.3, rasterized=True)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.axis("off")
    add_scale_bar(ax)
    ax.set_title(f"{sample} TA, 5 dpi: Coq8a and fusion-associated niches", fontsize=17, pad=12)
    handles = [
        Line2D([0], [0], marker="o", linestyle="", color="#7B2CBF", markersize=7, label="Coq8a+ bead"),
        Line2D([0], [0], marker="o", linestyle="", color="#F28E2B", markersize=7, label="Fusing Myocytes label"),
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor="none", markeredgecolor="#2F6DB0", markersize=7, label="MuSC label"),
        Line2D([0], [0], marker="o", linestyle="", color="#C9C9C9", markersize=7, label="Other injury bead"),
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        fontsize=9.5,
        bbox_to_anchor=(1.01, 0.10),
        loc="lower left",
        handletextpad=0.45,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_celltype_association(table: pd.DataFrame, output: Path):
    selected = [
        "Fusing Myocytes", "Myonuclei", "MuSCs", "FAPs",
        "Monocytes/Macrophages", "Dendritic cells", "Endothelial cells", "Tenocytes",
    ]
    names = ["Fusing Myocytes", "Myonuclei", "MuSC", "FAP", "Monocytes/Macrophages", "Dendritic", "Endothelial", "Tenocytes"]
    y = np.arange(len(selected))
    fig, ax = plt.subplots(figsize=(8.0, 5.7))
    for offset, sample in ((-0.12, "Young"), (0.12, "Geriatric")):
        data = table.loc[table["sample"].eq(sample)].set_index("cell_type").reindex(selected)
        values = np.log2(data["relative_enrichment"].replace(0, np.nan))
        ax.scatter(values, y + offset, s=58, color=SAMPLES[sample]["color"], label=sample, zorder=3)
        sig = data["p_adj_bh_within_section_15_labels"].lt(0.05).to_numpy()
        ax.scatter(values[sig], (y + offset)[sig], s=115, facecolors="none", edgecolors=SAMPLES[sample]["color"], linewidths=1.5)
    ax.axvline(0, color="#777777", lw=1, ls="--")
    ax.set_yticks(y, names)
    ax.invert_yaxis()
    ax.set_xlabel(r"log$_2$ enrichment among Coq8a-positive injury beads")
    ax.set_title("Coq8a-positive injury beads are enriched for Fusing Myocytes", loc="left", fontsize=15)
    ax.legend(frameon=False, loc="upper right")
    ax.text(0.99, 0.02, r"Open ring: BH-adjusted $p$ < 0.05", transform=ax.transAxes, ha="right", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=11)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_marker_audit(table: pd.DataFrame, output: Path):
    genes = ["Myh3", "Myh8", "Myog", "Mymk", "Mymx", "Ncam1"]
    y = np.arange(len(genes))
    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    for offset, sample in ((-0.11, "Young"), (0.11, "Geriatric")):
        data = table.loc[table["sample"].eq(sample)].set_index("gene").reindex(genes)
        ax.scatter(np.log2(data["detection_fold"]), y + offset, s=65, color=SAMPLES[sample]["color"], label=sample)
    ax.axvline(0, color="#777777", lw=1, ls="--")
    ax.set_yticks(y, genes)
    ax.invert_yaxis()
    ax.set_xlabel(r"log$_2$ detection fold, Fusing label vs other injury beads")
    ax.set_title("Fusing Myocytes regions express myogenic differentiation markers", loc="left", fontsize=15)
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=11)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_local_enrichment(table: pd.DataFrame, output: Path):
    fig, ax = plt.subplots(figsize=(7.6, 5.1))
    for sample in ("Young", "Geriatric"):
        data = table.loc[table["sample"].eq(sample)].sort_values("radius_um")
        ax.plot(data["radius_um"], data["observed_to_random_ratio"], marker="o", ms=6, lw=2, color=SAMPLES[sample]["color"], label=sample)
        sig = data["p_adj_bh_across_10_radius_tests"].lt(0.05)
        ax.scatter(data.loc[sig, "radius_um"], data.loc[sig, "observed_to_random_ratio"], s=105, facecolors="none", edgecolors=SAMPLES[sample]["color"], linewidths=1.5)
    ax.axhline(1, color="#777777", lw=1, ls="--")
    ax.set_xlabel("Radius around Fusing Myocytes foci (µm)")
    ax.set_ylabel("Coq8a-positive density, observed / random foci")
    ax.set_title("Coq8a is locally enriched around transcriptionally defined fusion foci", loc="left", fontsize=15)
    ax.legend(frameon=False)
    ax.text(0.99, 0.02, r"Open ring: BH-adjusted $p$ < 0.05", transform=ax.transAxes, ha="right", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=11)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_region_detection(table: pd.DataFrame, output: Path):
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    for sample in ("Young", "Geriatric"):
        data = table.loc[table["sample"].eq(sample)].set_index("region").reindex(["Outside", "Injury"])
        ax.plot([0, 1], data["coq8a_detection_pct"], marker="o", ms=7, lw=2, color=SAMPLES[sample]["color"], label=sample)
    ax.set_xticks([0, 1], ["Outside injury zone", "Injury zone"])
    ax.set_ylabel("Coq8a-positive beads (%)")
    ax.set_title("The injury zone is not globally Coq8a-high", loc="left", fontsize=15)
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=11)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_table(rows, path: Path):
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, float_format="%.8g")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-dir", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--permutations", type=int, default=2000)
    parser.add_argument("--block-permutations", type=int, default=5000)
    args = parser.parse_args()

    release = args.release_dir.resolve()
    data_dir = release / "data"
    figure_dir = release / "figures"
    table_dir = release / "tables"
    figure_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(data_dir / "GSE266933_compact_metadata.csv.gz")
    required = {"barcode", "sample", "x_um", "y_um", "max_pred_celltype", "injury_classification"}
    if not required.issubset(metadata.columns):
        raise ValueError(f"Missing metadata fields: {sorted(required - set(metadata.columns))}")

    rng = np.random.default_rng(SEED)
    frames = {}
    qc_rows = []
    count_rows = []
    region_rows = []
    celltype_rows = []
    marker_rows = []
    direct_rows = []
    for sample in ("Young", "Geriatric"):
        frame, qc = load_sample(data_dir, metadata, sample)
        frames[sample] = frame
        qc_rows.append(qc)
        count_rows.extend(celltype_counts(frame))
        region_rows.extend(coq_region_summary(frame))
        celltype_rows.extend(coq_celltype_association(frame))
        marker_rows.extend(fusing_marker_audit(frame))
        direct_rows.append(direct_musc_summary(frame))

    local_rows = []
    composition_rows = []
    for sample, frame in frames.items():
        for radius in RADII_UM:
            local_rows.append(local_coq_enrichment(frame, radius, args.permutations, rng))
            composition_rows.extend(neighborhood_composition(frame, radius, args.permutations, rng))

    local_table = pd.DataFrame(local_rows)
    local_table["p_adj_bh_across_10_radius_tests"] = multipletests(local_table["p_permutation"], method="fdr_bh")[1]
    composition_table = pd.DataFrame(composition_rows)
    composition_table["p_adj_bh_across_150_composition_tests"] = multipletests(composition_table["p_permutation"], method="fdr_bh")[1]

    block_rows = []
    for sample, frame in frames.items():
        for block_um in (100, 150, 200, 250):
            block_rows.append(spatial_block_association(frame, block_um, args.block_permutations, rng))
    block_table = pd.DataFrame(block_rows)
    block_table["p_adj_detection_bh_across_8_tests"] = multipletests(block_table["p_detection_permutation"], method="fdr_bh")[1]
    block_table["p_adj_expression_bh_across_8_tests"] = multipletests(block_table["p_expression_permutation"], method="fdr_bh")[1]

    sensitivity_rows = []
    for sample, frame in frames.items():
        filters = {
            "Published final beads": np.ones(len(frame), dtype=bool),
            "Raw counts ≥30": frame["raw_counts"].ge(30).to_numpy(),
            "Raw counts ≥50": frame["raw_counts"].ge(50).to_numpy(),
            "mtRNA ≤20%": frame["mt_pct_recomputed"].le(20).to_numpy(),
            "mtRNA ≤25%": frame["mt_pct_recomputed"].le(25).to_numpy(),
        }
        for label, mask in filters.items():
            row = local_coq_enrichment(frame, 50, max(1000, args.permutations // 2), rng, mask)
            row["qc_rule"] = label
            sensitivity_rows.append(row)

    qc_table = pd.DataFrame(qc_rows)
    counts_table = pd.DataFrame(count_rows)
    region_table = pd.DataFrame(region_rows)
    celltype_table = pd.DataFrame(celltype_rows)
    marker_table = pd.DataFrame(marker_rows)

    write_table(qc_table, table_dir / "dataset_qc.tsv")
    write_table(counts_table, table_dir / "celltype_counts.tsv")
    write_table(region_table, table_dir / "coq8a_region_detection.tsv")
    write_table(celltype_table, table_dir / "coq8a_celltype_association.tsv")
    write_table(marker_table, table_dir / "fusing_marker_audit.tsv")
    write_table(direct_rows, table_dir / "direct_musc_endpoint.tsv")
    write_table(local_table, table_dir / "coq8a_fusion_neighborhood.tsv")
    write_table(composition_table, table_dir / "fusion_neighborhood_composition.tsv")
    write_table(block_table, table_dir / "spatial_block_association.tsv")
    write_table(sensitivity_rows, table_dir / "coq8a_neighborhood_qc_sensitivity.tsv")

    plt.rcParams.update({"font.family": "Arial", "font.size": 11, "axes.labelsize": 12, "axes.titlesize": 15, "axes.linewidth": 1})
    plot_qc(qc_table, figure_dir / "Fig01_QC.png")
    plot_celltype_map(frames["Young"], "Young", figure_dir / "Fig02_Young_celltype_map.png")
    plot_celltype_map(frames["Geriatric"], "Geriatric", figure_dir / "Fig03_Geriatric_celltype_map.png")
    plot_coq_map(frames["Young"], "Young", figure_dir / "Fig04_Young_Coq8a_map.png")
    plot_coq_map(frames["Geriatric"], "Geriatric", figure_dir / "Fig05_Geriatric_Coq8a_map.png")
    plot_celltype_association(celltype_table, figure_dir / "Fig06_Coq8a_celltype_association.png")
    plot_marker_audit(marker_table, figure_dir / "Fig07_Fusing_marker_audit.png")
    plot_local_enrichment(local_table, figure_dir / "Fig08_Coq8a_fusion_neighborhood.png")
    plot_region_detection(region_table, figure_dir / "Fig09_Coq8a_injury_vs_outside.png")

    print(f"Completed GSE266933 spatial release at {release}")


if __name__ == "__main__":
    main()
