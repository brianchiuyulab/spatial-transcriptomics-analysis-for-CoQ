#!/usr/bin/env python
"""Create compact metadata from the authors' processed Dryad H5AD files."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


CELL_TYPES = (
    "B cells", "Dendritic cells", "Endothelial cells", "Erythrocytes", "FAPs",
    "Fusing Myocytes", "Monocytes/Macrophages", "MuSCs", "Myonuclei", "NK cells",
    "Neutrophils", "Pericytes and Smooth muscle cells", "Schwann and Neural/Glial cells",
    "T cells", "Tenocytes",
)


def celltype_score_frame(adata: ad.AnnData) -> pd.DataFrame:
    if all(name in adata.obs.columns for name in CELL_TYPES):
        return adata.obs.loc[:, CELL_TYPES].astype(float).copy()
    candidates = [
        "q05_cell_abundance_w_sf",
        "means_cell_abundance_w_sf",
        "cell_abundance_w_sf",
    ]
    for key in candidates:
        if key in adata.obsm:
            values = np.asarray(adata.obsm[key])
            names = adata.uns.get("mod", {}).get("factor_names") if isinstance(adata.uns.get("mod"), dict) else None
            if names is None:
                names = adata.uns.get("factor_names")
            if names is not None and values.shape[1] == len(names):
                frame = pd.DataFrame(values, index=adata.obs_names, columns=list(names))
                if all(name in frame.columns for name in CELL_TYPES):
                    return frame.loc[:, CELL_TYPES].astype(float)
    raise KeyError("Could not find the 15 deposited cell2location score columns.")


def extract(path: Path, sample: str, age_months: int) -> pd.DataFrame:
    adata = ad.read_h5ad(path, backed="r")
    obs = adata.obs.copy()
    obs.index = obs.index.astype(str)
    scores = celltype_score_frame(adata)
    scores.index = scores.index.astype(str)
    if "spatial" not in adata.obsm:
        raise KeyError(f"{path.name} does not contain obsm['spatial'] coordinates.")
    spatial = np.asarray(adata.obsm["spatial"])
    frame = pd.DataFrame({
        "barcode": obs.index,
        "sample": sample,
        "age_months": age_months,
        "x_um": spatial[:, 0],
        "y_um": spatial[:, 1],
    }, index=obs.index)
    frame = pd.concat([frame, scores], axis=1)
    frame["total_abundance"] = scores.sum(axis=1)
    frame["max_pred"] = scores.max(axis=1)
    frame["max_pred_celltype"] = scores.idxmax(axis=1)
    preferred = [
        "spatial_neighbors_100_true", "n_genes_by_counts", "log1p_n_genes_by_counts",
        "total_counts", "log1p_total_counts", "pct_counts_in_top_50_genes",
        "pct_counts_in_top_100_genes", "pct_counts_in_top_200_genes",
        "pct_counts_in_top_500_genes", "total_counts_mito", "log1p_total_counts_mito",
        "pct_counts_mito",
    ]
    preferred += [f"{name}_neighbors" for name in CELL_TYPES]
    preferred += ["Immune_neighbors", "injury_classification", "injury_neighbors", "Senescence_score"]
    for column in preferred:
        if column in obs.columns:
            frame[column] = obs[column].to_numpy()
    if "n_counts" in obs.columns:
        frame.insert(5, "n_counts", obs["n_counts"].to_numpy())
    elif "total_counts" in obs.columns:
        frame.insert(5, "n_counts", obs["total_counts"].to_numpy())
    return frame.reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--young", type=Path, required=True)
    parser.add_argument("--geriatric", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    combined = pd.concat([
        extract(args.young, "Young", 5),
        extract(args.geriatric, "Geriatric", 26),
    ], ignore_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output, index=False, compression="gzip")
    print(f"Wrote {len(combined):,} beads to {args.output}")


if __name__ == "__main__":
    main()
