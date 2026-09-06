#!/usr/bin/env python
"""Validate the final spatial release and its key numerical claims."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd
from PIL import Image


FIGURES = [
    "Fig01_QC.png", "Fig02_Young_celltype_map.png", "Fig03_Geriatric_celltype_map.png",
    "Fig04_Young_Coq8a_map.png", "Fig05_Geriatric_Coq8a_map.png",
    "Fig06_Coq8a_celltype_association.png", "Fig07_Fusing_marker_audit.png",
    "Fig08_Coq8a_fusion_neighborhood.png", "Fig09_Coq8a_injury_vs_outside.png",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-dir", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.release_dir.resolve()

    required_data = [
        "GSM8257020_Young_TA_ST_5dpi.h5ad",
        "GSM8257021_Geriatric_TA_ST_5dpi.h5ad",
        "GSE266933_compact_metadata.csv.gz",
    ]
    for name in required_data:
        path = root / "data" / name
        assert path.is_file() and path.stat().st_size > 0, f"Missing data file: {path}"

    for name in FIGURES:
        path = root / "figures" / name
        assert path.is_file() and path.stat().st_size > 50_000, f"Missing or small figure: {path}"
        with Image.open(path) as image:
            assert image.format == "PNG" and image.width >= 1500 and image.height >= 1000

    qc = pd.read_csv(root / "tables" / "dataset_qc.tsv", sep="\t").set_index("sample")
    assert int(qc.loc["Young", "final_beads"]) == 32104
    assert int(qc.loc["Geriatric", "final_beads"]) == 22920
    assert int(qc.loc["Young", "coq8a_positive"]) == 465
    assert int(qc.loc["Geriatric", "coq8a_positive"]) == 243

    association = pd.read_csv(root / "tables" / "coq8a_celltype_association.tsv", sep="\t")
    fusing = association.loc[association["cell_type"].eq("Fusing Myocytes")].set_index("sample")
    assert fusing.loc["Young", "p_adj_bh_within_section_15_labels"] < 0.05
    assert fusing.loc["Geriatric", "p_adj_bh_within_section_15_labels"] < 0.05
    assert fusing.loc["Young", "relative_enrichment"] > 1
    assert fusing.loc["Geriatric", "relative_enrichment"] > 1
    assert fusing["odds_ratio_plot"].gt(1).all()
    assert fusing["odds_ratio_ci_low"].gt(1).all()
    assert association[["odds_ratio_plot", "odds_ratio_ci_low", "odds_ratio_ci_high"]].notna().all().all()

    marker = pd.read_csv(root / "tables" / "fusing_marker_audit.tsv", sep="\t")
    assert marker["detection_fold"].gt(1).all()

    presentation = root / "presentation" / "spatial_analysis_summary_final.pptx"
    assert presentation.is_file() and presentation.stat().st_size > 100_000, f"Missing presentation: {presentation}"
    explanation = root / "presentation" / "spatial_analysis_explained_for_review_v2_final.pptx"
    assert explanation.is_file() and explanation.stat().st_size > 100_000, f"Missing explanation deck: {explanation}"

    manifest = root / "tables" / "release_manifest.tsv"
    records = []
    excluded_parts = {"data", ".git", "__pycache__"}
    for path in sorted(
        p for p in root.rglob("*")
        if p.is_file() and not excluded_parts.intersection(p.relative_to(root).parts)
    ):
        if path == manifest:
            continue
        records.append({"relative_path": path.relative_to(root).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(records).to_csv(manifest, sep="\t", index=False)
    print(f"Validated {len(FIGURES)} separate PNG figures and key numerical claims.")
    print(f"Manifest: {manifest}")


if __name__ == "__main__":
    main()
