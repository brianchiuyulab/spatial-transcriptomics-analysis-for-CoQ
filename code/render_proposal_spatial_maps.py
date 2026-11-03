"""Render proposal panels 5E/F with a shared field in micrometres.

Run from any directory: python code/render_proposal_spatial_maps.py --output PATH
The four PNGs use the same physical field and reserve a separate strip for the scale bar.
"""
import argparse
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from spatial_plot_config import CELL_COLORS, SAMPLES

ROOT = Path(__file__).resolve().parents[1]


def positive_barcodes(path):
    """Read just the Coq8a column of the deposited CSC count matrix."""
    with h5py.File(path) as raw:
        matrix = raw["X"]
        if matrix.attrs["encoding-type"] != "csc_matrix":
            raise ValueError("Expected the deposited CSC count matrix")
        columns = raw["var/_index"].asstr()[:].tolist()
        gene = columns.index("Coq8a")
        lo, hi = matrix["indptr"][gene:gene + 2]
        positive = matrix["indices"][lo:hi][matrix["data"][lo:hi] >= 1]
        barcodes = raw["obs/_index"].asstr()[:]
        return set(barcodes), set(barcodes[positive])


def render(frame, span, kind, path):
    fig = plt.figure(figsize=(4, 4.4), dpi=360)
    ax = fig.add_axes([0, .4 / 4.4, 1, 4 / 4.4])
    x, y = frame.x_um, frame.y_um
    xmid, ymid = (x.min() + x.max()) / 2, (y.min() + y.max()) / 2
    ax.set(xlim=(xmid - span / 2, xmid + span / 2),
           ylim=(ymid + span / 2, ymid - span / 2))
    ax.set_axis_off()

    def dots(mask, color, size):
        ax.scatter(x[mask], y[mask], c=color, s=size, linewidths=0)

    if kind == "E":
        for cell, color in CELL_COLORS.items():
            dots(frame.max_pred_celltype == cell, color, 1.1 if cell == "Myonuclei" else 2.25)
    else:
        injured = frame.injury_classification == "injury"
        dots(~injured, "#E8E8E8", .65)
        dots(injured, "#C9C9C9", .8)
        dots(injured & (frame.max_pred_celltype == "Fusing Myocytes"), "#F28E2B", 1.9)
        musc = injured & (frame.max_pred_celltype == "MuSCs")
        ax.scatter(x[musc], y[musc], facecolors="none", edgecolors="#2F6DB0", s=7, linewidths=.36)
        pos = frame.Coq8a_detected
        ax.scatter(x[pos], y[pos], c="#7B2CBF", edgecolors="white", s=5.5, linewidths=.12)

    strip = fig.add_axes([0, 0, 1, .4 / 4.4])
    strip.set(xlim=(0, span), ylim=(0, 1))
    strip.axis("off")
    right = span * .57
    strip.plot([right - 500, right], [.47, .47], c="black", lw=2.1, solid_capstyle="butt")
    strip.text(right + span * .03, .47, "500 µm", fontsize=16, va="center")
    fig.savefig(path, dpi=360)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(ROOT / "data/GSE266933_compact_metadata.csv.gz",
                           usecols=["sample", "barcode", "x_um", "y_um",
                                    "max_pred_celltype", "injury_classification"])
    if set(metadata["sample"]) != set(SAMPLES):
        raise ValueError("Unexpected samples in the deposited metadata")
    spans = metadata.groupby("sample")[["x_um", "y_um"]].agg(lambda v: v.max() - v.min())
    span = spans.to_numpy().max() * 1.08
    for sample, frame in metadata.groupby("sample"):
        cfg = SAMPLES[sample]
        raw, positive = positive_barcodes(ROOT / "data" / cfg["file"])
        if (len(frame) != cfg["final_beads"] or not frame.barcode.is_unique
                or not set(frame.barcode) <= raw or not frame.max_pred_celltype.isin(CELL_COLORS).all()
                or not np.isfinite(frame[["x_um", "y_um"]]).all().all()):
            raise ValueError(f"{sample}: unexpected metadata or barcode mismatch")
        frame = frame.assign(Coq8a_detected=frame.barcode.isin(positive))
        for kind in ("E", "F"):
            render(frame, span, kind, args.output / f"{kind}_{sample}.png")


if __name__ == "__main__":
    main()
