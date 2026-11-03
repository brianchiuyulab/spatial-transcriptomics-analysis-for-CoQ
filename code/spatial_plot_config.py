"""Shared sample definitions and cell-type colours for GSE266933."""

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
