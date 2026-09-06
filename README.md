# GSE266933 skeletal-muscle spatial transcriptomics: Coq8a analysis

This release reproduces a focused analysis of **Coq8a** in tibialis anterior muscle 5 days after notexin injury. The source study is Walter *et al.*, *Nature Aging* (2024), DOI: [10.1038/s43587-024-00756-3](https://doi.org/10.1038/s43587-024-00756-3). Raw count matrices are from [GEO GSE266933](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266933); processed spatial annotations were deposited in [Dryad](https://doi.org/10.5061/dryad.kkwh70sbv). The authors' analysis code is available at [ntekasi/ST_MuSCs](https://github.com/ntekasi/ST_MuSCs).

## Analysis scope

- Samples: female C57BL/6J mice, young (5 months) and geriatric (26 months), one spatial section per age, 5 days post injury.
- Technology: Curio Seeker/Slide-seq, approximately 10 µm bead resolution.
- Cell labels: the authors' deposited cell2location abundance estimates; each bead is displayed using the cell type with the largest deposited abundance.
- Primary endpoint: whether Coq8a-positive injury-zone beads are enriched for the Fusing Myocytes label within each section.
- Spatial endpoint: Coq8a-positive bead density around Fusing Myocytes foci compared with random injury-zone foci.
- Coq8a-positive definition: raw `Coq8a` count ≥ 1 on a bead; no percentile threshold is used.

The release does not retrain cell2location, infer new cell types, or compare ages statistically. Each age has one biological spatial sample. Cell-type and neighborhood results are within-section associations.

## Reproduce

Place the three input files in `data/` and run:

```bash
python code/run_spatial_analysis.py --release-dir .
```

Required Python packages are listed in `requirements.txt`. Randomization is deterministic (`seed = 266933`). The main neighborhood analysis uses 2,000 permutations; spatial-block checks use 5,000 permutations.

If starting from the authors' processed Dryad H5AD files, generate the compact metadata file first:

```bash
python code/prepare_author_metadata.py \
  --young Young_TA_ST_5dpi_pp.h5ad \
  --geriatric Geriatric_TA_ST_5dpi_pp.h5ad \
  --output data/GSE266933_compact_metadata.csv.gz
```

## Release contents

- `code/`: executable analysis and validation scripts.
- `config/`: fixed analysis parameters.
- `docs/`: manuscript-ready methods, interpretation, and figure guide.
- `figures/`: one journal figure per PNG file; no combined contact sheet.
- `tables/`: exact numerical outputs and the published-parameter audit.
- `presentation/`: submission-style one-figure-per-slide summary and a separate Chinese analysis-explanation deck.
- `data/`: local source data; excluded from the public code repository because of file size and source licensing.

## Main conclusion

Within both sections, Coq8a-positive beads in the injury zone are enriched for the Fusing Myocytes label. Fusing regions also show higher detection of myogenic differentiation/fusion markers, including *Myog*, *Mymk*, and *Mymx*. These observations support a local fusion-associated transcriptional niche, but they do not demonstrate cell-autonomous Coq8a expression, causality, histologic fusion index, or myotube width.

The release's spatial maps are newly rendered from the deposited author coordinates, injury labels, and cell2location abundances together with GEO counts; they are not unchanged figure files exported by the original authors.
