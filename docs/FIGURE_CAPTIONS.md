# Figure captions

## Fig01. Spatial transcriptomics quality control

Numbers of raw beads, beads with at least 20 transcripts, and final retained beads after the published spatial-neighbor filter. The right panel summarizes median raw counts, median mitochondrial-transcript percentage, and the percentage of retained beads with raw *Coq8a* count of at least one.

## Fig02. Young spatial cell-type labels

Dominant cell2location label for each retained bead in the young tibialis anterior section at 5 days post injury. The label is the largest of the 15 row-normalized q05-derived cell2location values deposited by the source study.

## Fig03. Geriatric spatial cell-type labels

Dominant cell2location label for each retained bead in the geriatric tibialis anterior section at 5 days post injury, using the same rule and color scale as Fig02.

## Fig04. Young spatial distribution of *Coq8a*

Spatial distribution of *Coq8a*-positive beads, Fusing Myocytes-labelled beads, and MuSC-labelled beads in the young section. *Coq8a*-positive denotes a raw count of at least one. Injury-zone classifications and cell2location labels are from the source study.

## Fig05. Geriatric spatial distribution of *Coq8a*

Spatial distribution of *Coq8a*-positive beads, Fusing Myocytes-labelled beads, and MuSC-labelled beads in the geriatric section, using the same definitions as Fig04.

## Fig06. Association between *Coq8a* detection and cell labels

Odds ratios and 95% confidence intervals for *Coq8a* detection across selected dominant cell2location labels within each section's injury zone. P values were calculated by two-sided Fisher exact tests and adjusted separately within each section across all 15 deposited cell labels using the Benjamini-Hochberg procedure. Magenta stars denote adjusted *p* < 0.05; the complete 15-label results are provided in `tables/coq8a_celltype_association.tsv`.

## Fig07. Myogenic marker profiles

Detection frequency and mean expression of myogenic differentiation and fusion markers in Fusing Myocytes-labelled and other injury-zone beads. Dot size indicates the percentage of beads with at least one detected transcript. Color indicates mean expression after total-count normalization to 10,000 followed by `log1p` transformation.

## Fig08. *Coq8a* enrichment around Fusing Myocytes-labelled regions

Ratio of observed *Coq8a*-positive density around Fusing Myocytes-labelled focal beads to the mean density around randomly selected injury-zone focal beads. The focal bead itself was excluded from each neighborhood. The null distribution used 2,000 random selections at each radius. Magenta stars denote Benjamini-Hochberg-adjusted *p* < 0.05 across the 10 sample-by-radius tests.

## Fig09. *Coq8a* detection by injury region

Percentage of retained beads with raw *Coq8a* count of at least one inside and outside the injury zone in each spatial section.
