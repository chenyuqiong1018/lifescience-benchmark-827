# Differential expression: Cisplatin+CBD combination (IC50) vs DMSO

## Objective and design
We compared cells treated with the cisplatin + CBD combination at IC50
(`Cisplatin_IC50_CBD_IC50`, samples 9-1, 9-2, 9-3) against DMSO controls
(3-1, 3-2, 3-3), with log2 fold changes expressed as combination over DMSO.
Only these six samples were used; all other experimental groups were
excluded. The model contained a single design term (`~condition`, i.e. Group).

## Filtering and method
From 63,677 raw Ensembl genes, 18,029 had a raw count greater than 10 in at
least one of the six selected samples and were retained before fitting.
PyDESeq2 0.5.0 then ran the full DESeq2 pipeline (median-of-ratios
normalization, genewise/trend/MAP dispersion estimation, Wald test) with
`refit_cooks=True` and the standard `DeseqStats` contrast
`["condition", "Cisplatin_IC50_CBD_IC50", "DMSO"]`. Because of DESeq2
independent filtering, 4,196 genes have no adjusted p-value; these are kept
as null (empty CSV cells, JSON `null`), never converted to zero.

## Results
Of 18,029 tested genes, **557** passed padj < 0.05, |log2FoldChange| > 0.5,
and baseMean > 10: **359 up-regulated** and **198 down-regulated** under the
combination treatment.

Prominent up-regulated genes include:

- **CDKN1A** (p21; log2FC ~ 1.20), a DNA-damage response consistent with
  cisplatin genotoxic stress;
- **CYP1A1** (~2.03) and **AHRR** (~2.13), xenobiotic/aryl-hydrocarbon
  pathway genes consistent with cannabinoid exposure;
- **ALDH1A3** (~1.62), **GDA** (~1.28), **NR2F2**, **SMAD3**, **FDXR**.

The strongest down-regulated gene is **ANKRD37** (log2FC ~ -2.0).

The full table is `differential_expression.csv` (keyed by Ensembl ID, with
gene symbols); counts and thresholds are in `summary.json`; the reproducible
script is `analysis.py`.

## Reproducibility note
The installed anndata 0.13.3 automatically reshapes 1-D arrays stored in
`varm`/`obsm` to shape (n, 1), which breaks PyDESeq2 0.5.0. `analysis.py`
applies a minimal, documented compatibility shim restoring the historical
1-D behavior; no other library code was modified.

## Association vs causation
These results quantify **statistical association**: expression differences
are associated with treatment assignment in this cell-line experiment. They
do not establish that any listed gene **causes** sensitivity or resistance
to the combination. Treatment-induced mRNA changes may be downstream
consequences of stress, cell-cycle arrest, or shifting cell-state
composition rather than mediators of drug response, and causal regulators
need not change their own expression. Establishing causation would require
perturbation experiments (knockdown/overexpression, dose-response, rescue
assays) and validation in independent models. Additionally, all six samples
were processed together, so batch effects cannot be separated from
treatment effects here.
