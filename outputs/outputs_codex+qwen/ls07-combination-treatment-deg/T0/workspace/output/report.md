# Differential expression: Cisplatin_IC50 + CBD_IC50 combination vs DMSO

## Objective and design

This analysis identifies genes whose expression differs between cells treated
with the cisplatin + CBD combination at their respective IC50 concentrations
(group `Cisplatin_IC50_CBD_IC50`, samples 9-1, 9-2, 9-3) and vehicle-treated
controls (`DMSO`, samples 3-1, 3-2, 3-3). Only these six samples were used;
all other experimental groups were excluded. Fold changes are reported as
combination treatment relative to DMSO (numerator / denominator), with
`Group` as the sole design term.

## Methods

Raw counts were pre-filtered to genes with a raw count greater than 10 in at
least one of the six selected samples, reducing 63,677 genes to 18,029.
Differential expression was performed with PyDESeq2 0.5.0
(`DeseqDataSet` with `refit_cooks=True`; standard `DeseqStats` contrast
`["Group", "Cisplatin_IC50_CBD_IC50", "DMSO"]`), implementing the DESeq2
workflow: size-factor normalization, negative-binomial GLM with shrinkage-free
Wald testing, and Benjamini-Hochberg adjustment with independent filtering.
No Cook's-distance outliers required replacement. Because independent
filtering excludes low-mean genes from adjustment, 4,545 genes have no
adjusted p-value; these are preserved as null (empty CSV cells / JSON null),
never coerced to zero.

## Results

- Genes tested: 18,029
- Genes passing all criteria (padj < 0.05, |log2FoldChange| > 0.5,
  baseMean > 10): 555 - 356 upregulated, 199 downregulated in the combination
  treatment relative to DMSO.

Strongest associations (by padj) include GDA (log2FC 1.28), AHRR (2.13),
ALDH1A3 (1.62), ADGRF1 (0.83), CYP1A1 (2.03), and NR2F2 (0.97). The prominent
up-regulation of xenobiotic-response genes such as AHRR and CYP1A1 is
consistent with activation of stress-response and drug-metabolism programs
under combined cytotoxic challenge. Full per-gene statistics are in
`output/differential_expression.csv`; counts and the significant gene list are
in `output/summary.json`.

## Interpretation and limitations

These results are statistical associations between the combination treatment
and steady-state mRNA abundance in this cell-based experiment. Differential
expression does not by itself establish causation: the observed changes may be
downstream consequences of drug-induced stress, cell-cycle arrest, or shifts
in cell-state composition rather than direct transcriptional effects of the
compounds, and the DMSO comparison does not separate combination-specific
effects from responses to either single agent. With only three replicates per
condition, effect-size estimates for individual genes are uncertain, and genes
failing independent filtering (null padj) cannot be declared non-significant -
they simply lack an adjusted value. Replication in independent experiments,
including single-agent arms and orthogonal assays, is needed before any
causal or mechanistic claims are made.
