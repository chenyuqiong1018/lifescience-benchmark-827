# Differential expression: Cisplatin_IC50_CBD_IC50 vs DMSO

## Objective

Identify genes whose expression differs between cells treated with the
cisplatin + CBD combination (each at its IC50) and the DMSO vehicle control,
using the frozen six-sample contrast from the BixBench `bix-43` capsule.

## Methods

- Samples: the three `Cisplatin_IC50_CBD_IC50` replicates (9-1, 9-2, 9-3) as
  numerator and the three `DMSO` replicates (3-1, 3-2, 3-3) as denominator.
  No samples from the other eight groups were used.
- Design: `~Group` as the only design term.
- Pre-fit filter: genes with raw count > 10 in at least one of the six
  selected samples were retained (18,029 of 63,677 genes; 45,648 dropped).
- Model: PyDESeq2 0.5.0 (`DeseqDataSet` with `refit_cooks=True`; standard
  `DeseqStats` contrast `["Group", "Cisplatin_IC50_CBD_IC50", "DMSO"]`),
  i.e., DESeq2's negative-binomial GLM with median-of-ratios size factors,
  shrunken dispersions, and Wald tests. No Cook's outliers required
  replacement.
- Significance rule: `padj < 0.05`, `|log2FoldChange| > 0.5`, and
  `baseMean > 10`. Adjusted p-values use Benjamini-Hochberg with DESeq2-style
  independent filtering on `baseMean`; genes excluded by that filter keep
  `padj` as null (empty CSV cell / JSON null), never zero.

## Results

- Genes tested: 18,029; genes with unavailable `padj`: 4,196.
- Passing genes: **557** (359 up-regulated, 198 down-regulated in the
  combination treatment relative to DMSO).
- Strongest associations (lowest `padj`): GDA, AHRR, ALDH1A3, ADGRF1, CYP1A1,
  NR2F2, SMAD3, CDKN1A, ANKRD37, FDXR. The direction is consistent with
  xenobiotic and DNA-damage stress responses (e.g., CYP1A1, AHRR, CDKN1A),
  which is biologically plausible for cisplatin exposure.

Full statistics are in `differential_expression.csv` (keyed by Ensembl ID,
with display gene names where mappable); run metadata are in `summary.json`.

## Association versus causation

These results are **statistical associations**, not causal claims. The model
estimates how strongly each gene's expression co-varies with treatment
assignment, after normalizing for sequencing depth. It does not by itself
show that the combination treatment mechanistically causes the observed
expression changes, nor that any single compound is responsible. Several
alternative explanations remain open: (i) the two drugs were applied jointly,
so individual drug effects and drug-drug interactions cannot be separated
from this design; (ii) with three replicates per group, low-abundance genes
have unstable estimates (4,196 genes were excluded by independent
filtering); (iii) unmeasured batch effects, differences in cell viability or
cell-cycle composition between treated and control wells, and library
preparation artifacts can induce treatment-correlated expression shifts; and
(iv) IC50-dosed conditions also reflect differential survival, so some
changes may be consequences of selection rather than regulation. Causal
interpretation would require additional experiments (single-agent arms,
dose-response time courses, randomized blocking, and orthogonal validation
such as qPCR or perturbation of candidate regulators).
