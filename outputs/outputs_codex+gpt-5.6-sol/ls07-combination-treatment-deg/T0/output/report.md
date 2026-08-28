# Combination treatment versus DMSO: differential expression

## Scope and method

The analysis used only the six requested samples: DMSO 3-1, 3-2, 3-3 and `Cisplatin_IC50_CBD_IC50` 9-1, 9-2, 9-3. PyDESeq2 0.5.0 fitted `~Group` with `refit_cooks=True`; the standard contrast used `Cisplatin_IC50_CBD_IC50` as numerator and DMSO as denominator. No other condition, covariate, transcript lookup, sequence, haplotype, or genome-track data entered the model.

Genes were retained before fitting when any selected sample had a raw count greater than 10. This retained 18,029 of 63,677 genes.

## Thresholded findings

Passing genes simultaneously satisfy `padj < 0.05`, `abs(log2FoldChange) > 0.5`, and `baseMean > 10`. A total of 555 genes pass: 356 up and 199 down in combination treatment relative to DMSO. For 4,545 retained genes, `padj` is unavailable and remains an empty CSV cell, never zero.

The complete filtered result table is keyed by Ensembl ID; the supplied mapping contributes display names without dropping unmapped genes.

These estimates identify statistical associations between treatment group and expression in this experiment. They do not establish that the combination treatment caused a specific mechanism or phenotype.
