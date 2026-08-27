# Combination-treatment differential expression

## Analysis

PyDESeq2 0.5.0 compared `Cisplatin_IC50_CBD_IC50` (numerator) with `DMSO` (denominator). The model used only `Group` (`~Group`) and exactly three replicates from each group: DMSO 3-1, 3-2, 3-3 and combination treatment 9-1, 9-2, 9-3. No other layout samples entered the count matrix or metadata.

Before fitting, genes were retained when at least one of these six samples had a raw count greater than 10. This retained 18,029 of 63,677 input genes. The fit used `refit_cooks=True` and the standard `DeseqStats` contrast `Group, Cisplatin_IC50_CBD_IC50, DMSO`.

## Results

A gene passes only when all three strict criteria hold: `padj < 0.05`, `abs(log2FoldChange) > 0.5`, and `baseMean > 10`. There are 555 passing genes: 356 higher and 199 lower in combination treatment relative to DMSO. Adjusted p-values are unavailable for 4,545 retained genes; these remain empty in the CSV and are counted explicitly rather than converted to zero.

`differential_expression.csv` contains all filtered-gene statistics keyed by Ensembl ID, with capsule-supplied display names where available, plus pass and direction fields.

These differential-expression results are statistical associations for this experiment. They do not by themselves establish that the combination treatment caused a particular molecular mechanism or downstream phenotype.
