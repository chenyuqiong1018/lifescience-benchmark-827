# Differential-expression report

PyDESeq2 0.5.0 was used to compare `Cisplatin_IC50_CBD_IC50` (numerator; samples 9-1, 9-2, 9-3) with `DMSO` (denominator; samples 3-1, 3-2, 3-3). No other groups were used. The sole design term was `Group` (`~Group`), `refit_cooks=True`, and the standard `DeseqStats` contrast was `['Group', 'Cisplatin_IC50_CBD_IC50', 'DMSO']`.

Before model fitting, genes were retained when at least one of the six selected samples had a raw count greater than 10. This retained 18,029 of 63,677 input genes. A gene passed only when adjusted p-value < 0.05, absolute log2 fold change > 0.5, and base mean > 10. 555 genes passed: 356 had positive and 199 had negative log2 fold change. Adjusted p-values were unavailable for 4,545 retained genes and remain empty in the CSV and null in the JSON summary; they were not treated as zero.

The supplied Ensembl-to-gene-name table was used only to add labels; statistical rows remain keyed by Ensembl gene ID. No external biomarker or transcript database was queried. These results identify expression changes statistically associated with the treatment contrast under this experiment. They do not establish that the combination treatment causally regulates any individual gene, nor do they isolate interaction effects between the two agents.
