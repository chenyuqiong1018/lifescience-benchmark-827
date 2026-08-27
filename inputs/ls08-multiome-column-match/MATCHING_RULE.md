# Frozen multiome matching rule

Use `ensembl112_gene_coordinates.tsv` (Ensembl Genes release 112). Retain gene symbols that occur exactly once in the annotation, are present in the RNA matrix, and lie on chromosomes 1–22, X, or Y. Map each retained gene's strand-aware transcription start site to its containing 10 kb ATAC bin. Apply `log1p` to both RNA TPM and the mapped ATAC-bin value, select the 2,000 mapped genes with highest variance across RNA columns, compute the 8×8 Pearson correlation matrix across genes, and choose the maximum-total-correlation one-to-one assignment with the Hungarian algorithm. `runner_up_score` is the second-highest row-wise correlation. Population and column identifiers remain the supplied strings.

This evaluator-frozen rule makes the local full-artifact extension reproducible; CompBioBench itself asks only for the hidden permutation.
