# Frozen composition rule

For each matrix, transpose to cells by genes, divide every cell by its library size, multiply by 10,000, and apply `log1p`. For each cell type in `MARKER_PANEL.tsv`, calculate the arithmetic mean across its listed marker genes. Assign each cell to the type with the largest mean; ties follow the row order in the marker panel. Do not remove cells or genes. Fractions use the matrix column count as denominator. The depleted call is the listed type with the smallest sample-2/sample-1 fraction ratio among types having at least 1% in sample 1.

This simple marker-score rule is frozen for reproducibility of the benchmark artifact; it is not a substitute for a full biological annotation workflow.
