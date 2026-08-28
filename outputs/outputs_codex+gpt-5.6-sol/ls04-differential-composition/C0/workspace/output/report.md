# Retinal differential composition

## Frozen workflow

Both input files are integer, coordinate-format Matrix Market count matrices with genes in rows and cells in columns. No cell or gene was removed. For every cell, all gene counts contributed to its library size; counts were divided by that library size, multiplied by 10,000, and transformed with `log1p`. For each of the 16 marker-panel rows, the score is the arithmetic mean of the transformed values of exactly the listed markers. Each cell was assigned to the largest score, with marker-panel row order resolving exact ties.

Fractions use the full declared matrix column count. The depleted population is selected exactly as specified: among types with sample-1 fraction >= 1%, choose the smallest sample-2/sample-1 fraction ratio.

## Quality control

- Sample 1: 36,601 genes x 6,295 cells; 9,246,637 stored entries; 0 empty libraries. Library sizes range from 205 to 45431 counts (median 1999.0).
- Sample 2: 36,601 genes x 5,004 cells; 9,953,348 stored entries; 0 empty libraries. Library sizes range from 204 to 48129 counts (median 3197.0).
- All listed marker symbols were present exactly once in the shared gene table. Parsed entry counts matched both Matrix Market headers.

## Depleted population

**horizontal cell** is the frozen-rule depleted call. It changes from 237/6,295 cells (0.037649) in sample 1 to 49/5,004 cells (0.009792) in sample 2, for a sample-2/sample-1 fraction ratio of 0.260092. Full counts and fractions, including zero-count types, are in `composition.csv`.

## Annotation evidence and uncertainty

The evidence for every label is limited to relative expression of the marker sets in `MARKER_PANEL.tsv` after the frozen normalization. In sample 1, 172 cells had an exact top-score tie and the median top-versus-runner-up score margin was 0.993560; in sample 2, the corresponding values were 26 and 2.176783. Ties were retained and resolved by panel order, as required.

This deterministic marker rule is not a full biological annotation workflow. It does not model batch effects, ambient RNA, doublets, donor variability, uncertainty in marker specificity, or sampling uncertainty in the composition ratio. The depleted call is therefore a reproducible description under the supplied rule, not proof of biological loss or a causal effect. Confirmatory work would ordinarily inspect broader expression programs, technical covariates, replicate structure, and independent retinal annotations.
