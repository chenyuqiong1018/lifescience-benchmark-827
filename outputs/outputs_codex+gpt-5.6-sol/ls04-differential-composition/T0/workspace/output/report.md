# Differential retinal cell composition

## Reproducible method

The matrices were used exactly as supplied: 36,601 genes x 6,295 cells in sample 1 and 36,601 genes x 5,004 cells in sample 2. No cells or genes were filtered. Each cell was divided by its full library size, scaled to 10,000, and transformed with `log1p`. For each marker-panel cell type, the arithmetic mean across its listed transformed markers was computed. The maximum score defined the label, and panel row order broke ties. Fractions use all matrix columns.

The frozen depleted-call rule was then applied without post hoc changes: retain types with sample-1 fraction >= 0.01 and choose the smallest sample-2/sample-1 fraction ratio.

## QC

- Sample 1: 9,246,637 stored entries; library size 205–45,431, median 1999.0; 0 empty cells.
- Sample 2: 9,953,348 stored entries; library size 204–48,129, median 3197.0; 0 empty cells.
- Every listed marker appeared exactly once, and parsed stored-entry counts matched both headers.

## Result

The depleted call is **horizontal cell**: 237/6,295 cells in sample 1 (fraction 0.037649) versus 49/5,004 in sample 2 (fraction 0.009792). The sample-2/sample-1 fraction ratio is 0.260092, meaning the sample-2 fraction is about 26.0% of the sample-1 fraction under the frozen annotation.

For scale only, cell-level binomial approximations give 95% Wilson intervals [0.033222, 0.042640] and [0.007415, 0.012921], and a log-ratio approximation [0.191656, 0.352964]. These are descriptive conditional calculations, not valid sample-level biological inference when cells share donors or preparations.

## Annotation evidence and uncertainty

Labels are supported only by the supplied marker panels after the fixed transformation. Sample 1 has 172 exact top-score ties and a median winning-score margin of 0.993560; sample 2 has 26 ties and median margin 2.176783. Required tie handling was retained.

The marker rule is deterministic but biologically limited. It does not account for doublets, ambient RNA, marker overlap, continuous states, batch effects, or donor-level replication. With one aggregate matrix per sample, uncertainty between biological replicates and causality cannot be established. The call should therefore be interpreted as composition depletion under the supplied scoring rule, not proof that the population was biologically eliminated.
