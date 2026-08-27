# ATAC-seq organ-label swap analysis

## Call

**Swap detected: true — Cloaca and Stomach** (high confidence).

The unique top reciprocal marker-coherence gain is 1.2527; the runner-up is GallBladder–Stomach at 0.6748, a margin of 0.5780. Cloaca markers rise from 0.994 in their labeled column to 1.976 in Stomach; Stomach markers rise from 0.847 to 1.118 in Cloaca.

## Method and safeguards

Official GTF gene TSS coordinates were converted from whole chromosome arms to sequential `_1/_2/...` ATAC fragments using the supplied chromosome sizes. Counts in the TSS bin and its two neighbors were divided by full-table library totals, so the decision is not based on library size. Conserved organ-marker promoter profiles were scaled across samples, and every one of the 105 unordered pairs was scored by the improvement in reciprocal diagonal coherence. `sample_similarity.csv` contains all pairs, ordered by decreasing finite score.

## Limits

Marker sets are curated conserved vertebrate markers and are not a substitute for matched axolotl reference epigenomes. Confidence therefore reflects the unique internal reciprocal signal and its margin, not independent sample provenance. The decision rule requires a positive top gain and margin >0.25; otherwise the script returns `swap_detected=false`.
