Done. All deliverables are in `output/`, verified by an independent reimplementation.

## Result
**The severely depleted population in sample 2 is the horizontal cell** — 3.76% of cells in sample 1 (237/6,295) vs 0.98% in sample 2 (49/5,004), a sample2/sample1 fraction ratio of **0.260**, the smallest among all types with ≥1% in sample 1 per the frozen rule. Raw-marker evidence agrees: LHX1⁺ cells drop 122→4, ONECUT1⁺ 275→103, ONECUT2⁺ 272→73. Astrocytes are a close second (ratio 0.274) and are flagged in the report's uncertainty section.

## Method
Followed `inputs/ANNOTATION_RULE.md` exactly: MTX → cells×genes, CP10K + log1p, per-cell mean marker score per panel type, argmax assignment with panel-order tie-break, no cell/gene removal, fractions over matrix column counts. QC confirmed matrix orientation against the 36,601-row gene table and that all 81 marker genes are present and uniquely mapped.

## Deliverables
- `output/composition.csv` — 32 rows (16 types × 2 samples) with sample, cell_type, n_cells, fraction; counts sum exactly to 6,295 / 5,004
- `output/depleted_call.json` — the call, criterion, full ratio ranking among eligible types, and ineligible (<1%) types
- `output/analysis.py` — deterministic pipeline reproducing all artifacts (plus `qc.json` with all QC metrics)
- `output/report.md` — QC tables (library size: medians 1,999 vs 3,197 UMIs; genes/cell 1,166 vs 1,672), normalization, annotation evidence, and uncertainty (assignment-margin stats: 6.9%/4.3% of cells with top-2 score gap <0.1; horizontal-vs-astrocyte near-tie; shared markers RLBP1/SLC1A3; no doublet filtering by design; depth differences absorbed by CP10K)

**Verification:** a separate dense float64 numpy reimplementation of the frozen rule reproduced all 32 cell counts exactly (PASS).