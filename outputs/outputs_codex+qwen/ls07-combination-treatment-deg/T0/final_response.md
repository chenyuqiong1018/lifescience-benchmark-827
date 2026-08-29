Done. All four deliverables are in `output/` and validated.

**Pipeline** (`output/analysis.py`)
- Selected exactly the 6 samples: `Cisplatin_IC50_CBD_IC50` (9-1/9-2/9-3, numerator) vs `DMSO` (3-1/3-2/3-3, denominator); all other groups excluded.
- Prefilter: kept genes with raw count > 10 in ≥1 of the six samples → 63,677 → **18,029 genes**, all tested.
- PyDESeq2 0.5.0: `DeseqDataSet(design="~Group", refit_cooks=True)` + standard `DeseqStats` contrast `["Group", "Cisplatin_IC50_CBD_IC50", "DMSO"]`. One environment fix was needed: anndata ≥0.13 reshapes 1-D arrays in `varm` to 2-D, which breaks PyDESeq2 0.5.0; analysis.py includes a documented shim restoring the pre-0.13 behavior (verified by an end-to-end synthetic smoke test).

**Results**
- 555 genes pass (padj < 0.05, |log2FC| > 0.5, baseMean > 10): **356 up / 199 down** in combination vs DMSO; top hits include GDA, AHRR, ALDH1A3, CYP1A1.
- 4,545 genes have unavailable padj (independent filtering); preserved as empty CSV cells / JSON `null` — zero coerced zeroes. 0 Cook's outliers replaced.

**Artifacts**
- `output/differential_expression.csv` — 18,029 rows; columns `gene_id, gene_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, significant`
- `output/summary.json` — counts, thresholds, sample lists, significant gene details (strict-valid JSON, no NaN literals)
- `output/analysis.py`, `output/report.md` (345 words, explicitly distinguishes statistical association from causation)

**Validation** — one independent pass re-derived the prefilter count from raw inputs (18,029 ✓), re-applied the pass rule against the CSV (matches `significant` column ✓), checked row count, column schema, unique index, null-vs-zero padj handling, and summary/CSV consistency. All checks passed.