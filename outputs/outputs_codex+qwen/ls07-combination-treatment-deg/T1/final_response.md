Done. All four deliverables are in `output/`, validated against independent recomputation.

## What was run
- **Samples**: only `Cisplatin_IC50_CBD_IC50` (9-1/9-2/9-3, numerator) vs `DMSO` (3-1/3-2/3-3, denominator) — no other groups touched.
- **Pre-fit filter**: genes with raw count > 10 in ≥1 of the six samples → 18,029 of 63,677 genes retained.
- **Model**: PyDESeq2 **0.5.0**, design `~Group` (only term), `refit_cooks=True`, standard `DeseqStats` contrast `["Group", "Cisplatin_IC50_CBD_IC50", "DMSO"]`. No Cook's outliers needed replacement.
- **Pass rule**: `padj < 0.05` & `|log2FC| > 0.5` & `baseMean > 10` → **557 DE genes** (359 up, 198 down); top hits GDA, AHRR, ALDH1A3, CYP1A1, CDKN1A.

## Deliverables
| File | Content |
|---|---|
| `output/differential_expression.csv` | 18,029 rows keyed by Ensembl ID: gene_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, passed. The 4,196 independent-filtering nulls are empty cells, never zero. |
| `output/summary.json` | Contrast, samples, thresholds, gene counts (tested 18,029 / padj-null 4,196 / passed 557 = 359 up + 198 down), software versions. |
| `output/analysis.py` | Reproducible script (documented anndata compatibility shim included). |
| `output/report.md` | 423 words; explicitly separates statistical association from causation (joint dosing, 3-replicate limits, batch/viability confounds). |

## Validation (independent path)
- Prefilter recount from raw matrix: 18,029 ✓
- BH recomputed independently over the 13,833-gene family matched pydeseq2 `padj` to 3×10⁻¹⁵; the null cutoff cleanly separates at baseMean ≈ 24.89 ✓
- `passed` flag recomputed from the rule matches exactly; no padj cell is `"0"`; sample IDs match the layout ✓

**One environment note**: the installed anndata 0.13 reshapes 1-D `varm` arrays to 2-D, which breaks PyDESeq2 0.5.0 (`IndexError` in `fit_genewise_dispersions`). `analysis.py` applies a minimal, documented monkeypatch restoring pre-0.11 storage semantics; PyDESeq2 itself was used unmodified at exactly version 0.5.0.