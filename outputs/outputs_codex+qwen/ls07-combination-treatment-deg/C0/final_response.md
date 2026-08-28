Done. All four deliverables are in `output/`:

**Results** — Cisplatin_IC50_CBD_IC50 (samples 9-1/9-2/9-3) vs DMSO (3-1/3-2/3-3), combination as numerator:
- 63,677 genes → **18,029 retained** by the pre-fit rule (raw count > 10 in ≥1 of the 6 samples)
- PyDESeq2 0.5.0, design `~condition` only, `refit_cooks=True`, contrast `["condition", "Cisplatin_IC50_CBD_IC50", "DMSO"]`
- **557 genes pass** (padj < 0.05, |log2FC| > 0.5, baseMean > 10): 359 up / 198 down; top hits include CDKN1A, CYP1A1, AHRR, ALDH1A3, GDA; strongest down: ANKRD37
- 4,196 genes lack padj due to DESeq2 independent filtering — preserved as null (empty CSV cells / JSON `null`), never zeroed

**Environment fix**: the installed anndata 0.13.3 reshapes 1-D `varm`/`obsm` arrays to (n,1), which crashes PyDESeq2 0.5.0. `analysis.py` contains a minimal, documented shim restoring 1-D round-tripping; the library itself is unmodified.

**Validation** (independent, lightweight):
- BH re-adjustment of the 13,833 non-null p-values reproduces `padj` to 3×10⁻¹⁵ max deviation
- Significance rule re-derived from columns matches `summary.json` (557 = 359 + 198)
- Contrast direction confirmed from raw counts (GDA: DMSO mean 2,048 vs combination 4,333 → positive LFC)
- Prefilter count independently verified; final assertion confirmed exact filenames, 18,029-row schema, unique Ensembl keys, no `NaN` strings

`report.md` is 382 words and explicitly distinguishes statistical association from causation (downstream stress/cell-state effects, batch confounding, need for perturbation experiments).