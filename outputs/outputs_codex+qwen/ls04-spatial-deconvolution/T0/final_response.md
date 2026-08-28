Done. **Spot_710-1 is a three-way mixture**: Endothelial (~36%), B_Cell (~32%), Macrophage (~32%).

**Deliverables** (all in `output/`, validated):
- `spot_710_composition.csv` — exact columns `cell_type,weight,evidence`; 3 rows; weights 0.359 / 0.321 / 0.320, sum = 1.000, all nonnegative
- `analysis.py` — reproducible pipeline with built-in schema assertion
- `report.md` — methods, evidence, caveats

**How it was determined**
- **Primary method:** NNLS deconvolution of the spot's CPM-normalized profile against per-cell-type mean signatures from the 1,200-cell reference (6 types). Result: Endothelial 0.354, B_Cell 0.317, Macrophage 0.317, T_Cell 0.012, others 0 (fit R² = 0.61).
- **Independent validation:** marker-mass cell-count equivalents (spot counts over each type's top-20 fold-change markers ÷ per-cell reference marker mass): Endothelial 1.27, Macrophage 1.20, B_Cell 1.04 vs ~0.25 background for the absent types — matching NNLS within ~0.05. All 20/20 top markers of the three called types are expressed (mean counts 6–8) while the other types' blocks sit at background (~1.5).

**Key caveat resolved:** cosine nearest-neighbour voting gave 30/30 Endothelial, but ~99% of the spot's 3,693 UMIs come from a shared expression program common to all cells, which dominates full-transcript similarity. The sparse marker blocks — which carry the cell-type signal — clearly show three co-expressed programs, so a single-type call was rejected in favor of the mixture.