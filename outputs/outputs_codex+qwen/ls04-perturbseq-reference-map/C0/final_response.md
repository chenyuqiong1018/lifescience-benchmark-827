Done. All deliverables are in `output/` and validated.

## Results

| Target gene | Query guide | Score (Pearson r) | Runner-up | Confidence |
|---|---|---|---|---|
| **PABPC1** | **guide18** | 0.661 | 0.443 | +0.217 (high) |
| **NUDT21** | **guide13** | 0.711 | 0.608 | +0.103 (high) |
| **LEO1** | **guide14** | 0.335 | 0.419 (PAF1) | −0.084 (flagged ambiguous) |

## Method (single primary pipeline)
Query and reference are separated by a strong cell-state shift (~10.4 pooled SDs on PC1; raw-profile centroid mapping is non-discriminatory, mean best-r ≈ 0.90 for everything). So the pipeline maps **NT-referenced perturbation signatures**: log1p-CPM(10k) group means minus the same-dataset NT-control baseline, scored by Pearson correlation on the top 1,500 reference signal-to-noise (F-ratio) features, with a pre-specified NNLS mixture-model refinement for guides hitting near-collinear protein-complex families (only guide14 was refined: PAF1 → LEO1).

**LEO1 case**: guide14 is unambiguously a PAF1-complex perturbation, and within-complex identity is at the data's resolution limit (gene-specific partial correlations |r| < 0.06). Pearson ranking favors PAF1, but three independent lines of evidence support LEO1: NNLS mixture weights (LEO1 0.25 vs PAF1 0.03), decisive gene-pull margin (guide14 0.335 vs next guide 0.176), and the globally optimal one-to-one assignment placing LEO1 on guide14 at no cost (13.366 vs 13.344). The negative confidence in the CSV explicitly flags this.

**Leakage prevention**: asserted in code — the query's only metadata column is `guide` (anonymized), it carries no target-gene column, no query ID lexically matches any reference symbol, reference labels are used only to build the reference atlas, and no joint embedding/label transfer is used.

**Ambiguity quantification**: permutation p-values (31/32 guides significant at p=0.005; guide32 with 3 cells is non-significant and flagged unreliable), margins/near-ties per guide, low-n flags, attractor analysis (SCAF1 claims 7 weak guides; 15/33 reference genes unreachable by any guide's top hit), subsample stability (81.3%) and Spearman robustness (78.1%).

Deliverables: `output/guide_mapping.csv` (32 rows, exact requested schema), `output/analysis.py` (deterministic, seed 0, ~30 s runtime), `output/report.md`, plus supplementary `output/diagnostics.json`.