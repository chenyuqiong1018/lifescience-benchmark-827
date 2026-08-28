All deliverables are complete and validated. Summary:

**Ranking (per the stable tuple in `inputs/SCORING_RULE.md`):**

| rank | model_id | global_score | interface_score | critical_residue_risk (Å) | decision |
|---|---|---|---|---|---|
| 1 | model_A | 0.81 | 1 | 4.8 | preferred |
| 2 | model_C | 0.62 | 1 | 4.1 | alternate |
| 3 | model_B | 0.74 | 0 | 8.5 | reject_incomplete_mapping |

**Key reasoning:**
- **Chain-mapping completeness is the dominant criterion**, so model_B ranks last and is rejected despite having the best lddt (0.79) and rmsd_a (2.0) — an incomplete mapping makes its other metrics insufficient for selection. Completeness is exposed as `interface_score` (1/0) per the rule.
- Among complete models, higher `tm_score` decides: model_A (0.81) → `preferred`, model_C (0.62) → `alternate`.
- **Critical-region uncertainty** (residues 181–240, length-weighted mean error) is criterion 6: a tie-breaker applied only after completeness, tm_score, lddt, rmsd_a, and coverage. No ties occurred in this dataset, so it did not alter the order, but it corroborates the decisions (model_B is worst at 8.5 Å; model_C's low 4.1 Å correctly does not override its weaker global similarity).

**Artifacts in `output/`:**
- `model_ranking.csv` — exact required header `rank,model_id,global_score,interface_score,critical_residue_risk,decision`, exactly 3 data rows, no extras.
- `analysis.py` — primary method; reads only `inputs/`, implements the full 7-key tuple and length-weighted risk computation.
- `report.md` — method, derivation, completeness/uncertainty discussion, and caveats (no coordinate-level, interface, or experimental claims; `global_score` is exactly `tm_score`, no invented composite).

Validation: one independent recomputation pass asserted schema exactness, row count, allowed decision values, value cross-checks against both input files, and re-derived the ordering from the rule tuple — all passed.