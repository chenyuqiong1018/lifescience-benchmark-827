All deliverables are complete and validated.

**Ranking result** (`output/model_ranking.csv`, schema-asserted: exact header, 3 rows, allowed decision values):

| rank | model_id | global_score | interface_score | critical_residue_risk | decision |
|---|---|---|---|---|---|
| 1 | model_A | 0.81 | 1 | 4.8 | preferred |
| 2 | model_C | 0.62 | 1 | 4.1 | alternate |
| 3 | model_B | 0.74 | 0 | 8.5 | reject_incomplete_mapping |

**Method (single primary run):** `output/analysis.py` reads only `inputs/model_metrics.csv` and `inputs/residue_errors.csv`, computes `critical_residue_risk` as the length-weighted mean of `mean_error_a` over rows overlapping residues 181–240 (each model's single row spans the whole interval, so risk equals that row's value), and sorts by the rule's stable tuple: completeness → higher TM-score → higher lDDT → lower RMSD → more aligned residues → lower critical risk → lexical model_id. `global_score = tm_score` (no composite), `interface_score` is the completeness indicator per the rule.

**Key findings explained in `output/report.md`:**
- **Chain-mapping completeness** is the dominant key: model_B has the best RMSD (2.0 Å), highest lDDT (0.79), and a higher TM-score than model_C, yet ranks last and is rejected solely because its mapping is incomplete. Among complete models, TM-score alone separates A (0.81) from C (0.62).
- **Critical-region risk** is tie-break key 6, so it does not change this dataset's ordering — model_C's lower risk (4.1 vs 4.8 Å) cannot overcome TM-score; model_B's high risk (8.5 Å) compounds but does not cause its rejection. Risk uncertainty is nil here since every critical row covers 181–240 exactly; the report notes how length-weighting would matter for partial overlaps.

No coordinate-level, interface, or experimental claims are made beyond the supplied metrics.