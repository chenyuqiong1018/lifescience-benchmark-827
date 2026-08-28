Task complete. **Verdict: the nominated pocket is not reliable (`pocket_reliable = false`, `prioritize_mutations = false`)** — all pocket candidates are deferred to structure validation.

**Deliverables (all in `output/`):**
- `mutation_priorities.csv` — header `rank,mutation,pocket_support,confidence_penalty,decision`, 4 data rows:

| rank | mutation | pocket_support | confidence_penalty | decision |
|---|---|---|---|---|
| 1 | Y218F | unsupported_low_confidence | plddt_lt_50;pae_gt_10A | defer_structure_validation |
| 2 | R244A | unsupported_low_confidence | plddt_lt_50;pae_gt_10A | defer_structure_validation |
| 3 | L265W | unsupported_low_confidence | plddt_lt_50;pae_gt_10A | defer_structure_validation |
| 4 | S301A | cautious_support | none | out_of_scope_non_pocket |

- `pocket_assessment.json` — exactly the three required keys; both booleans `false` with a reason citing the failed thresholds.
- `analysis.py` — reproducible implementation of the frozen rule (stdlib only; join → rank → label → assess).
- `report.md` — full assessment with uncertainty propagation and references.

**Key findings per SCORING_RULE.md:**
- Both pocket intervals fail reliability (need pLDDT ≥ 70 & PAE ≤ 10 Å everywhere): 210–230 has pLDDT 47 / PAE 14.2 Å; 231–270 has pLDDT 43 / PAE 17.8 Å. pLDDT < 50 means coordinates should not be interpreted; PAE > 10 Å means relative placement to the core is uncertain — so pocket geometry cannot support mutation prioritization.
- Ranking: Y218F first on higher pLDDT (47 vs 43); R244A beats L265W on lower |predicted ΔΔG| (2.8 vs 4.1 kcal/mol) after the pLDDT/PAE tie; non-pocket S301A (core, high-confidence interval 271–340) placed last as out of scope.
- Predicted ΔΔG/activity changes are reported strictly as model hypotheses; the fixture contains no measured binding or activity data, and nothing is described as a measured effect.

**Validation:** one production run of `analysis.py` plus a lightweight schema assertion — exact header, 4 rows, ranks 1–4, mutation order, allowed value sets for all three categorical columns, JSON key set, and the frozen `false/false` booleans all passed (`ALL SCHEMA ASSERTIONS PASSED`).