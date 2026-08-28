Done. All three deliverables are in `output/` and validated.

**Ranking** (per the rule's stable tuple: complete mapping → TM-score → lDDT → RMSD → aligned residues → critical-region error → lexical ID):

| rank | model | global_score | interface_score | critical_residue_risk (Å) | decision |
|---|---|---|---|---|---|
| 1 | model_A | 0.81 | 1 | 4.8 | preferred |
| 2 | model_C | 0.62 | 1 | 4.1 | alternate |
| 3 | model_B | 0.74 | 0 | 8.5 | reject_incomplete_mapping |

**Key drivers**
- **Chain-mapping completeness is decisive:** it's the first sort key, so model_B — despite the best RMSD (2.0 Å), highest lDDT (0.79), and a TM-score (0.74) above model_C's (0.62) — is demoted to rank 3 as `reject_incomplete_mapping` because its mapping is incomplete. Only mapping-complete models can be `preferred`/`alternate`.
- **Critical-region uncertainty (residues 181–240):** the 6th sort key, and no tie required it here. It still qualifies the result: preferred model_A has the largest critical-region error of its own regions (4.8 Å vs 1.4 Å elsewhere), while alternate model_C has the lowest critical-region error (4.1 Å), making it the fallback if that region matters most.

**Deliverables:** `output/model_ranking.csv` (exact 6-column schema, 3 data rows — schema assertion passed), `output/analysis.py` (implements the rule, reproducible), `output/report.md` (full explanation, no coordinate-level/interface/experimental claims beyond the inputs).