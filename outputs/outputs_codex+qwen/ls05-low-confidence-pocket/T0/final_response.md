All deliverables are complete and validated in `output/`:

**Assessment conclusion:** The nominated pocket is **not reliable** — both pocket intervals fail the thresholds (210–230: pLDDT 47, PAE 14.2 Å; 231–270: pLDDT 43, PAE 17.8 Å vs. required pLDDT ≥ 70 and PAE ≤ 10 Å). Both fall below pLDDT 50, which AlphaFold guidance says should not be interpreted, so `pocket_reliable = false` and `prioritize_mutations = false`.

**`mutation_priorities.csv`** (frozen ranking: higher pLDDT → lower PAE → lower |predicted ΔΔG| → lexical; non-pocket last):
| rank | mutation | pocket_support | confidence_penalty | decision |
|---|---|---|---|---|
| 1 | Y218F | unsupported_low_confidence | plddt_lt_50;pae_gt_10A | defer_structure_validation |
| 2 | R244A | unsupported_low_confidence | plddt_lt_50;pae_gt_10A | defer_structure_validation |
| 3 | L265W | unsupported_low_confidence | plddt_lt_50;pae_gt_10A | defer_structure_validation |
| 4 | S301A | cautious_support | none | out_of_scope_non_pocket |

**Other artifacts:**
- `pocket_assessment.json` — exactly the three required keys; both booleans `false`, with a reason citing the failing intervals and the pLDDT/PAE guidance.
- `analysis.py` — reads only `inputs/`, applies SCORING_RULE.md exactly, regenerates the CSV/JSON.
- `report.md` — propagates the pLDDT/PAE uncertainty into the decisions and states explicitly that predicted ΔΔG/activity values are model hypotheses conditional on a low-confidence structure, not measured binding or activity effects.

Schema assertion passed: exact header, 4 rows, rank order, allowed values, and JSON keys/types all verified.