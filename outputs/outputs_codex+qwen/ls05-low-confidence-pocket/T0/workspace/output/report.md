# Pocket Reliability Assessment and Mutation Prioritization

## Question
Is the nominated pocket reliable enough to prioritize mutations?

**Answer: No.** `pocket_reliable = false`, `prioritize_mutations = false`. All three
pocket mutation candidates are deferred to structure validation; no mutation should be
prioritized for experimental follow-up on the current structural evidence.

## Inputs used (only files under `inputs/`)
- `inputs/confidence.csv` - per-interval pLDDT and PAE-to-core with pocket membership.
- `inputs/mutation_candidates.csv` - candidate mutations with predicted ddG and activity
  change (model hypotheses only).
- `inputs/SCORING_RULE.md` - frozen rule applied exactly.

## Confidence evidence and uncertainty propagation
| Interval | pLDDT | PAE to core (A) | Pocket member |
|---|---|---|---|
| 1-209 | 88 | 2.1 | false |
| 210-230 | 47 | 14.2 | true |
| 231-270 | 43 | 17.8 | true |
| 271-340 | 84 | 3.4 | false |

Reliability requires **every** pocket interval to have pLDDT >= 70 and PAE <= 10 A.
Both pocket intervals fail on both criteria:

- **210-230:** pLDDT 47 (< 50) and PAE 14.2 A (> 10 A).
- **231-270:** pLDDT 43 (< 50) and PAE 17.8 A (> 10 A).

Per AlphaFold DB guidance, pLDDT < 50 marks coordinates that **should not be
interpreted**, and high PAE means the relative placement of regions is uncertain
(https://alphafold.ebi.ac.uk/faq). AlphaFold predictions are hypotheses and do not
replace experimental structure determination (Terwilliger et al., 2024,
DOI 10.1038/s41592-023-02087-4).

Because the pocket coordinates themselves are unreliable, this uncertainty propagates to
every structure-based quantity derived from them: the predicted ddG and activity-change
values are **model hypotheses conditional on a low-confidence structure**, not measured
binding or activity effects. They are reported below only as such, and they cannot
justify experimental prioritization until the structure is validated.

## Join of mutations to confidence intervals
| Mutation | Residue | Joined interval | pLDDT | PAE (A) | Pocket member |
|---|---|---|---|---|---|
| Y218F | 218 | 210-230 | 47 | 14.2 | true |
| R244A | 244 | 231-270 | 43 | 17.8 | true |
| L265W | 265 | 231-270 | 43 | 17.8 | true |
| S301A | 301 | 271-340 | 84 | 3.4 | false |

## Ranking (frozen rule)
Pocket candidates ranked by higher pLDDT, then lower PAE, then lower absolute predicted
ddG (hypothesis values), then lexical mutation; non-pocket candidates placed after:

1. **Y218F** - pLDDT 47 (highest among pocket candidates).
2. **R244A** - tied with L265W on pLDDT 43 and PAE 17.8 A; lower |predicted ddG| (2.8
   vs 4.1 kcal/mol, model hypotheses).
3. **L265W** - same interval as R244A; higher |predicted ddG| hypothesis value.
4. **S301A** - non-pocket (core), placed after pocket candidates.

## Decisions written to `output/mutation_priorities.csv`
- All three pocket candidates join intervals with pLDDT < 50, so
  `pocket_support = unsupported_low_confidence`,
  `confidence_penalty = plddt_lt_50;pae_gt_10A` (both flags hold), and
  `decision = defer_structure_validation`.
- S301A joins the well-supported core interval (pLDDT 84, PAE 3.4 A), so
  `pocket_support = cautious_support`, `confidence_penalty = none`, and
  `decision = out_of_scope_non_pocket`.

## Assessment (`output/pocket_assessment.json`)
- `pocket_reliable: false` - at least one pocket interval (here, both) violates
  pLDDT >= 70 and PAE <= 10 A.
- `prioritize_mutations: false` - mutation prioritization is withheld until the pocket
  structure is validated experimentally or replaced by a higher-confidence model.

## Recommended next step (deferred, not executed here)
Obtain or commission experimental structural validation of the pocket region (or a
higher-confidence model/assembly) before revisiting any mutation prioritization. Until
then, predicted ddG/activity values remain unvalidated hypotheses.

## Reproducibility
`output/analysis.py` reads only `inputs/confidence.csv` and
`inputs/mutation_candidates.csv`, applies the frozen rule from `inputs/SCORING_RULE.md`,
and regenerates `output/mutation_priorities.csv` and `output/pocket_assessment.json`.
