# Low-confidence pocket assessment

## Overall assessment

The nominated pocket is **not reliable enough to prioritize mutations**. This is a confidence-limited structural hypothesis, not evidence that any mutation changes binding, stability, or activity.

## Evidence audit

**Strengths.** Both input schemas are complete, the four confidence intervals are ordered and non-overlapping, every mutation maps to exactly one interval, and each `region` label agrees with the joined `pocket_member` flag. The ranking rule is prespecified and deterministic.

**Critical concern.** Both pocket intervals fail both confidence criteria. Residues 210–230 have pLDDT 47 and PAE 14.2 Å; residues 231–270 have pLDDT 43 and PAE 17.8 Å. Low local confidence weakens coordinate interpretation, while high PAE weakens confidence in placement relative to the core. The combined limitation is severe enough to defer every pocket candidate.

## Ranked hypotheses

| Rank | Mutation | Region | Joined interval | pLDDT | PAE to core (Å) | Predicted ΔΔG hypothesis (kcal/mol) | Activity hypothesis | Decision |
|---:|---|---|---|---:|---:|---:|---|---|
| 1 | Y218F | pocket | 210–230 | 47 | 14.2 | 0.4 | increase | defer_structure_validation |
| 2 | R244A | pocket | 231–270 | 43 | 17.8 | 2.8 | increase | defer_structure_validation |
| 3 | L265W | pocket | 231–270 | 43 | 17.8 | 4.1 | unknown | defer_structure_validation |
| 4 | S301A | core | 271–340 | 84 | 3.4 | 0.2 | no_change | out_of_scope_non_pocket |

Y218F ranks first because its joined interval has higher pLDDT and lower PAE than the interval for R244A and L265W. R244A and L265W tie on confidence, so the lower absolute predicted ΔΔG hypothesis places R244A first. S301A is non-pocket and is placed last by rule.

## Claim boundary and recommendation

The ΔΔG and activity columns are unvalidated model hypotheses. They are not measured effects and cannot rescue an uninterpretable pocket. The proportionate next step is structure validation or a better-supported structural model before using these hypotheses to prioritize mutagenesis; experimental measurements would then be needed to establish actual stability, binding, or activity effects.
