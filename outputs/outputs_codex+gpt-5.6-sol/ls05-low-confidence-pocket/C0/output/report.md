# Low-confidence pocket assessment

## Conclusion

The nominated pocket is **not reliable enough to prioritize mutations**. Both pocket intervals fail the required pLDDT ≥ 70 and PAE ≤ 10 Å criteria. Every pocket candidate therefore receives `unsupported_low_confidence`, both confidence penalties, and `defer_structure_validation`.

## Ranked candidates

| Rank | Mutation | Region | Joined pLDDT | PAE to core (Å) | Predicted ΔΔG hypothesis (kcal/mol) | Predicted activity hypothesis | Decision |
|---:|---|---|---:|---:|---:|---|---|
| 1 | Y218F | pocket | 47 | 14.2 | 0.4 | increase | defer_structure_validation |
| 2 | R244A | pocket | 43 | 17.8 | 2.8 | increase | defer_structure_validation |
| 3 | L265W | pocket | 43 | 17.8 | 4.1 | unknown | defer_structure_validation |
| 4 | S301A | core | 84 | 3.4 | 0.2 | no_change | out_of_scope_non_pocket |

Pocket candidates are ordered exactly by higher joined pLDDT, then lower PAE, then lower absolute predicted ΔΔG, then mutation name. Y218F ranks first because its interval has pLDDT 47 and PAE 14.2 Å, compared with pLDDT 43 and PAE 17.8 Å for R244A and L265W. R244A precedes L265W because their confidence values tie and |2.8| < |4.1| kcal/mol. S301A is non-pocket and is placed after all pocket candidates.

## Uncertainty propagation

pLDDT below 50 means the pocket coordinates should not be interpreted, while PAE above 10 Å indicates uncertain placement relative to the core. The simultaneous failures apply to all three pocket candidates and prevent mutation prioritization from this structure hypothesis. The supplied ΔΔG and activity fields are model hypotheses used only as prescribed tie-breaker/context; they are not measured binding, stability, or activity effects.
