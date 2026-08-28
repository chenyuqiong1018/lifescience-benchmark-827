# Pocket Reliability Assessment and Mutation Prioritization

**Verdict: the nominated pocket is NOT reliable. `pocket_reliable = false`, `prioritize_mutations = false`.**
All three pocket-residue mutation candidates are deferred to structure validation; none should be
prioritized for experimental follow-up on the basis of this predicted structure alone.

## 1. Inputs and method

Only files under `inputs/` were used:

- `confidence.csv` — per-residue-interval confidence estimates (pLDDT, PAE-to-core) with pocket membership flags.
- `mutation_candidates.csv` — mutation hypotheses with predicted ΔΔG and predicted activity change (model output only).
- `SCORING_RULE.md` — the frozen local-extension rule, implemented exactly in `output/analysis.py`.

Method (per the frozen rule):

1. Each mutation was joined to the interval in `confidence.csv` containing its residue number.
2. The three `region=pocket` candidates were ranked by higher pLDDT, then lower PAE, then lower
   absolute predicted ΔΔG, then lexical mutation name; the non-pocket candidate follows them.
3. `pocket_support`, `confidence_penalty`, and `decision` use the exact output values prescribed
   by the rule (see `output/mutation_priorities.csv`).
4. The pocket is deemed reliable only if **every** pocket interval has pLDDT >= 70 and PAE <= 10 Å.

## 2. Confidence data and joins

| mutation | residue | joined interval | pLDDT | PAE to core (Å) | interval pocket member |
|----------|---------|-----------------|-------|------------------|------------------------|
| Y218F    | 218     | 210–230         | 47    | 14.2             | true                   |
| R244A    | 244     | 231–270         | 43    | 17.8             | true                   |
| L265W    | 265     | 231–270         | 43    | 17.8             | true                   |
| S301A    | 301     | 271–340         | 84    | 3.4              | false (core)           |

Both pocket intervals (210–230 and 231–270) fail the reliability requirement
(pLDDT 47 / 43 < 70 and PAE 14.2 / 17.8 Å > 10 Å), so `pocket_reliable = false`.

## 3. Ranking and decisions (`output/mutation_priorities.csv`)

| rank | mutation | pocket_support | confidence_penalty | decision |
|------|----------|----------------|--------------------|----------|
| 1 | Y218F | unsupported_low_confidence | plddt_lt_50;pae_gt_10A | defer_structure_validation |
| 2 | R244A | unsupported_low_confidence | plddt_lt_50;pae_gt_10A | defer_structure_validation |
| 3 | L265W | unsupported_low_confidence | plddt_lt_50;pae_gt_10A | defer_structure_validation |
| 4 | S301A | cautious_support | none | out_of_scope_non_pocket |

Ranking rationale: Y218F ranks first on higher joined-interval pLDDT (47 vs 43). R244A and L265W tie
on pLDDT (43) and PAE (17.8 Å), so the tie breaks on lower absolute predicted ΔΔG (2.8 vs 4.1 kcal/mol),
placing R244A second. S301A is a `core` (non-pocket) candidate and is placed after all pocket candidates.
Its joined interval (271–340, pLDDT 84, PAE 3.4 Å) is high-confidence, so under the frozen rule its
`pocket_support` label is the `cautious_support` fallback and no penalty applies, but as a non-pocket
candidate it is out of scope for pocket-based prioritization (`out_of_scope_non_pocket`).

## 4. Propagation of pLDDT/PAE uncertainty

- **pLDDT < 50 (pocket intervals 47 and 43):** AlphaFold DB guidance classifies pLDDT < 50 as
  coordinates that should not be interpreted. Side-chain placements in the pocket are therefore not
  trustworthy, and any geometry derived from them (distances, contacts, steric complementarity) is
  unreliable.
- **PAE > 10 Å (14.2 and 17.8 Å):** high PAE means the relative placement of the pocket with respect
  to the protein core is uncertain. Even locally plausible features cannot be positioned reliably
  relative to the rest of the structure.
- Because every pocket candidate maps to a low-confidence interval, the penalty
  `plddt_lt_50;pae_gt_10A` and the decision `defer_structure_validation` are applied to all three.
  The ranking among them is only a relative ordering under the frozen rule; it does **not** imply any
  of them is ready for prioritization (`prioritize_mutations = false`).
- The well-confidence core interval (271–340) does not rescue the pocket: pocket reliability requires
  **every** pocket interval to pass the thresholds, and neither does.

## 5. Predicted effects are model hypotheses, not measurements

`predicted_ddg_kcal_mol` and `predicted_activity_change` in `mutation_candidates.csv` are
computational model hypotheses. This fixture contains **no experimental binding or activity
measurements**, and nothing in the outputs or this report treats predicted ΔΔG or predicted activity
changes as measured effects. AlphaFold-style predictions are hypotheses and do not replace experimental
structure determination (Terwilliger et al., 2024).

## 6. Conclusion and recommended next step

The pocket fails every reliability threshold (pLDDT 47/43 vs required >= 70; PAE 14.2/17.8 Å vs
required <= 10 Å). Mutation prioritization is therefore **deferred**: obtain experimental structural
validation (e.g., cryo-EM/X-ray of the pocket region, or orthogonal evidence) before investing in
mutagenesis based on this pocket. The core candidate S301A is out of scope for pocket-based
prioritization.

## References

- AlphaFold Protein Structure Database FAQ — pLDDT/PAE interpretation: https://alphafold.ebi.ac.uk/faq
- Terwilliger, D. W. et al. (2024). AlphaFold predictions are valuable hypotheses and do not replace
  experimental structure determination. *Nature Methods*. DOI: 10.1038/s41592-023-02087-4

## Reproducibility

`output/analysis.py` (Python 3, standard library only) reads `inputs/confidence.csv` and
`inputs/mutation_candidates.csv`, applies the frozen rule from `inputs/SCORING_RULE.md`, and writes
`output/mutation_priorities.csv` and `output/pocket_assessment.json`. Run with:

```
python output/analysis.py
```