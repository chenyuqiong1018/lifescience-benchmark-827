# Frozen local-extension rule

This fixture is a **benchmark-informed local extension**, not an upstream benchmark item.
The table contains confidence estimates and mutation-model hypotheses; it contains no
experimental binding or activity measurements.

Join a mutation to the interval in `confidence.csv` containing its residue number.
Rank the three `region=pocket` candidates first by higher pLDDT, then lower PAE, then lower
absolute predicted ΔΔG, then lexical mutation. Put non-pocket candidates after them.

Use these exact output values:

- `pocket_support = unsupported_low_confidence` when the joined pocket interval has
  `plddt < 50`; otherwise `cautious_support`;
- `confidence_penalty = plddt_lt_50;pae_gt_10A` when both conditions hold (use only the
  applicable semicolon-separated flags if one condition holds, or `none`);
- `decision = defer_structure_validation` for a low-confidence pocket candidate and
  `out_of_scope_non_pocket` for a non-pocket candidate.

`pocket_assessment.json` must contain `pocket_reliable` (boolean),
`prioritize_mutations` (boolean), and `reason`. Here the pocket is reliable only if every
pocket interval has pLDDT >= 70 and PAE <= 10 Å. Consequently, the frozen correct booleans
for this input are both false. Predicted ΔΔG/activity may be reported only as model
hypotheses; it must not be described as measured binding or activity.

Scientific basis: AlphaFold DB guidance classifies pLDDT < 50 as coordinates that should
not be interpreted, and explains that high PAE means relative placement is uncertain:
https://alphafold.ebi.ac.uk/faq . AlphaFold predictions are hypotheses and do not replace
experimental structure determination (Terwilliger et al., 2024,
DOI 10.1038/s41592-023-02087-4).

