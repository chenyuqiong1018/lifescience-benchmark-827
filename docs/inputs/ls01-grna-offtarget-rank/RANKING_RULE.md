# Frozen guide-ranking rule

This local fixture uses the following participant-visible deterministic rule.

1. A guide is **rejected** when any annotated off-target is a coding exon with at most two mismatches, when `expression_risk=high`, or when `coding_bystander_count>0`.
2. Risk class is `critical` for a high-expression-risk hit, at least two coding bystanders, or a coding-exon hit with one mismatch; otherwise it is `high` for a coding-exon hit with at most two mismatches or one coding bystander, and `low` otherwise.
3. Rank non-rejected guides before rejected guides. Within the non-rejected group rank by decreasing `activity_score`, then `guide_id`. Within rejected guides rank `high` before `critical`, then by decreasing `activity_score`, then `guide_id`.
4. Label rank 1 `recommend`, other non-rejected guides `acceptable`, and rejected guides `reject`.

The rule is a transparent local benchmark policy, not a claim of universal CRISPR design practice.
