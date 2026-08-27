# gRNA Candidate Ranking Report

**Date:** 2026-08-27
**Inputs:** `inputs/candidates.csv` (6 guides), `inputs/off_targets.csv` (6 annotated off-target hits)
**Rule:** `inputs/RANKING_RULE.md` (frozen local benchmark policy)
**Reproducible via:** `output/analysis.py` -> `output/ranked_guides.csv`

## 1. Method

All six candidate guides were scored using only the two supplied CSVs (no external
data fetched, inputs unmodified). The frozen rule was applied exactly:

1. **Rejection:** a guide is rejected if any annotated off-target is a **coding exon
   with <= 2 mismatches**, or has **`expression_risk=high`**, or the guide has
   **`coding_bystander_count > 0`**.
2. **Risk class:** `critical` for a high-expression-risk hit, >= 2 coding bystanders,
   or a coding-exon hit with exactly 1 mismatch; otherwise `high` for a coding-exon
   hit with <= 2 mismatches or exactly 1 coding bystander; otherwise `low`.
3. **Ordering:** non-rejected guides first, by decreasing `activity_score` then
   `guide_id`; rejected guides after, with risk class `high` before `critical`, then
   decreasing `activity_score`, then `guide_id`.
4. **Labels:** rank 1 = `recommend`; other non-rejected = `acceptable`; rejected = `reject`.

Coding/exonic near-matches and mismatch counts are treated as **hard safety evidence**:
they can veto a guide regardless of how high its on-target activity is.

## 2. Final ranking

| Rank | Guide | On-target score | Risk class | Decision   | Key evidence |
|-----:|-------|----------------:|------------|------------|--------------|
| 1    | g02   | 0.67            | low        | recommend  | Intron hit, 3 mismatches, low expression risk, 0 bystanders |
| 2    | g06   | 0.64            | low        | acceptable | Intergenic hit, 3 mismatches, low expression risk, 0 bystanders |
| 3    | g04   | 0.59            | low        | acceptable | Intergenic hit, 4 mismatches, low expression risk, 0 bystanders |
| 4    | g03   | 0.74            | high       | reject     | Coding-exon hit with 2 mismatches (medium expression risk) |
| 5    | g01   | 0.82            | critical   | reject     | Coding-exon hit with 1 mismatch + high expression risk + 1 bystander |
| 6    | g05   | 0.78            | critical   | reject     | Coding-exon hit with 1 mismatch + high expression risk + 2 bystanders |

Ranks are unique (1-6) and every input guide appears exactly once.

## 3. Per-guide rationale

### g02 — recommend (rank 1, score 0.67, risk low)
Single annotated off-target in an **intron** with **3 mismatches** and low expression
risk; zero coding bystanders. It clears every safety gate, so it wins the non-rejected
group on activity.

### g06 — acceptable (rank 2, score 0.64, risk low)
Intergenic off-target with 3 mismatches and low expression risk; zero bystanders.
Safe, but activity 0.64 < g02's 0.67, so it is acceptable rather than recommended.

### g04 — acceptable (rank 3, score 0.59, risk low)
Intergenic off-target with **4 mismatches** (the largest mismatch count in the set)
and low expression risk; zero bystanders. The safest mismatch profile, but the lowest
activity among the non-rejected guides.

### g03 — reject (rank 4, score 0.74, risk high)
The off-target lands in a **coding exon with 2 mismatches** (<= 2 threshold) with
medium expression risk. Under the frozen rule this is disqualifying, so g03 is
rejected despite outscoring every non-rejected guide on activity. Risk class `high`
(coding-exon hit <= 2 mismatches, no critical trigger).

### g01 — reject (rank 5, score 0.82, risk critical)
**Highest on-target activity in the entire set (0.82)**, yet it trips **all three**
safety gates: a coding-exon off-target with only **1 mismatch**, `expression_risk=high`,
and 1 coding bystander. Risk class `critical` (high-expression hit, coding-exon hit
with 1 mismatch). Rejected.

### g05 — reject (rank 6, score 0.78, risk critical)
Also trips **all three** safety gates: coding-exon off-target with **1 mismatch**,
`expression_risk=high`, and **2 coding bystanders**. Risk class `critical`. Within the
rejected group it ranks after g01 because both are `critical` and 0.78 < 0.82.

## 4. Explicit trade-off statement

The recommended guide **g02 (activity 0.67) is only the 4th most active guide** of six.
The three more active guides — g01 (0.82), g05 (0.78), g03 (0.74) — are all rejected
because their off-targets are **coding-exonic near-matches (1-2 mismatches)**, with
high expression risk and/or coding bystanders for g01 and g05. Under the frozen rule,
coding/exonic near-matches and low mismatch counts are treated as hard safety evidence
that **vetoes activity**: we deliberately give up ~0.15 of on-target activity (0.82 ->
0.67) to avoid cleaving sites that are 1-2 mismatches away from coding sequence in
high-expression contexts. This trade-off is intentional and stated here rather than
hidden; a policy that weighted activity above coding-exon safety would rank g01 first.

Secondary trade-off among non-rejected guides: g04 has the safest mismatch profile
(4 mismatches, intergenic) but the lowest activity; the rule prefers activity within
the safe group, so g02 > g06 > g04.

## 5. Notes and limitations

- Each guide has exactly one annotated off-target in `off_targets.csv`; the code
  generalizes to multiple hits per guide (any hit can trigger rejection/risk).
- `coding_bystander_count` is taken as annotated in `inputs/candidates.csv`.
- Risk class is assigned per the frozen rule for all guides, including rejected ones.
- This ranking implements a transparent local benchmark policy; it is not a claim of
  universal CRISPR design practice.
