# Guide Ranking Report: On-Target Activity vs. Annotated Off-Target Risk

Date: 2026-08-27
Scope: local workspace only; inputs read from `./inputs` (unchanged); deliverables written to `./output`. No external data was fetched.

## 1. Inputs

| File | Contents |
|---|---|
| `inputs/candidates.csv` | 6 candidate guides with `guide_id`, `spacer`, `activity_score` (on-target activity), `coding_bystander_count` |
| `inputs/off_targets.csv` | One annotated off-target hit per guide with `locus`, `mismatches`, `region`, `expression_risk` |
| `inputs/RANKING_RULE.md` | Frozen, participant-visible deterministic ranking policy (implemented exactly) |

## 2. Method (frozen rule, implemented in `output/analysis.py`)

1. **Rejection**: a guide is rejected if any annotated off-target is a **coding exon with <= 2 mismatches**, or has `expression_risk=high`, or `coding_bystander_count > 0`.
2. **Risk class**:
   - `critical`: any high-expression-risk hit, >= 2 coding bystanders, or a coding-exon hit with exactly 1 mismatch;
   - `high`: otherwise, a coding-exon hit with <= 2 mismatches or exactly 1 coding bystander;
   - `low`: otherwise.
3. **Ordering**: non-rejected guides first (decreasing `activity_score`, then `guide_id`); rejected guides after (risk `high` before `critical`, then decreasing `activity_score`, then `guide_id`).
4. **Decision labels**: rank 1 = `recommend`; other non-rejected = `acceptable`; rejected = `reject`.

Coding/exonic near matches and mismatch counts are treated as **primary safety evidence** throughout: they drive both rejection and risk classification.

## 3. Final ranking (`output/ranked_guides.csv`)

| Rank | Guide | On-target score | Risk class | Decision | Key evidence |
|---|---|---|---|---|---|
| 1 | g02 | 0.67 | low | recommend | intronic hit, 3 mismatches, low expression risk, 0 coding bystanders |
| 2 | g06 | 0.64 | low | acceptable | intergenic hit, 3 mismatches, low expression risk, 0 coding bystanders |
| 3 | g04 | 0.59 | low | acceptable | intergenic hit, 4 mismatches, low expression risk, 0 coding bystanders |
| 4 | g03 | 0.74 | high | reject | coding-exon hit with 2 mismatches (medium expression risk) |
| 5 | g01 | 0.82 | critical | reject | coding-exon hit with 1 mismatch, high expression risk, 1 coding bystander |
| 6 | g05 | 0.78 | critical | reject | coding-exon hit with 1 mismatch, high expression risk, 2 coding bystanders |

All 6 input guides appear exactly once; ranks 1-6 are unique.

## 4. Per-guide assessment

### g02 - rank 1, recommend (activity 0.67, risk low)
Single off-target `chr5:220010-220029` is intronic with 3 mismatches and low expression risk; no coding bystanders. No safety liability under the frozen rule. Selected as the recommendation: best activity among guides with a clean safety profile.

### g06 - rank 2, acceptable (activity 0.64, risk low)
Off-target `chr14:660130-660149` is intergenic with 3 mismatches and low expression risk; no coding bystanders. Safe, slightly less active than g02.

### g04 - rank 3, acceptable (activity 0.59, risk low)
Off-target `chr9:440070-440089` is intergenic with 4 mismatches (the most divergent hit in the set) and low expression risk; no coding bystanders. Safe but lowest activity among kept guides.

### g03 - rank 4, reject (activity 0.74, risk high)
Off-target `chr7:330040-330059` is a **coding exon with 2 mismatches** (medium expression risk). This triggers rejection under the <= 2-mismatch coding-exon rule and sets risk class `high`. **Trade-off stated:** its 0.74 activity is better than every kept guide, but a coding-exon near match with only 2 mismatches is direct safety evidence of a plausible coding-region edit; under the frozen policy that liability outweighs the activity advantage.

### g01 - rank 5, reject (activity 0.82, risk critical)
Off-target `chr2:100120-100139` is a **coding exon with 1 mismatch** and **high expression risk**, plus `coding_bystander_count=1`. All three rejection triggers fire; risk is `critical`. **Trade-off stated:** g01 has the **highest on-target activity in the entire set (0.82)**, so rejecting it is a real efficacy cost - but a 1-mismatch coding-exon hit in a highly expressed locus is the strongest off-target liability class, and the frozen rule does not allow activity to compensate for it.

### g05 - rank 6, reject (activity 0.78, risk critical)
Off-target `chr12:550100-550119` is a **coding exon with 1 mismatch** and **high expression risk**, plus `coding_bystander_count=2` (the maximum in the set). All three rejection triggers fire; risk is `critical`. **Trade-off stated:** 0.78 activity is second-best overall, yet the combination of a 1-mismatch coding-exon hit, high expression risk, and two coding bystanders makes it the least safe guide; it ranks last despite strong activity.

## 5. Explicit trade-off summary

The central trade-off in this dataset is between activity and safety, and it is surfaced rather than hidden:

- The **two most active guides (g01 = 0.82, g05 = 0.78) are rejected** because each carries a 1-mismatch coding-exon off-target with high expression risk (and coding bystanders). Under the frozen rule, coding/exonic near matches and low mismatch counts are decisive safety evidence that activity cannot offset.
- The **recommended guide g02 (0.67) is only mid-range on activity**. Choosing it explicitly sacrifices ~0.15 activity score versus g01 in exchange for eliminating coding-exon near-match and high-expression-risk liabilities.
- Among kept guides, ordering is purely by decreasing on-target activity because all three have `low` risk; among rejected guides, g03 (`high` risk) is ordered ahead of the two `critical` guides per the rule, then by activity.
- Within the rejected group the rule intentionally ranks `high` risk before `critical` risk (a lesser safety liability ranks closer to the acceptable set), which is why g03 (activity 0.74) outranks g01 (0.82) and g05 (0.78) despite lower activity.

## 6. Verification

- [x] Ranks 1-6 unique; every input guide appears exactly once (asserted in `output/analysis.py`).
- [x] Rejection/risk logic matches `inputs/RANKING_RULE.md` verbatim.
- [x] Coding/exonic near matches and mismatch counts used as explicit safety evidence in both decisions and rationales.
- [x] Activity-vs-safety trade-offs stated explicitly for g01, g03, g05 and for the recommendation of g02.
- [x] Inputs unmodified; no external data fetched; all deliverables in `./output` (`ranked_guides.csv`, `analysis.py`, `report.md`).

Caveat: the frozen rule is a transparent local benchmark policy, not a claim of universal CRISPR design practice.
