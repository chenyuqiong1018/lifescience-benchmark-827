# Guide Ranking Report — On-Target Activity vs. Annotated Off-Target Risk

Date: 2026-08-27 · Workspace: `ls01-grna-offtarget-rank/C0` · Scope: local inputs only, no external data fetched, inputs unmodified.

## 1. Task and inputs

Rank all 6 candidate guides from `inputs/candidates.csv` by on-target activity and
annotated off-target risk, using the frozen deterministic policy in
`inputs/RANKING_RULE.md`, and emit a decision per guide.

- `inputs/candidates.csv` — 6 guides with `spacer`, `activity_score` (on-target), `coding_bystander_count`.
- `inputs/off_targets.csv` — one annotated off-target per guide with `locus`, `mismatches`, `region`, `expression_risk`.

## 2. Policy applied (from RANKING_RULE.md)

1. **Reject** a guide if any annotated off-target is a coding exon with ≤ 2 mismatches, or has `expression_risk=high`, or if `coding_bystander_count > 0`.
2. **Risk class**: `critical` = high-expression-risk hit, ≥ 2 coding bystanders, or coding-exon hit with exactly 1 mismatch; else `high` = coding-exon hit with ≤ 2 mismatches or exactly 1 coding bystander; else `low`.
3. **Order**: non-rejected before rejected. Non-rejected sorted by `activity_score` desc, then `guide_id`. Rejected sorted with risk `high` before `critical`, then `activity_score` desc, then `guide_id`.
4. **Decision**: rank 1 → `recommend`; other non-rejected → `acceptable`; rejected → `reject`.

Coding/exonic near matches and mismatch counts are treated throughout as **safety evidence**: a coding-exon off-target with few mismatches is presumed capable of being cleaved and, because it lands in expressed coding sequence, of producing a damaging edit.

## 3. Evidence per guide

| guide | activity | off-target locus | region | mismatches | expression risk | coding bystanders | rejection triggers |
|---|---|---|---|---|---|---|---|
| g01 | 0.82 | chr2:100120-100139 | coding_exon | 1 | high | 1 | exon ≤2 mm; high expr; bystanders >0 |
| g02 | 0.67 | chr5:220010-220029 | intron | 3 | low | 0 | none |
| g03 | 0.74 | chr7:330040-330059 | coding_exon | 2 | medium | 0 | exon ≤2 mm |
| g04 | 0.59 | chr9:440070-440089 | intergenic | 4 | low | 0 | none |
| g05 | 0.78 | chr12:550100-550119 | coding_exon | 1 | high | 2 | exon ≤2 mm; high expr; bystanders >0 |
| g06 | 0.64 | chr14:660130-660149 | intergenic | 3 | low | 0 | none |

Risk classes derived per rule 2: g01 `critical`, g02 `low`, g03 `high`, g04 `low`, g05 `critical`, g06 `low`.

## 4. Final ranking (also written to `output/ranked_guides.csv`)

| rank | guide_id | on_target_score | risk_class | decision |
|---|---|---|---|---|
| 1 | g02 | 0.67 | low | recommend |
| 2 | g06 | 0.64 | low | acceptable |
| 3 | g04 | 0.59 | low | acceptable |
| 4 | g03 | 0.74 | high | reject |
| 5 | g01 | 0.82 | critical | reject |
| 6 | g05 | 0.78 | critical | reject |

Ordering check: non-rejected (g02, g06, g04) by activity desc 0.67 > 0.64 > 0.59; rejected with `high` before `critical` (g03 first), then criticals by activity desc (g01 0.82 > g05 0.78).

## 5. Trade-offs, stated explicitly

- **The most active guides are rejected.** g01 (0.82) and g05 (0.78) are the two strongest on-target performers, but each carries a **coding-exon off-target with only 1 mismatch in a high-expression-risk locus**, plus coding bystanders (1 and 2 respectively). Under the frozen rule these are `critical` safety liabilities, so the ranking deliberately sacrifices ≈0.15–0.23 on-target activity to avoid likely deleterious coding edits.
- **g03 is the clearest activity-vs-safety trade-off.** Its activity (0.74) beats every usable guide (best usable = g02 at 0.67), yet a single coding-exon off-target at **2 mismatches** (medium expression risk) is enough to reject it. The rule values exon near-match safety evidence above that activity advantage.
- **The recommended guide is not the panel's most active guide.** g02 (0.67) is only the 4th-highest activity of six, but it is the most active guide with a clean safety profile: its sole off-target is intronic at 3 mismatches with low expression risk and zero coding bystanders. Recommending it accepts weaker on-target activity in exchange for no coding/exonic near-match liability.
- **Mismatch count matters.** All three safe guides have ≥ 3 mismatches and non-coding placement; all three rejected guides have ≤ 2 mismatches inside coding exons. The threshold at 2 mismatches in coding sequence is the decisive separator in this panel.
- **Within the rejected group**, ordering is `high` before `critical` per the frozen rule — this is why g03 (risk `high`, activity 0.74) outranks the more dangerous but higher-activity g01 and g05; among `critical` guides, g01 (0.82) outranks g05 (0.78).

## 6. Reproducibility and validation

- `output/analysis.py` is self-contained and deterministic: it reads only `inputs/candidates.csv` and `inputs/off_targets.csv`, applies the frozen rule, and re-running it reproduces `output/ranked_guides.csv` byte-for-byte (verified by file hash).
- Validated invariants: ranks are unique and sequential 1–6; every input guide appears exactly once; header is exactly `rank,guide_id,on_target_score,risk_class,decision,rationale`; `on_target_score` values match the input `activity_score`; rank 1 is `recommend`; all rejected guides are `reject` with non-empty rationales.

## 7. Scope caveat

Per `RANKING_RULE.md`, this policy is a transparent local benchmark rule, not a claim of universal CRISPR design practice. Real-world guide selection would additionally consider chromatin context, SNP overlap, guide chemistry, and empirical validation.
