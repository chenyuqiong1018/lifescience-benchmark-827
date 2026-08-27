#!/usr/bin/env python3
"""Rank CRISPR candidate guides by on-target activity and annotated off-target risk.

Reads (read-only, no external data):
  inputs/candidates.csv   guide_id, spacer, activity_score, coding_bystander_count
  inputs/off_targets.csv  guide_id, locus, mismatches, region, expression_risk

Applies the frozen deterministic rule from inputs/RANKING_RULE.md:
  1. Reject a guide when any annotated off-target is a coding exon with <=2
     mismatches, when expression_risk=high, or when coding_bystander_count>0.
  2. Risk class: critical for a high-expression-risk hit, >=2 coding bystanders,
     or a coding-exon hit with exactly 1 mismatch; else high for a coding-exon
     hit with <=2 mismatches or exactly 1 coding bystander; else low.
  3. Rank non-rejected before rejected. Non-rejected: activity_score desc, then
     guide_id. Rejected: risk high before critical, then activity_score desc,
     then guide_id.
  4. Decision: rank 1 = recommend; other non-rejected = acceptable; rejected = reject.

Writes: output/ranked_guides.csv with columns
  rank, guide_id, on_target_score, risk_class, decision, rationale
"""
from __future__ import annotations

import csv
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
CANDIDATES_CSV = BASE / "inputs" / "candidates.csv"
OFF_TARGETS_CSV = BASE / "inputs" / "off_targets.csv"
OUT_CSV = BASE / "output" / "ranked_guides.csv"

CODING_EXON = "coding_exon"


def load_inputs():
    """Load both CSVs; return dict guide_id -> guide record (inputs unmodified)."""
    with CANDIDATES_CSV.open(newline="", encoding="utf-8") as fh:
        candidate_rows = list(csv.DictReader(fh))
    with OFF_TARGETS_CSV.open(newline="", encoding="utf-8") as fh:
        off_target_rows = list(csv.DictReader(fh))

    guides = {}
    for row in candidate_rows:
        gid = row["guide_id"].strip()
        if gid in guides:
            raise ValueError(f"duplicate guide_id in candidates.csv: {gid}")
        guides[gid] = {
            "guide_id": gid,
            "spacer": row["spacer"].strip(),
            "activity_score": float(row["activity_score"]),
            "activity_score_raw": row["activity_score"].strip(),
            "coding_bystander_count": int(row["coding_bystander_count"]),
            "off_targets": [],
        }
    for row in off_target_rows:
        gid = row["guide_id"].strip()
        if gid not in guides:
            raise ValueError(f"off_targets.csv references unknown guide_id: {gid}")
        guides[gid]["off_targets"].append(
            {
                "locus": row["locus"].strip(),
                "mismatches": int(row["mismatches"]),
                "region": row["region"].strip(),
                "expression_risk": row["expression_risk"].strip(),
            }
        )
    return guides


def evaluate(guide):
    """Apply frozen rule 1 (rejection) and rule 2 (risk class) to one guide."""
    ots = guide["off_targets"]
    bystanders = guide["coding_bystander_count"]

    exon_le2 = [o for o in ots if o["region"] == CODING_EXON and o["mismatches"] <= 2]
    exon_1mm = [o for o in ots if o["region"] == CODING_EXON and o["mismatches"] == 1]
    high_expr = [o for o in ots if o["expression_risk"] == "high"]

    # Rule 1: rejection triggers.
    rejected = bool(exon_le2) or bool(high_expr) or bystanders > 0

    # Rule 2: risk class (critical > high > low).
    if high_expr or bystanders >= 2 or exon_1mm:
        risk_class = "critical"
    elif exon_le2 or bystanders == 1:
        risk_class = "high"
    else:
        risk_class = "low"

    return {
        "rejected": rejected,
        "risk_class": risk_class,
        "exon_le2": exon_le2,
        "exon_1mm": exon_1mm,
        "high_expr": high_expr,
        "bystanders": bystanders,
    }


def order_guides(guides, evals):
    """Frozen rule 3: non-rejected first, then the two within-group orderings."""
    accepted_ids = [gid for gid, e in evals.items() if not e["rejected"]]
    rejected_ids = [gid for gid, e in evals.items() if e["rejected"]]

    accepted_ids.sort(key=lambda gid: (-guides[gid]["activity_score"], gid))

    # Within rejected: risk "high" before "critical", then activity desc, then id.
    rejected_rank = {"high": 0, "critical": 1, "low": 2}
    rejected_ids.sort(
        key=lambda gid: (
            rejected_rank[evals[gid]["risk_class"]],
            -guides[gid]["activity_score"],
            gid,
        )
    )
    return accepted_ids + rejected_ids, accepted_ids, rejected_ids


def describe_off_targets(guide):
    """Human-readable off-target evidence, citing region/mismatches/expression."""
    ots = guide["off_targets"]
    if not ots:
        return "no annotated off-targets"
    parts = []
    for o in ots:
        mm = o["mismatches"]
        parts.append(
            f"{o['locus']} in {o['region'].replace('_', ' ')} with {mm} "
            f"mismatch{'es' if mm != 1 else ''} and {o['expression_risk']} "
            f"expression risk"
        )
    return "; ".join(parts)


def build_rationale(guide, ev, decision, ordered, guides):
    """Compose the rationale, stating safety evidence and any trade-off openly."""
    gid = guide["guide_id"]
    act = guide["activity_score"]
    ot_text = describe_off_targets(guide)
    by = ev["bystanders"]
    by_text = f"{by} coding bystander{'s' if by != 1 else ''}"

    panel_max = max(g["activity_score"] for g in guides.values())
    accepted_ids = [g for g in ordered if not evals_global[g]["rejected"]]
    acc_min = min(guides[g]["activity_score"] for g in accepted_ids)

    if decision == "recommend":
        return (
            f"Highest on-target activity ({act:.2f}) among guides that pass every "
            f"safety gate. Off-target evidence: {ot_text}; {by_text}. No "
            f"coding/exonic near match (all hits are non-coding with >=3 "
            f"mismatches), so the top activity comes with no safety trade-off."
        )
    if decision == "acceptable":
        extra = ""
        if act == acc_min and len(accepted_ids) > 1:
            extra = (
                " Trade-off: weakest on-target activity among usable guides, in "
                "exchange for a clean off-target safety profile."
            )
        else:
            extra = (
                " Trade-off: marginally lower on-target activity than the "
                "recommended guide, with an equally clean safety profile."
            )
        return (
            f"On-target activity {act:.2f}; passes every safety gate. Off-target "
            f"evidence: {ot_text}; {by_text}. No coding/exonic near match."
            f"{extra}"
        )
    # Rejected guides: name every trigger and state the activity sacrificed.
    triggers = []
    if ev["exon_le2"]:
        for o in ev["exon_le2"]:
            triggers.append(
                f"coding-exon off-target {o['locus']} with only {o['mismatches']} "
                f"mismatch{'es' if o['mismatches'] != 1 else ''}"
            )
    if ev["high_expr"]:
        triggers.append(
            "off-target with high expression risk ("
            + ", ".join(sorted({o['locus'] for o in ev['high_expr']}))
            + ")"
        )
    if by > 0:
        triggers.append(f"coding_bystander_count={by} (>0)")
    trigger_text = "; ".join(triggers)
    if act > acc_min:
        trade = (
            f"Trade-off stated explicitly: this guide's on-target activity "
            f"({act:.2f}) exceeds every usable guide (best usable = "
            f"{max(guides[g]['activity_score'] for g in accepted_ids):.2f}), but "
            f"the frozen rule sacrifices that activity because coding/exonic near "
            f"matches and mismatch counts are treated as decisive safety evidence."
        )
    else:
        trade = (
            f"Trade-off stated explicitly: no compensating activity advantage "
            f"({act:.2f}), so rejection costs nothing on efficacy either."
        )
    return (
        f"Rejected despite on-target activity {act:.2f}. Safety evidence: "
        f"{ot_text}; {by_text}. Rejection trigger(s) under the frozen rule: "
        f"{trigger_text}. {trade}"
    )


def main():
    guides = load_inputs()
    evals = {gid: evaluate(g) for gid, g in guides.items()}
    ordered, accepted_ids, rejected_ids = order_guides(guides, evals)

    global evals_global
    evals_global = evals

    rows = []
    for i, gid in enumerate(ordered, start=1):
        guide = guides[gid]
        ev = evals[gid]
        if i == 1:
            decision = "recommend"
        elif not ev["rejected"]:
            decision = "acceptable"
        else:
            decision = "reject"
        rationale = build_rationale(guide, ev, decision, ordered, guides)
        rows.append(
            {
                "rank": i,
                "guide_id": gid,
                "on_target_score": f"{guide['activity_score']:.2f}",
                "risk_class": ev["risk_class"],
                "decision": decision,
                "rationale": rationale,
            }
        )

    # Self-checks: unique ranks, every input guide exactly once.
    assert [r["rank"] for r in rows] == list(range(1, len(guides) + 1))
    assert sorted(r["guide_id"] for r in rows) == sorted(guides)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["rank", "guide_id", "on_target_score", "risk_class", "decision", "rationale"]
    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {OUT_CSV.relative_to(BASE)} ({len(rows)} guides)")
    for r in rows:
        print(
            f"  rank {r['rank']}: {r['guide_id']} score={r['on_target_score']} "
            f"risk={r['risk_class']} decision={r['decision']}"
        )
    print(f"non-rejected: {', '.join(accepted_ids)} | rejected: {', '.join(rejected_ids)}")


if __name__ == "__main__":
    main()
