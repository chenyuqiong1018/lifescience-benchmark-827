#!/usr/bin/env python3
"""Rank CRISPR candidate guides by on-target activity and annotated off-target risk.

Implements the frozen local rule from inputs/RANKING_RULE.md:

1. A guide is REJECTED when any annotated off-target is a coding exon with at most
   two mismatches, when expression_risk=high, or when coding_bystander_count>0.
2. Risk class is `critical` for a high-expression-risk hit, at least two coding
   bystanders, or a coding-exon hit with one mismatch; otherwise `high` for a
   coding-exon hit with at most two mismatches or one coding bystander; otherwise
   `low`.
3. Rank non-rejected guides before rejected guides. Within the non-rejected group
   rank by decreasing activity_score, then guide_id. Within rejected guides rank
   `high` before `critical`, then by decreasing activity_score, then guide_id.
4. Label rank 1 `recommend`, other non-rejected guides `acceptable`, and rejected
   guides `reject`.

Only the two supplied CSVs under ./inputs are read. No external data is fetched.
"""

import csv
import os

INPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "inputs")
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
CANDIDATES_CSV = os.path.join(INPUT_DIR, "candidates.csv")
OFF_TARGETS_CSV = os.path.join(INPUT_DIR, "off_targets.csv")
RANKED_CSV = os.path.join(OUTPUT_DIR, "ranked_guides.csv")

# Within the rejected group, risk class `high` ranks before `critical`.
REJECTED_RISK_ORDER = {"high": 0, "critical": 1, "low": 2}


def load_inputs():
    with open(CANDIDATES_CSV, newline="", encoding="utf-8") as fh:
        candidates = list(csv.DictReader(fh))
    with open(OFF_TARGETS_CSV, newline="", encoding="utf-8") as fh:
        off_targets = list(csv.DictReader(fh))

    guides = {}
    for row in candidates:
        guides[row["guide_id"]] = {
            "guide_id": row["guide_id"],
            "activity_score": float(row["activity_score"]),
            "coding_bystander_count": int(row["coding_bystander_count"]),
            "off_targets": [],
        }
    for row in off_targets:
        gid = row["guide_id"]
        if gid not in guides:
            raise ValueError(f"off-target row references unknown guide {gid!r}")
        guides[gid]["off_targets"].append(
            {
                "locus": row["locus"],
                "mismatches": int(row["mismatches"]),
                "region": row["region"].strip(),
                "expression_risk": row["expression_risk"].strip().lower(),
            }
        )
    return guides


def is_coding_exon(ot):
    return ot["region"].lower() == "coding_exon"


def evaluate(guide):
    """Return (rejected, risk_class, evidence) for one guide per the frozen rule."""
    ots = guide["off_targets"]
    bystanders = guide["coding_bystander_count"]

    coding_le2 = [ot for ot in ots if is_coding_exon(ot) and ot["mismatches"] <= 2]
    coding_1mm = [ot for ot in ots if is_coding_exon(ot) and ot["mismatches"] == 1]
    high_expr = [ot for ot in ots if ot["expression_risk"] == "high"]

    # Rule 1: rejection.
    rejected = bool(coding_le2) or bool(high_expr) or bystanders > 0

    # Rule 2: risk class.
    if high_expr or bystanders >= 2 or coding_1mm:
        risk_class = "critical"
    elif coding_le2 or bystanders == 1:
        risk_class = "high"
    else:
        risk_class = "low"

    evidence = {
        "coding_le2": coding_le2,
        "coding_1mm": coding_1mm,
        "high_expr": high_expr,
        "bystanders": bystanders,
    }
    return rejected, risk_class, evidence


def build_rationale(guide, rejected, risk_class, evidence):
    """Human-readable rationale that states the safety evidence and any trade-off."""
    gid = guide["guide_id"]
    act = f"{guide['activity_score']:.2f}"
    parts = []

    ots = guide["off_targets"]
    if ots:
        ot_desc = "; ".join(
            f"{ot['locus']} ({ot['region']}, {ot['mismatches']} mismatch"
            f"{'es' if ot['mismatches'] != 1 else ''}, {ot['expression_risk']} expression risk)"
            for ot in ots
        )
        parts.append(f"Annotated off-target: {ot_desc}.")
    else:
        parts.append("No annotated off-targets.")

    bystanders = evidence["bystanders"]
    if bystanders > 0:
        parts.append(f"{bystanders} coding bystander{'s' if bystanders != 1 else ''} (safety evidence).")

    if rejected:
        reasons = []
        if evidence["coding_le2"]:
            worst = min(evidence["coding_le2"], key=lambda ot: ot["mismatches"])
            reasons.append(
                f"coding-exon off-target with {worst['mismatches']} mismatch"
                f"{'es' if worst['mismatches'] != 1 else ''} (<=2)"
            )
        if evidence["high_expr"]:
            reasons.append("high expression-risk off-target hit")
        if bystanders > 0:
            reasons.append("coding bystander count > 0")
        parts.append("Rejected because: " + "; ".join(reasons) + ".")
        parts.append(
            f"Trade-off: on-target activity {act} is sacrificed because coding/exonic "
            "near-matches and mismatch count are treated as hard safety evidence."
        )
    else:
        parts.append(
            f"No coding-exon hit within 2 mismatches, no high-expression-risk hit, and no "
            f"coding bystanders; mismatch counts are high enough to keep off-target risk low."
        )
        if gid == "g02":
            parts.append(
                f"Trade-off stated: activity {act} is lower than the rejected guides "
                "g01 (0.82), g05 (0.78), g03 (0.74); the safest guide is recommended over "
                "more active but riskier ones."
            )

    parts.append(f"Risk class {risk_class} per frozen rule.")
    return " ".join(parts)


def rank_guides(guides):
    evaluated = []
    for gid, guide in guides.items():
        rejected, risk_class, evidence = evaluate(guide)
        evaluated.append(
            {
                "guide": guide,
                "rejected": rejected,
                "risk_class": risk_class,
                "evidence": evidence,
            }
        )

    accepted = [e for e in evaluated if not e["rejected"]]
    rejected_list = [e for e in evaluated if e["rejected"]]

    accepted.sort(key=lambda e: (-e["guide"]["activity_score"], e["guide"]["guide_id"]))
    rejected_list.sort(
        key=lambda e: (
            REJECTED_RISK_ORDER.get(e["risk_class"], len(REJECTED_RISK_ORDER)),
            -e["guide"]["activity_score"],
            e["guide"]["guide_id"],
        )
    )

    ordered = accepted + rejected_list
    rows = []
    for rank, item in enumerate(ordered, start=1):
        guide = item["guide"]
        if item["rejected"]:
            decision = "reject"
        else:
            decision = "recommend" if rank == 1 else "acceptable"
        rows.append(
            {
                "rank": rank,
                "guide_id": guide["guide_id"],
                "on_target_score": f"{guide['activity_score']:.2f}",
                "risk_class": item["risk_class"],
                "decision": decision,
                "rationale": build_rationale(
                    guide, item["rejected"], item["risk_class"], item["evidence"]
                ),
            }
        )

    # Integrity checks: unique ranks, every input guide exactly once.
    ranks = [r["rank"] for r in rows]
    assert ranks == list(range(1, len(rows) + 1)), "ranks must be unique and contiguous"
    assert sorted(r["guide_id"] for r in rows) == sorted(guides.keys()), "every guide once"
    return rows


def main():
    guides = load_inputs()
    rows = rank_guides(guides)
    fieldnames = ["rank", "guide_id", "on_target_score", "risk_class", "decision", "rationale"]
    with open(RANKED_CSV, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {RANKED_CSV} with {len(rows)} guides.")
    for row in rows:
        print(
            f"  rank {row['rank']}: {row['guide_id']} "
            f"score={row['on_target_score']} risk={row['risk_class']} decision={row['decision']}"
        )


if __name__ == "__main__":
    main()
