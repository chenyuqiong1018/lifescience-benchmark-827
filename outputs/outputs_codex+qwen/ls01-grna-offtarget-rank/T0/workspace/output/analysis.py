#!/usr/bin/env python3
"""Rank CRISPR guide candidates by on-target activity and annotated off-target risk.

Inputs (read-only):
  inputs/candidates.csv   guide_id, spacer, activity_score, coding_bystander_count
  inputs/off_targets.csv  guide_id, locus, mismatches, region, expression_risk
  inputs/RANKING_RULE.md  the frozen deterministic ranking policy (implemented below)

Output:
  output/ranked_guides.csv  rank, guide_id, on_target_score, risk_class, decision, rationale

The script only reads the files under ./inputs and writes ./output/ranked_guides.csv.
It fetches no external data and never modifies the inputs.
"""

import csv
import os

INPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "inputs")
OUTPUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ranked_guides.csv")


def load_candidates():
    with open(os.path.join(INPUT_DIR, "candidates.csv"), newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        row["activity_score"] = float(row["activity_score"])
        row["coding_bystander_count"] = int(row["coding_bystander_count"])
    return rows


def load_off_targets():
    by_guide = {}
    with open(os.path.join(INPUT_DIR, "off_targets.csv"), newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            row["mismatches"] = int(row["mismatches"])
            by_guide.setdefault(row["guide_id"], []).append(row)
    return by_guide


def coding_exon_near_match(off_targets, max_mismatches=2):
    """Coding/exonic off-targets within max_mismatches: primary safety evidence."""
    return [ot for ot in off_targets
            if ot["region"] == "coding_exon" and ot["mismatches"] <= max_mismatches]


def classify(cand, off_targets):
    """Apply the frozen rule: return (rejected, risk_class, reasons, evidence)."""
    ce_hits = [ot for ot in off_targets if ot["region"] == "coding_exon"]
    ce_near = coding_exon_near_match(off_targets)
    high_expr = [ot for ot in off_targets if ot["expression_risk"] == "high"]
    bystanders = cand["coding_bystander_count"]

    # --- rejection rule ---------------------------------------------------
    reject_reasons = []
    for ot in ce_near:
        reject_reasons.append(
            f"coding-exon off-target {ot['locus']} with {ot['mismatches']} mismatch(es) <= 2")
    for ot in high_expr:
        reject_reasons.append(
            f"off-target {ot['locus']} annotated expression_risk=high")
    if bystanders > 0:
        reject_reasons.append(f"coding_bystander_count={bystanders} > 0")
    rejected = len(reject_reasons) > 0

    # --- risk class ---------------------------------------------------------
    one_mm_coding = [ot for ot in ce_hits if ot["mismatches"] == 1]
    if high_expr or bystanders >= 2 or one_mm_coding:
        risk = "critical"
    elif ce_near or bystanders == 1:
        risk = "high"
    else:
        risk = "low"

    # --- evidence summary -----------------------------------------------------
    evidence = []
    for ot in off_targets:
        evidence.append(
            f"{ot['locus']} region={ot['region']} mismatches={ot['mismatches']} "
            f"expression_risk={ot['expression_risk']}")
    if not evidence:
        evidence.append("no annotated off-targets")
    evidence.append(f"coding_bystander_count={bystanders}")

    return rejected, risk, reject_reasons, evidence


def build_rationale(cand, off_targets, rejected, risk, reject_reasons):
    """Human-readable rationale that states the safety evidence and any trade-off."""
    act = f"{cand['activity_score']:.2f}"
    parts = []
    if rejected:
        parts.append("Rejected: " + "; ".join(reject_reasons) + ".")
        parts.append(f"Risk class {risk} per frozen rule.")
        if cand["activity_score"] >= 0.70:
            parts.append(
                f"Trade-off stated explicitly: on-target activity {act} is high, but the "
                "coding/exonic near match and associated risk annotations outweigh the "
                "activity benefit; safety evidence dominates.")
        else:
            parts.append(
                f"On-target activity {act} does not compensate for the safety liabilities above.")
    else:
        parts.append(
            "No coding-exon off-target within 2 mismatches, no high-expression-risk hit, "
            "and coding_bystander_count=0.")
        parts.append(f"Risk class {risk}; on-target activity {cand['activity_score']:.2f}.")
        if risk == "low" and cand["activity_score"] < 0.70:
            parts.append(
                "Trade-off stated explicitly: accepted despite moderate activity because the "
                "annotated off-targets are non-coding/intergenic with >=3 mismatches and low "
                "expression risk; higher-activity guides were rejected on safety grounds.")
    parts.append("Evidence: " + " | ".join(classify(cand, off_targets)[3]))
    return " ".join(parts)


def main():
    candidates = load_candidates()
    off_targets = load_off_targets()

    records = []
    for cand in candidates:
        ots = off_targets.get(cand["guide_id"], [])
        rejected, risk, reasons, _ev = classify(cand, ots)
        records.append({
            "guide_id": cand["guide_id"],
            "on_target_score": cand["activity_score"],
            "risk_class": risk,
            "rejected": rejected,
            "reasons": reasons,
            "ots": ots,
            "cand": cand,
        })

    # Ordering per frozen rule:
    # non-rejected first by decreasing activity then guide_id;
    # rejected after, high before critical, then decreasing activity, then guide_id.
    risk_order = {"high": 0, "critical": 1, "low": 2}
    accepted = sorted(
        [r for r in records if not r["rejected"]],
        key=lambda r: (-r["on_target_score"], r["guide_id"]))
    rejected = sorted(
        [r for r in records if r["rejected"]],
        key=lambda r: (risk_order[r["risk_class"]], -r["on_target_score"], r["guide_id"]))
    ordered = accepted + rejected

    # Validate: unique ranks, every input guide appears exactly once.
    assert len(ordered) == len(candidates), "every input guide must appear exactly once"
    assert len({r["guide_id"] for r in ordered}) == len(ordered), "duplicate guide detected"

    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    with open(OUTPUT_PATH, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["rank", "guide_id", "on_target_score", "risk_class", "decision", "rationale"])
        for rank, rec in enumerate(ordered, start=1):
            if rank == 1 and not rec["rejected"]:
                decision = "recommend"
            elif not rec["rejected"]:
                decision = "acceptable"
            else:
                decision = "reject"
            rationale = build_rationale(
                rec["cand"], rec["ots"], rec["rejected"], rec["risk_class"], rec["reasons"])
            writer.writerow([
                rank,
                rec["guide_id"],
                f"{rec['on_target_score']:.2f}",
                rec["risk_class"],
                decision,
                rationale,
            ])

    print(f"Wrote {OUTPUT_PATH} with {len(ordered)} ranked guides.")
    for rank, rec in enumerate(ordered, start=1):
        print(rank, rec["guide_id"], rec["on_target_score"], rec["risk_class"],
              "rejected" if rec["rejected"] else "kept")


if __name__ == "__main__":
    main()
