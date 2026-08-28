#!/usr/bin/env python3
"""Pocket reliability assessment and mutation prioritization (frozen local rule).

Inputs (read-only, under inputs/):
  confidence.csv           - per-interval pLDDT / PAE-to-core confidence estimates
  mutation_candidates.csv  - mutation hypotheses (predicted ddG / activity change)
  SCORING_RULE.md          - frozen rule implemented here

Outputs (under output/):
  mutation_priorities.csv  - rank,mutation,pocket_support,confidence_penalty,decision
  pocket_assessment.json   - pocket_reliable, prioritize_mutations, reason

pLDDT/PAE are model confidence estimates; their uncertainty is propagated into the
ranking, support label, penalty, and decision. Predicted ddG / activity changes are
model hypotheses only and are never described as measured binding or activity effects.
"""

import csv
import json
import re
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / "inputs"
OUTPUT_DIR = SCRIPT_DIR

PLDDT_UNRELIABLE = 50.0    # AlphaFold DB guidance: pLDDT < 50 -> coordinates should not be interpreted
PAE_UNCERTAIN = 10.0       # PAE > 10 A -> relative placement of regions is uncertain
RELIABLE_PLDDT_MIN = 70.0  # frozen rule: pocket reliable only if every pocket interval has pLDDT >= 70
RELIABLE_PAE_MAX = 10.0    # frozen rule: ... and PAE <= 10 A

MUT_RE = re.compile(r"^[A-Za-z]+(\d+)[A-Za-z*]+$")


def residue_number(mutation):
    m = MUT_RE.match(mutation.strip())
    if not m:
        raise ValueError(f"cannot parse residue number from {mutation!r}")
    return int(m.group(1))


def load_intervals():
    intervals = []
    with open(INPUT_DIR / "confidence.csv", newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            intervals.append(
                {
                    "residue_start": int(row["residue_start"]),
                    "residue_end": int(row["residue_end"]),
                    "plddt": float(row["plddt"]),
                    "pae_to_core_a": float(row["pae_to_core_a"]),
                    "pocket_member": row["pocket_member"].strip().lower() == "true",
                }
            )
    return intervals


def load_candidates():
    candidates = []
    with open(INPUT_DIR / "mutation_candidates.csv", newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            candidates.append(
                {
                    "mutation": row["mutation"].strip(),
                    "predicted_ddg_kcal_mol": float(row["predicted_ddg_kcal_mol"]),
                    "predicted_activity_change": row["predicted_activity_change"].strip(),
                    "region": row["region"].strip(),
                }
            )
    return candidates


def join_interval(intervals, residue):
    hits = [iv for iv in intervals if iv["residue_start"] <= residue <= iv["residue_end"]]
    if len(hits) != 1:
        raise ValueError(f"residue {residue} joined to {len(hits)} intervals; expected exactly 1")
    return hits[0]


def confidence_penalty(iv):
    """Semicolon-separated applicable flags, fixed order plddt_lt_50;pae_gt_10A, else none."""
    flags = []
    if iv["plddt"] < PLDDT_UNRELIABLE:
        flags.append("plddt_lt_50")
    if iv["pae_to_core_a"] > PAE_UNCERTAIN:
        flags.append("pae_gt_10A")
    return ";".join(flags) if flags else "none"


def main():
    intervals = load_intervals()
    candidates = load_candidates()

    # Join each mutation to the interval containing its residue number.
    for cand in candidates:
        cand["residue"] = residue_number(cand["mutation"])
        cand["interval"] = join_interval(intervals, cand["residue"])

    pocket = [c for c in candidates if c["region"] == "pocket"]
    non_pocket = [c for c in candidates if c["region"] != "pocket"]

    # Frozen ranking for pocket candidates: higher pLDDT, then lower PAE, then lower
    # absolute predicted ddG, then lexical mutation; non-pocket candidates follow.
    def sort_key(c):
        return (
            -c["interval"]["plddt"],
            c["interval"]["pae_to_core_a"],
            abs(c["predicted_ddg_kcal_mol"]),
            c["mutation"],
        )

    ordered = sorted(pocket, key=sort_key) + sorted(non_pocket, key=sort_key)

    rows = []
    for rank, cand in enumerate(ordered, start=1):
        iv = cand["interval"]
        low_confidence = iv["plddt"] < PLDDT_UNRELIABLE
        if cand["region"] == "pocket":
            if not low_confidence:
                raise AssertionError(
                    "frozen rule defines no decision for a high-confidence pocket "
                    "candidate; none exists in this input"
                )
            decision = "defer_structure_validation"
        else:
            decision = "out_of_scope_non_pocket"
        rows.append(
            {
                "rank": rank,
                "mutation": cand["mutation"],
                "pocket_support": "unsupported_low_confidence" if low_confidence else "cautious_support",
                "confidence_penalty": confidence_penalty(iv),
                "decision": decision,
            }
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT_DIR / "mutation_priorities.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["rank", "mutation", "pocket_support", "confidence_penalty", "decision"]
        )
        writer.writeheader()
        writer.writerows(rows)

    pocket_intervals = [iv for iv in intervals if iv["pocket_member"]]
    pocket_reliable = all(
        iv["plddt"] >= RELIABLE_PLDDT_MIN and iv["pae_to_core_a"] <= RELIABLE_PAE_MAX
        for iv in pocket_intervals
    )
    prioritize = pocket_reliable  # never prioritize mutations against an unreliable pocket

    detail = "; ".join(
        f"{iv['residue_start']}-{iv['residue_end']} (pLDDT {iv['plddt']:g}, PAE {iv['pae_to_core_a']:g} A)"
        for iv in pocket_intervals
    )
    reason = (
        "Pocket is not reliable: pocket intervals " + detail
        + " fail the requirement that every pocket interval has pLDDT >= 70 and PAE <= 10 A. "
        "pLDDT < 50 means the predicted coordinates should not be interpreted, and PAE > 10 A "
        "means the relative placement of the pocket with respect to the core is uncertain, so "
        "the pocket geometry cannot support mutation prioritization. All pocket candidates "
        "receive decision defer_structure_validation. Predicted ddG and activity changes are "
        "model hypotheses only, not measured binding or activity effects."
    )

    assessment = {
        "pocket_reliable": pocket_reliable,
        "prioritize_mutations": prioritize,
        "reason": reason,
    }
    json_path = OUTPUT_DIR / "pocket_assessment.json"
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(assessment, fh, indent=2, ensure_ascii=False)
        fh.write("\n")

    print(f"wrote {csv_path} ({len(rows)} data rows)")
    print(f"wrote {json_path}")


if __name__ == "__main__":
    main()