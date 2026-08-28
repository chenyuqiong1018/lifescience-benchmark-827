#!/usr/bin/env python3
"""Pocket reliability assessment and mutation prioritization.

Implements inputs/SCORING_RULE.md exactly:
  * Join each mutation to the confidence.csv interval containing its residue number.
  * Rank region=pocket candidates by higher pLDDT, then lower PAE, then lower
    absolute predicted ddG, then lexical mutation; non-pocket candidates follow.
  * pocket_support / confidence_penalty / decision use the frozen allowed values.

pLDDT/PAE uncertainty is propagated into the decisions, and predicted ddG /
activity changes are treated strictly as model hypotheses, never as measured
binding or activity effects.
"""

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INPUTS = ROOT / "inputs"
OUTPUT = ROOT / "output"

# Frozen rule thresholds.
PLDDT_RELIABLE_MIN = 70.0   # pocket reliable only if every pocket interval >= this
PAE_RELIABLE_MAX = 10.0     # ... and PAE <= this
PLDDT_LOW = 50.0            # below: coordinates should not be interpreted
PAE_HIGH = 10.0             # above: relative placement uncertain


def parse_residue(mutation: str) -> int:
    m = re.match(r"^[A-Za-z]+(\d+)[A-Za-z*]+$", mutation.strip())
    if not m:
        raise ValueError(f"cannot parse residue from mutation '{mutation}'")
    return int(m.group(1))


def load_confidence(path: Path):
    intervals = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            intervals.append({
                "residue_start": int(row["residue_start"]),
                "residue_end": int(row["residue_end"]),
                "plddt": float(row["plddt"]),
                "pae_to_core_a": float(row["pae_to_core_a"]),
                "pocket_member": row["pocket_member"].strip().lower() == "true",
            })
    return intervals


def join_interval(residue: int, intervals):
    hits = [iv for iv in intervals
            if iv["residue_start"] <= residue <= iv["residue_end"]]
    if len(hits) != 1:
        raise ValueError(f"residue {residue} matched {len(hits)} intervals")
    return hits[0]


def pocket_support_for(plddt: float) -> str:
    return "unsupported_low_confidence" if plddt < PLDDT_LOW else "cautious_support"


def confidence_penalty_for(plddt: float, pae: float) -> str:
    flags = []
    if plddt < PLDDT_LOW:
        flags.append("plddt_lt_50")
    if pae > PAE_HIGH:
        flags.append("pae_gt_10A")
    return ";".join(flags) if flags else "none"


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    intervals = load_confidence(INPUTS / "confidence.csv")

    with open(INPUTS / "mutation_candidates.csv", newline="") as fh:
        mutations = list(csv.DictReader(fh))

    records = []
    for row in mutations:
        residue = parse_residue(row["mutation"])
        iv = join_interval(residue, intervals)
        records.append({
            "mutation": row["mutation"].strip(),
            "region": row["region"].strip(),
            "residue": residue,
            "predicted_ddg_kcal_mol": float(row["predicted_ddg_kcal_mol"]),
            "predicted_activity_change": row["predicted_activity_change"].strip(),
            "interval": iv,
        })

    pocket = [r for r in records if r["region"] == "pocket"]
    non_pocket = [r for r in records if r["region"] != "pocket"]

    # Frozen ranking: higher pLDDT, lower PAE, lower |predicted ddG|, lexical.
    pocket.sort(key=lambda r: (-r["interval"]["plddt"],
                               r["interval"]["pae_to_core_a"],
                               abs(r["predicted_ddg_kcal_mol"]),
                               r["mutation"]))
    ordered = pocket + non_pocket

    # All pocket candidates in this input join low-confidence intervals; the
    # frozen rule defines a decision only for that case.
    assert all(r["interval"]["plddt"] < PLDDT_LOW for r in pocket), \
        "unexpected high-confidence pocket candidate; frozen rule undefined"

    csv_rows = []
    for rank, rec in enumerate(ordered, start=1):
        iv = rec["interval"]
        decision = ("defer_structure_validation" if rec["region"] == "pocket"
                    else "out_of_scope_non_pocket")
        csv_rows.append({
            "rank": rank,
            "mutation": rec["mutation"],
            "pocket_support": pocket_support_for(iv["plddt"]),
            "confidence_penalty": confidence_penalty_for(iv["plddt"], iv["pae_to_core_a"]),
            "decision": decision,
        })

    with open(OUTPUT / "mutation_priorities.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["rank", "mutation", "pocket_support",
                                                "confidence_penalty", "decision"],
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(csv_rows)

    pocket_intervals = [iv for iv in intervals if iv["pocket_member"]]
    pocket_reliable = bool(pocket_intervals) and all(
        iv["plddt"] >= PLDDT_RELIABLE_MIN and iv["pae_to_core_a"] <= PAE_RELIABLE_MAX
        for iv in pocket_intervals)
    prioritize_mutations = pocket_reliable  # never prioritize on an unreliable pocket

    detail = "; ".join(
        f"{iv['residue_start']}-{iv['residue_end']} pLDDT {iv['plddt']:g}, "
        f"PAE {iv['pae_to_core_a']:g} A" for iv in pocket_intervals)
    reason = (
        "Every pocket interval fails the reliability thresholds "
        "(pLDDT >= 70 and PAE <= 10 A): " + detail + ". "
        "Both pocket intervals fall below pLDDT 50, the range AlphaFold guidance "
        "says should not be interpreted, and their high PAE means the relative "
        "placement of the pocket is uncertain. All pocket mutation candidates "
        "therefore receive defer_structure_validation. Predicted ddG and activity "
        "changes are model hypotheses conditional on this unreliable structure; "
        "they are not measured binding or activity effects.")

    assessment = {
        "pocket_reliable": pocket_reliable,
        "prioritize_mutations": prioritize_mutations,
        "reason": reason,
    }
    with open(OUTPUT / "pocket_assessment.json", "w", encoding="utf-8") as fh:
        json.dump(assessment, fh, indent=2, ensure_ascii=False)
        fh.write("\n")

    print("wrote", OUTPUT / "mutation_priorities.csv")
    print("wrote", OUTPUT / "pocket_assessment.json")
    for r in csv_rows:
        print(r)
    print(json.dumps(assessment, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
