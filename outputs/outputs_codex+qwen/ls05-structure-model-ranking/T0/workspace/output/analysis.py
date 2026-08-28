#!/usr/bin/env python3
"""Rank structural models exactly per inputs/SCORING_RULE.md.

This fixture is a benchmark-informed local extension: the supplied metrics are
already-computed comparison values against one frozen reference. Only files under
./inputs are read. No coordinate-level, interface, or experimental properties are
claimed; chain_mapping_complete is used only as the mapping-completeness flag and
is exposed as interface_score per the required output schema.
"""

import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
WORKSPACE = os.path.dirname(HERE)
INPUT_DIR = os.path.join(WORKSPACE, "inputs")
OUTPUT_DIR = os.path.join(WORKSPACE, "output")

# Critical region (residues 181-240) from SCORING_RULE.md.
CRITICAL_START = 181
CRITICAL_END = 240


def read_model_metrics(path):
    models = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            models[row["model_id"]] = {
                "model_id": row["model_id"],
                "tm_score": float(row["tm_score"]),
                "lddt": float(row["lddt"]),
                "rmsd_a": float(row["rmsd_a"]),
                "aligned_residues": int(row["aligned_residues"]),
                "chain_mapping_complete": row["chain_mapping_complete"].strip().lower() == "true",
            }
    return models


def read_residue_errors(path):
    errors = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            errors.setdefault(row["model_id"], []).append(
                (int(row["residue_start"]), int(row["residue_end"]), float(row["mean_error_a"]))
            )
    return errors


def critical_region_risk(rows, start=CRITICAL_START, end=CRITICAL_END):
    """Length-weighted mean of mean_error_a over rows overlapping [start, end].

    Each critical-region row in this fixture spans the whole interval, so the
    reported risk equals that row's value; the general length-weighted form is
    retained for correctness.
    """
    num = 0.0
    den = 0
    for r_start, r_end, mean_error in rows:
        lo = max(r_start, start)
        hi = min(r_end, end)
        if hi >= lo:
            weight = hi - lo + 1
            num += weight * mean_error
            den += weight
    if den == 0:
        return None
    return num / den


def fmt(value, places=6):
    text = f"{value:.{places}f}".rstrip("0").rstrip(".")
    return text if text else "0"


def rank_models(models, errors):
    records = []
    for model_id in sorted(models):
        m = models[model_id]
        risk = critical_region_risk(errors.get(model_id, []))
        records.append((m, risk))

    def sort_key(item):
        m, risk = item
        return (
            0 if m["chain_mapping_complete"] else 1,       # 1. complete before incomplete
            -m["tm_score"],                                # 2. higher tm_score
            -m["lddt"],                                    # 3. higher lddt
            m["rmsd_a"],                                   # 4. lower rmsd_a
            -m["aligned_residues"],                        # 5. more aligned_residues
            risk if risk is not None else float("inf"),    # 6. lower critical-region error
            m["model_id"],                                 # 7. lexical model_id tie break
        )

    ordered = sorted(records, key=sort_key)
    ranked = []
    for rank, (m, risk) in enumerate(ordered, start=1):
        if not m["chain_mapping_complete"]:
            decision = "reject_incomplete_mapping"
        elif rank == 1:
            decision = "preferred"
        else:
            decision = "alternate"
        ranked.append({
            "rank": rank,
            "model_id": m["model_id"],
            "global_score": fmt(m["tm_score"], 4),
            "interface_score": 1 if m["chain_mapping_complete"] else 0,
            "critical_residue_risk": fmt(risk, 6) if risk is not None else "",
            "decision": decision,
        })
    return ranked


def main():
    models = read_model_metrics(os.path.join(INPUT_DIR, "model_metrics.csv"))
    errors = read_residue_errors(os.path.join(INPUT_DIR, "residue_errors.csv"))
    ranked = rank_models(models, errors)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    out_path = os.path.join(OUTPUT_DIR, "model_ranking.csv")
    fieldnames = ["rank", "model_id", "global_score", "interface_score",
                  "critical_residue_risk", "decision"]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in ranked:
            writer.writerow(row)

    for row in ranked:
        print(row)
    print("WROTE", out_path)


if __name__ == "__main__":
    main()
