"""Rank structural models per inputs/SCORING_RULE.md (frozen local-extension rule).

Reads precomputed comparison metrics against one frozen reference
(inputs/model_metrics.csv, inputs/residue_errors.csv). No coordinate-level,
interface, or experimental claims are made; only the supplied metrics are used.

Outputs: output/model_ranking.csv
Columns: rank,model_id,global_score,interface_score,critical_residue_risk,decision
"""

import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(HERE, "..", "inputs")
METRICS_CSV = os.path.join(INPUT_DIR, "model_metrics.csv")
RESIDUE_ERRORS_CSV = os.path.join(INPUT_DIR, "residue_errors.csv")
OUT_CSV = os.path.join(HERE, "model_ranking.csv")

CRITICAL_START = 181
CRITICAL_END = 240


def load_rows(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def critical_risk(model_id, residue_rows):
    """Length-weighted mean of mean_error_a over rows overlapping residues 181-240."""
    overlaps = []
    for r in residue_rows:
        if r["model_id"] != model_id:
            continue
        start, end = int(r["residue_start"]), int(r["residue_end"])
        ov_start, ov_end = max(start, CRITICAL_START), min(end, CRITICAL_END)
        length = ov_end - ov_start + 1
        if length > 0:
            overlaps.append((length, float(r["mean_error_a"])))
    if not overlaps:
        raise ValueError(f"no critical-region rows for {model_id}")
    total_len = sum(l for l, _ in overlaps)
    return sum(l * e for l, e in overlaps) / total_len


def main():
    metrics = load_rows(METRICS_CSV)
    residue_rows = load_rows(RESIDUE_ERRORS_CSV)

    records = []
    for m in metrics:
        complete = m["chain_mapping_complete"].strip().lower() == "true"
        records.append(
            {
                "model_id": m["model_id"],
                "complete": complete,
                "tm_score": float(m["tm_score"]),
                "lddt": float(m["lddt"]),
                "rmsd_a": float(m["rmsd_a"]),
                "aligned_residues": int(m["aligned_residues"]),
                "risk": critical_risk(m["model_id"], residue_rows),
            }
        )

    # Stable ranking tuple:
    # 1) complete mapping first; 2) higher tm_score; 3) higher lddt;
    # 4) lower rmsd_a; 5) more aligned_residues; 6) lower critical-region
    # mean error; 7) lexical model_id.
    records.sort(
        key=lambda r: (
            0 if r["complete"] else 1,
            -r["tm_score"],
            -r["lddt"],
            r["rmsd_a"],
            -r["aligned_residues"],
            r["risk"],
            r["model_id"],
        )
    )

    os.makedirs(HERE, exist_ok=True)
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "rank",
                "model_id",
                "global_score",
                "interface_score",
                "critical_residue_risk",
                "decision",
            ]
        )
        for rank, r in enumerate(records, start=1):
            if not r["complete"]:
                decision = "reject_incomplete_mapping"
            elif rank == 1:
                decision = "preferred"
            else:
                decision = "alternate"
            writer.writerow(
                [
                    rank,
                    r["model_id"],
                    f"{r['tm_score']:g}",
                    1 if r["complete"] else 0,
                    f"{r['risk']:g}",
                    decision,
                ]
            )

    print(f"wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
