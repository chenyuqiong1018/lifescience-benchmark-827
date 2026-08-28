"""Rank structural models per inputs/SCORING_RULE.md (frozen local-extension rule).

Reads inputs/model_metrics.csv and inputs/residue_errors.csv, ranks models by the
stable tuple defined in the rule, and writes output/model_ranking.csv.

No coordinate-level, interface, or experimental claims are made: only the supplied
comparison metrics are used.
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CRITICAL_START, CRITICAL_END = 181, 240  # residues 181--240 per rule


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def critical_region_risk(model_id, residue_rows):
    """Length-weighted mean of mean_error_a over rows overlapping residues 181-240.

    In this fixture each model has exactly one critical-region row covering the
    whole interval, so the reported risk equals that row's value.
    """
    overlaps = []
    for r in residue_rows:
        if r["model_id"] != model_id:
            continue
        start, end = int(r["residue_start"]), int(r["residue_end"])
        ov = min(end, CRITICAL_END) - max(start, CRITICAL_START) + 1
        if ov > 0:
            overlaps.append((ov, float(r["mean_error_a"])))
    if not overlaps:
        raise ValueError(f"no critical-region rows for {model_id}")
    total_len = sum(w for w, _ in overlaps)
    return sum(w * v for w, v in overlaps) / total_len


def main():
    metrics = read_csv(os.path.join(ROOT, "inputs", "model_metrics.csv"))
    errors = read_csv(os.path.join(ROOT, "inputs", "residue_errors.csv"))

    for m in metrics:
        m["tm_score"] = float(m["tm_score"])
        m["lddt"] = float(m["lddt"])
        m["rmsd_a"] = float(m["rmsd_a"])
        m["aligned_residues"] = int(m["aligned_residues"])
        m["chain_mapping_complete"] = m["chain_mapping_complete"].strip().lower() == "true"
        m["critical_risk"] = critical_region_risk(m["model_id"], errors)

    # Stable tuple from SCORING_RULE.md:
    # 1 complete chain mapping first; 2 higher tm_score; 3 higher lddt;
    # 4 lower rmsd_a; 5 more aligned_residues; 6 lower critical-region mean error;
    # 7 lexical model_id.
    ranked = sorted(
        metrics,
        key=lambda m: (
            not m["chain_mapping_complete"],
            -m["tm_score"],
            -m["lddt"],
            m["rmsd_a"],
            -m["aligned_residues"],
            m["critical_risk"],
            m["model_id"],
        ),
    )

    out_path = os.path.join(HERE, "model_ranking.csv")
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["rank", "model_id", "global_score", "interface_score",
                    "critical_residue_risk", "decision"])
        for i, m in enumerate(ranked, start=1):
            if not m["chain_mapping_complete"]:
                decision = "reject_incomplete_mapping"
            elif i == 1:
                decision = "preferred"
            else:
                decision = "alternate"
            w.writerow([
                i,
                m["model_id"],
                m["tm_score"],
                1 if m["chain_mapping_complete"] else 0,
                m["critical_risk"],
                decision,
            ])
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
