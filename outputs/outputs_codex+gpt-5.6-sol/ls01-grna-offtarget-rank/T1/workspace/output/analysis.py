from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


OUT = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[6]
DATA = ROOT / "inputs" / "ls01-grna-offtarget-rank"


def load(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    candidates = load(DATA / "candidates.csv")
    hits = load(DATA / "off_targets.csv")
    candidate_ids = [row["guide_id"] for row in candidates]
    if len(candidate_ids) != len(set(candidate_ids)):
        raise ValueError("duplicate guide_id in candidates")
    for row in candidates:
        spacer = row["spacer"].upper()
        if len(spacer) != 20 or set(spacer) - set("ACGT"):
            raise ValueError(f"invalid DNA spacer: {row['guide_id']}")

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for hit in hits:
        grouped[hit["guide_id"]].append(hit)
    if set(grouped) != set(candidate_ids):
        raise ValueError("off-target annotations must cover every guide exactly by ID set")

    rows: list[dict[str, object]] = []
    for candidate in candidates:
        guide_hits = grouped[candidate["guide_id"]]
        bystanders = int(candidate["coding_bystander_count"])
        coding_near = [
            hit
            for hit in guide_hits
            if hit["region"] == "coding_exon" and int(hit["mismatches"]) <= 2
        ]
        high_expression = [
            hit for hit in guide_hits if hit["expression_risk"] == "high"
        ]
        rejected = bool(coding_near or high_expression or bystanders > 0)

        if (
            high_expression
            or bystanders >= 2
            or any(int(hit["mismatches"]) == 1 for hit in coding_near)
        ):
            risk = "critical"
        elif coding_near or bystanders == 1:
            risk = "high"
        else:
            risk = "low"

        activity = float(candidate["activity_score"])
        if coding_near:
            safety = "; ".join(
                f"coding-exon hit {hit['locus']} has {hit['mismatches']} mismatch(es)"
                for hit in coding_near
            )
        else:
            nearest = min(guide_hits, key=lambda hit: int(hit["mismatches"]))
            safety = (
                f"nearest hit {nearest['locus']} has {nearest['mismatches']} mismatches "
                f"in {nearest['region']}"
            )
        extras: list[str] = []
        if high_expression:
            extras.append("high expression risk")
        if bystanders:
            extras.append(f"{bystanders} coding bystander(s)")
        tradeoff = (
            f"activity {activity:.2f} does not override safety rejection"
            if rejected
            else f"activity {activity:.2f} is ranked after safety clearance"
        )
        rationale = "; ".join([safety, *extras, tradeoff]) + "."
        rows.append(
            {
                "guide_id": candidate["guide_id"],
                "on_target_score": activity,
                "risk_class": risk,
                "rejected": rejected,
                "rationale": rationale,
            }
        )

    rejected_risk = {"high": 0, "critical": 1}
    rows.sort(
        key=lambda row: (
            bool(row["rejected"]),
            rejected_risk.get(str(row["risk_class"]), 0)
            if row["rejected"]
            else 0,
            -float(row["on_target_score"]),
            str(row["guide_id"]),
        )
    )
    for rank, row in enumerate(rows, 1):
        row["rank"] = rank
        row["decision"] = (
            "reject"
            if row["rejected"]
            else ("recommend" if rank == 1 else "acceptable")
        )

    fields = [
        "rank",
        "guide_id",
        "on_target_score",
        "risk_class",
        "decision",
        "rationale",
    ]
    with (OUT / "ranked_guides.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
