from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


OUTPUT_DIR = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[6]
INPUT_DIR = REPO_ROOT / "inputs" / "ls01-grna-offtarget-rank"
DNA_ALPHABET = frozenset("ACGT")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def validate_inputs(
    candidates: list[dict[str, str]], off_targets: list[dict[str, str]]
) -> None:
    guide_ids = [row["guide_id"] for row in candidates]
    if len(guide_ids) != len(set(guide_ids)):
        raise ValueError("candidate guide IDs must be unique")
    for row in candidates:
        spacer = row["spacer"].upper()
        if len(spacer) != 20 or not set(spacer).issubset(DNA_ALPHABET):
            raise ValueError(f"invalid 20-nt DNA spacer for {row['guide_id']}")
    annotated_ids = {row["guide_id"] for row in off_targets}
    if set(guide_ids) != annotated_ids:
        raise ValueError("candidate and off-target guide sets must match exactly")


def evaluate(
    candidate: dict[str, str], hits: list[dict[str, str]]
) -> dict[str, object]:
    coding_bystanders = int(candidate["coding_bystander_count"])
    coding_near_hits = [
        hit
        for hit in hits
        if hit["region"] == "coding_exon" and int(hit["mismatches"]) <= 2
    ]
    high_expression_hits = [
        hit for hit in hits if hit["expression_risk"] == "high"
    ]
    rejected = bool(
        coding_near_hits or high_expression_hits or coding_bystanders > 0
    )

    if (
        high_expression_hits
        or coding_bystanders >= 2
        or any(int(hit["mismatches"]) == 1 for hit in coding_near_hits)
    ):
        risk_class = "critical"
    elif coding_near_hits or coding_bystanders == 1:
        risk_class = "high"
    else:
        risk_class = "low"

    evidence: list[str] = []
    if coding_near_hits:
        details = ", ".join(
            f"{hit['locus']} ({hit['mismatches']} mismatch(es), coding exon)"
            for hit in coding_near_hits
        )
        evidence.append(f"annotated near match: {details}")
    else:
        nearest = min(hits, key=lambda hit: int(hit["mismatches"]))
        evidence.append(
            "nearest annotated hit: "
            f"{nearest['locus']} ({nearest['mismatches']} mismatches, {nearest['region']})"
        )
    if high_expression_hits:
        evidence.append("expression risk is high")
    if coding_bystanders:
        evidence.append(f"coding bystander count is {coding_bystanders}")

    activity = float(candidate["activity_score"])
    tradeoff = (
        f"on-target activity {activity:.2f} cannot override a safety rejection"
        if rejected
        else f"on-target activity {activity:.2f} is ranked after safety clearance"
    )
    return {
        "guide_id": candidate["guide_id"],
        "on_target_score": activity,
        "risk_class": risk_class,
        "rejected": rejected,
        "rationale": "; ".join(evidence + [tradeoff]) + ".",
    }


def main() -> None:
    candidates = read_csv(INPUT_DIR / "candidates.csv")
    off_targets = read_csv(INPUT_DIR / "off_targets.csv")
    validate_inputs(candidates, off_targets)

    by_guide: dict[str, list[dict[str, str]]] = defaultdict(list)
    for hit in off_targets:
        by_guide[hit["guide_id"]].append(hit)
    rows = [evaluate(candidate, by_guide[candidate["guide_id"]]) for candidate in candidates]

    rejected_risk_order = {"high": 0, "critical": 1}
    rows.sort(
        key=lambda row: (
            bool(row["rejected"]),
            rejected_risk_order.get(str(row["risk_class"]), 0)
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
    with (OUTPUT_DIR / "ranked_guides.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
