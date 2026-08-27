from __future__ import annotations

import csv
from pathlib import Path


OUTPUT_DIR = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[6]
INPUT_DIR = REPO_ROOT / "inputs" / "ls01-grna-offtarget-rank"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def classify(candidate: dict[str, str], hits: list[dict[str, str]]) -> tuple[bool, str, str]:
    coding_bystanders = int(candidate["coding_bystander_count"])
    coding_near = [
        hit
        for hit in hits
        if hit["region"] == "coding_exon" and int(hit["mismatches"]) <= 2
    ]
    high_expression = [hit for hit in hits if hit["expression_risk"] == "high"]
    rejected = bool(coding_near or high_expression or coding_bystanders > 0)

    critical = (
        bool(high_expression)
        or coding_bystanders >= 2
        or any(int(hit["mismatches"]) == 1 for hit in coding_near)
    )
    if critical:
        risk_class = "critical"
    elif coding_near or coding_bystanders == 1:
        risk_class = "high"
    else:
        risk_class = "low"

    evidence: list[str] = []
    if coding_near:
        mismatch_text = "/".join(str(hit["mismatches"]) for hit in coding_near)
        evidence.append(f"coding-exon near match at {mismatch_text} mismatch(es)")
    else:
        min_mismatch = min(int(hit["mismatches"]) for hit in hits)
        regions = "/".join(sorted({hit["region"] for hit in hits}))
        evidence.append(f"nearest annotated hit has {min_mismatch} mismatches in {regions}")
    if high_expression:
        evidence.append("high expression-risk annotation")
    if coding_bystanders:
        evidence.append(f"{coding_bystanders} coding bystander(s)")

    activity = float(candidate["activity_score"])
    if rejected:
        tradeoff = (
            f"Activity {activity:.2f} does not outweigh the participant-visible safety rejection rule."
        )
    else:
        tradeoff = (
            f"Activity {activity:.2f} is compared only after clearing the annotated safety filters."
        )
    rationale = "; ".join(evidence) + ". " + tradeoff
    return rejected, risk_class, rationale


def main() -> None:
    candidates = read_csv(INPUT_DIR / "candidates.csv")
    off_targets = read_csv(INPUT_DIR / "off_targets.csv")
    hits_by_guide: dict[str, list[dict[str, str]]] = {}
    for hit in off_targets:
        hits_by_guide.setdefault(hit["guide_id"], []).append(hit)

    rows: list[dict[str, object]] = []
    for candidate in candidates:
        rejected, risk_class, rationale = classify(
            candidate, hits_by_guide[candidate["guide_id"]]
        )
        rows.append(
            {
                "guide_id": candidate["guide_id"],
                "on_target_score": float(candidate["activity_score"]),
                "risk_class": risk_class,
                "rejected": rejected,
                "rationale": rationale,
            }
        )

    risk_order = {"high": 0, "critical": 1}
    rows.sort(
        key=lambda row: (
            bool(row["rejected"]),
            risk_order.get(str(row["risk_class"]), 0) if row["rejected"] else 0,
            -float(row["on_target_score"]),
            str(row["guide_id"]),
        )
    )

    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
        row["decision"] = (
            "reject"
            if row["rejected"]
            else ("recommend" if rank == 1 else "acceptable")
        )

    fieldnames = [
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
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
