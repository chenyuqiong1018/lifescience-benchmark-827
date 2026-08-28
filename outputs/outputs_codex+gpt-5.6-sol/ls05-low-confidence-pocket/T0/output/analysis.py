from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path


OUTPUT_DIR = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[5]
INPUT_DIR = REPO / "inputs" / "ls05-low-confidence-pocket"

CONFIDENCE_FIELDS = ("residue_start", "residue_end", "plddt", "pae_to_core_a", "pocket_member")
CANDIDATE_FIELDS = ("mutation", "predicted_ddg_kcal_mol", "predicted_activity_change", "region")


@dataclass(frozen=True)
class Interval:
    start: int
    end: int
    plddt: Decimal
    pae: Decimal
    pocket_member: bool


@dataclass(frozen=True)
class Candidate:
    mutation: str
    residue: int
    ddg_hypothesis: Decimal
    activity_hypothesis: str
    region: str
    interval: Interval

    def rank_key(self) -> tuple[object, ...]:
        if self.region == "pocket":
            return (0, -self.interval.plddt, self.interval.pae, abs(self.ddg_hypothesis), self.mutation)
        return (1, self.mutation)


def read_rows(path: Path, expected: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != expected:
            raise ValueError(f"Unexpected schema in {path.name}: {reader.fieldnames}")
        rows = list(reader)
    if not rows or any(value == "" for row in rows for value in row.values()):
        raise ValueError(f"Missing data in {path.name}")
    return rows


confidence_rows = read_rows(INPUT_DIR / "confidence.csv", CONFIDENCE_FIELDS)
candidate_rows = read_rows(INPUT_DIR / "mutation_candidates.csv", CANDIDATE_FIELDS)

intervals = sorted(
    (
        Interval(
            start=int(row["residue_start"]),
            end=int(row["residue_end"]),
            plddt=Decimal(row["plddt"]),
            pae=Decimal(row["pae_to_core_a"]),
            pocket_member=row["pocket_member"].lower() == "true",
        )
        for row in confidence_rows
    ),
    key=lambda interval: interval.start,
)
for index, interval in enumerate(intervals):
    if interval.start > interval.end:
        raise ValueError(f"Reversed interval: {interval}")
    if index and intervals[index - 1].end >= interval.start:
        raise ValueError("Confidence intervals overlap")

candidates: list[Candidate] = []
for row in candidate_rows:
    match = re.fullmatch(r"[A-Za-z](\d+)[A-Za-z]", row["mutation"])
    if not match:
        raise ValueError(f"Unsupported mutation: {row['mutation']}")
    residue = int(match.group(1))
    joined = [interval for interval in intervals if interval.start <= residue <= interval.end]
    if len(joined) != 1:
        raise ValueError(f"Expected one interval for {row['mutation']}; found {len(joined)}")
    candidate = Candidate(
        mutation=row["mutation"],
        residue=residue,
        ddg_hypothesis=Decimal(row["predicted_ddg_kcal_mol"]),
        activity_hypothesis=row["predicted_activity_change"],
        region=row["region"],
        interval=joined[0],
    )
    if (candidate.region == "pocket") != candidate.interval.pocket_member:
        raise ValueError(f"Region/pocket membership mismatch for {candidate.mutation}")
    candidates.append(candidate)

ranked = sorted(candidates, key=Candidate.rank_key)
records: list[dict[str, object]] = []
for rank, candidate in enumerate(ranked, start=1):
    penalties = []
    if candidate.interval.plddt < 50:
        penalties.append("plddt_lt_50")
    if candidate.interval.pae > 10:
        penalties.append("pae_gt_10A")
    records.append(
        {
            "rank": rank,
            "mutation": candidate.mutation,
            "pocket_support": (
                "unsupported_low_confidence" if candidate.interval.plddt < 50 else "cautious_support"
            ),
            "confidence_penalty": ";".join(penalties) if penalties else "none",
            "decision": (
                "out_of_scope_non_pocket"
                if candidate.region != "pocket"
                else "defer_structure_validation"
                if candidate.interval.plddt < 50
                else "cautious_support"
            ),
        }
    )

headers = ["rank", "mutation", "pocket_support", "confidence_penalty", "decision"]
with (OUTPUT_DIR / "mutation_priorities.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=headers)
    writer.writeheader()
    writer.writerows(records)

pocket_intervals = [interval for interval in intervals if interval.pocket_member]
reliable = bool(pocket_intervals) and all(
    interval.plddt >= 70 and interval.pae <= 10 for interval in pocket_intervals
)
assessment = {
    "pocket_reliable": reliable,
    "prioritize_mutations": reliable,
    "reason": (
        "The two pocket intervals both fail the required confidence criteria: 210-230 has pLDDT 47 "
        "and PAE 14.2 A; 231-270 has pLDDT 43 and PAE 17.8 A. Mutation hypotheses should "
        "remain deferred until the pocket structure is validated."
    ),
}
(OUTPUT_DIR / "pocket_assessment.json").write_text(
    json.dumps(assessment, indent=2) + "\n", encoding="utf-8"
)

table_rows = "\n".join(
    f"| {rank} | {candidate.mutation} | {candidate.region} | {candidate.interval.start}–{candidate.interval.end} | "
    f"{candidate.interval.plddt} | {candidate.interval.pae} | {candidate.ddg_hypothesis} | "
    f"{candidate.activity_hypothesis} | {records[rank - 1]['decision']} |"
    for rank, candidate in enumerate(ranked, 1)
)
report = f"""# Low-confidence pocket assessment

## Overall assessment

The nominated pocket is **not reliable enough to prioritize mutations**. This is a confidence-limited structural hypothesis, not evidence that any mutation changes binding, stability, or activity.

## Evidence audit

**Strengths.** Both input schemas are complete, the four confidence intervals are ordered and non-overlapping, every mutation maps to exactly one interval, and each `region` label agrees with the joined `pocket_member` flag. The ranking rule is prespecified and deterministic.

**Critical concern.** Both pocket intervals fail both confidence criteria. Residues 210–230 have pLDDT 47 and PAE 14.2 Å; residues 231–270 have pLDDT 43 and PAE 17.8 Å. Low local confidence weakens coordinate interpretation, while high PAE weakens confidence in placement relative to the core. The combined limitation is severe enough to defer every pocket candidate.

## Ranked hypotheses

| Rank | Mutation | Region | Joined interval | pLDDT | PAE to core (Å) | Predicted ΔΔG hypothesis (kcal/mol) | Activity hypothesis | Decision |
|---:|---|---|---|---:|---:|---:|---|---|
{table_rows}

Y218F ranks first because its joined interval has higher pLDDT and lower PAE than the interval for R244A and L265W. R244A and L265W tie on confidence, so the lower absolute predicted ΔΔG hypothesis places R244A first. S301A is non-pocket and is placed last by rule.

## Claim boundary and recommendation

The ΔΔG and activity columns are unvalidated model hypotheses. They are not measured effects and cannot rescue an uninterpretable pocket. The proportionate next step is structure validation or a better-supported structural model before using these hypotheses to prioritize mutagenesis; experimental measurements would then be needed to establish actual stability, binding, or activity effects.
"""
(OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")

print(" > ".join(candidate.mutation for candidate in ranked))
