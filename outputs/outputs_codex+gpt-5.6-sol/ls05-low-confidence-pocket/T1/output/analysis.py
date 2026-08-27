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
PLDDT_RELIABLE = Decimal("70")
PAE_RELIABLE = Decimal("10")


@dataclass(frozen=True)
class ConfidenceInterval:
    start: int
    end: int
    plddt: Decimal
    pae_to_core: Decimal
    pocket_member: bool

    @property
    def reliable(self) -> bool:
        return self.plddt >= PLDDT_RELIABLE and self.pae_to_core <= PAE_RELIABLE


@dataclass(frozen=True)
class MutationHypothesis:
    mutation: str
    residue: int
    predicted_ddg: Decimal
    predicted_activity: str
    region: str
    confidence: ConfidenceInterval

    def ordering_key(self) -> tuple[object, ...]:
        if self.region == "pocket":
            return (
                0,
                -self.confidence.plddt,
                self.confidence.pae_to_core,
                abs(self.predicted_ddg),
                self.mutation,
            )
        return (1, self.mutation)


def load_csv(path: Path, exact_fields: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != exact_fields:
            raise ValueError(f"Schema mismatch in {path.name}: {reader.fieldnames}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"No data in {path.name}")
    if any(value.strip() == "" for row in rows for value in row.values()):
        raise ValueError(f"Missing value in {path.name}")
    return rows


confidence_rows = load_csv(
    INPUT_DIR / "confidence.csv",
    ("residue_start", "residue_end", "plddt", "pae_to_core_a", "pocket_member"),
)
candidate_rows = load_csv(
    INPUT_DIR / "mutation_candidates.csv",
    ("mutation", "predicted_ddg_kcal_mol", "predicted_activity_change", "region"),
)

intervals: list[ConfidenceInterval] = []
for row in confidence_rows:
    member = row["pocket_member"].lower()
    if member not in {"true", "false"}:
        raise ValueError("Invalid pocket_member flag")
    interval = ConfidenceInterval(
        start=int(row["residue_start"]),
        end=int(row["residue_end"]),
        plddt=Decimal(row["plddt"]),
        pae_to_core=Decimal(row["pae_to_core_a"]),
        pocket_member=member == "true",
    )
    if interval.start > interval.end:
        raise ValueError(f"Reversed interval: {interval}")
    intervals.append(interval)
intervals.sort(key=lambda item: item.start)
if any(left.end >= right.start for left, right in zip(intervals, intervals[1:])):
    raise ValueError("Confidence intervals overlap")

hypotheses: list[MutationHypothesis] = []
for row in candidate_rows:
    match = re.fullmatch(r"[A-Za-z](\d+)[A-Za-z]", row["mutation"])
    if not match:
        raise ValueError(f"Unsupported mutation notation: {row['mutation']}")
    residue = int(match.group(1))
    matches = [interval for interval in intervals if interval.start <= residue <= interval.end]
    if len(matches) != 1:
        raise ValueError(f"Confidence join failed for {row['mutation']}: {len(matches)} matches")
    hypothesis = MutationHypothesis(
        mutation=row["mutation"],
        residue=residue,
        predicted_ddg=Decimal(row["predicted_ddg_kcal_mol"]),
        predicted_activity=row["predicted_activity_change"],
        region=row["region"],
        confidence=matches[0],
    )
    if (hypothesis.region == "pocket") != hypothesis.confidence.pocket_member:
        raise ValueError(f"Region and pocket membership disagree for {hypothesis.mutation}")
    hypotheses.append(hypothesis)
if len({item.mutation for item in hypotheses}) != len(hypotheses):
    raise ValueError("Duplicate mutation")

ranked = sorted(hypotheses, key=MutationHypothesis.ordering_key)
output_rows: list[dict[str, object]] = []
for rank, item in enumerate(ranked, start=1):
    penalties: list[str] = []
    if item.confidence.plddt < 50:
        penalties.append("plddt_lt_50")
    if item.confidence.pae_to_core > 10:
        penalties.append("pae_gt_10A")
    output_rows.append(
        {
            "rank": rank,
            "mutation": item.mutation,
            "pocket_support": (
                "unsupported_low_confidence" if item.confidence.plddt < 50 else "cautious_support"
            ),
            "confidence_penalty": ";".join(penalties) if penalties else "none",
            "decision": (
                "out_of_scope_non_pocket"
                if item.region != "pocket"
                else "defer_structure_validation"
                if item.confidence.plddt < 50
                else "cautious_support"
            ),
        }
    )

fields = ["rank", "mutation", "pocket_support", "confidence_penalty", "decision"]
with (OUTPUT_DIR / "mutation_priorities.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(output_rows)

pocket_intervals = [interval for interval in intervals if interval.pocket_member]
pocket_reliable = bool(pocket_intervals) and all(interval.reliable for interval in pocket_intervals)
assessment = {
    "pocket_reliable": pocket_reliable,
    "prioritize_mutations": pocket_reliable,
    "reason": (
        "Pocket reliability requires every pocket interval to have pLDDT >= 70 and PAE <= 10 A. "
        "Residues 210-230 have pLDDT 47 and PAE 14.2 A; residues 231-270 have pLDDT 43 "
        "and PAE 17.8 A. Both intervals fail both criteria, so mutation prioritization is deferred."
    ),
}
(OUTPUT_DIR / "pocket_assessment.json").write_text(
    json.dumps(assessment, indent=2) + "\n", encoding="utf-8"
)

rank_rows = "\n".join(
    f"| {rank} | {item.mutation} | {item.region} | {item.confidence.start}–{item.confidence.end} | "
    f"{item.confidence.plddt} | {item.confidence.pae_to_core} | "
    f"{item.confidence.plddt - PLDDT_RELIABLE:+} | {PAE_RELIABLE - item.confidence.pae_to_core:+} | "
    f"{item.predicted_ddg} | {item.predicted_activity} | {output_rows[rank - 1]['decision']} |"
    for rank, item in enumerate(ranked, start=1)
)
report = f"""# Low-confidence pocket assessment

## Answer

The nominated pocket is **not reliable enough to prioritize mutations**. Both pocket intervals fail the pLDDT and PAE requirements, so all three pocket candidates are deferred for structure validation.

## Candidate-level uncertainty propagation

| Rank | Mutation | Region | Joined interval | pLDDT | PAE (Å) | pLDDT margin to 70 | PAE margin to ≤10 Å | Predicted ΔΔG hypothesis (kcal/mol) | Activity hypothesis | Decision |
|---:|---|---|---|---:|---:|---:|---:|---:|---|---|
{rank_rows}

Negative margins indicate failure. Y218F is 23 pLDDT points below and 4.2 Å above the reliability limits. R244A and L265W are each 27 pLDDT points below and 7.8 Å above the limits. These are not borderline failures: every pocket candidate inherits two confidence penalties.

## Deterministic ranking trace

1. All pocket candidates precede the non-pocket S301A.
2. Y218F precedes R244A and L265W because its joined interval has higher pLDDT (47 versus 43), then lower PAE (14.2 versus 17.8 Å).
3. R244A and L265W tie on pLDDT and PAE; the lower absolute predicted ΔΔG hypothesis ranks R244A before L265W (2.8 versus 4.1 kcal/mol).

## Evidence-quality assessment

**What is supported:** the supplied interval join, pLDDT/PAE thresholds, deterministic ordering, and conclusion that the pocket fails the prescribed reliability rule.

**What is not available:** no PDB coordinates, UniProt identifier, sequence, ligand, pocket-prediction result, or atom-level quality data were supplied. Therefore AlphaFold download, fpocket/P2Rank, structural quality, visualization, and binding-site characterization tools were not run; manufacturing such results would exceed the evidence.

**Claim boundary:** predicted ΔΔG and predicted activity are hypotheses, not measured stability, binding, or activity effects. Their values cannot overcome the low-confidence structure. The proportionate recommendation is to validate the pocket structure or obtain a better-supported model before prioritizing mutagenesis, followed by experiments to establish actual effects.
"""
(OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")

print({"order": [item.mutation for item in ranked], "pocket_reliable": pocket_reliable})
