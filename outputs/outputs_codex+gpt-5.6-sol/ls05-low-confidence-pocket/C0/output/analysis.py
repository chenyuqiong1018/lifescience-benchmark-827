from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path


OUTPUT_DIR = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[5]
INPUT_DIR = REPO / "inputs" / "ls05-low-confidence-pocket"


@dataclass(frozen=True)
class JoinedCandidate:
    mutation: str
    residue: int
    ddg: float
    activity_hypothesis: str
    region: str
    plddt: float
    pae: float
    pocket_member: bool

    @property
    def sort_key(self) -> tuple[object, ...]:
        if self.region == "pocket":
            return (0, -self.plddt, self.pae, abs(self.ddg), self.mutation)
        return (1, self.mutation)


def residue_number(mutation: str) -> int:
    match = re.fullmatch(r"[A-Za-z](\d+)[A-Za-z]", mutation)
    if not match:
        raise ValueError(f"Unsupported mutation format: {mutation}")
    return int(match.group(1))


with (INPUT_DIR / "confidence.csv").open(newline="", encoding="utf-8") as handle:
    confidence_rows = list(csv.DictReader(handle))
with (INPUT_DIR / "mutation_candidates.csv").open(newline="", encoding="utf-8") as handle:
    mutation_rows = list(csv.DictReader(handle))

joined: list[JoinedCandidate] = []
for row in mutation_rows:
    residue = residue_number(row["mutation"])
    matches = [
        interval
        for interval in confidence_rows
        if int(interval["residue_start"]) <= residue <= int(interval["residue_end"])
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one confidence interval for {row['mutation']}, found {len(matches)}")
    interval = matches[0]
    joined.append(
        JoinedCandidate(
            mutation=row["mutation"],
            residue=residue,
            ddg=float(row["predicted_ddg_kcal_mol"]),
            activity_hypothesis=row["predicted_activity_change"],
            region=row["region"],
            plddt=float(interval["plddt"]),
            pae=float(interval["pae_to_core_a"]),
            pocket_member=interval["pocket_member"].lower() == "true",
        )
    )

ranked = sorted(joined, key=lambda candidate: candidate.sort_key)
records: list[dict[str, object]] = []
for rank, candidate in enumerate(ranked, start=1):
    flags: list[str] = []
    if candidate.plddt < 50:
        flags.append("plddt_lt_50")
    if candidate.pae > 10:
        flags.append("pae_gt_10A")
    records.append(
        {
            "rank": rank,
            "mutation": candidate.mutation,
            "pocket_support": "unsupported_low_confidence" if candidate.plddt < 50 else "cautious_support",
            "confidence_penalty": ";".join(flags) if flags else "none",
            "decision": (
                "defer_structure_validation"
                if candidate.region == "pocket" and candidate.plddt < 50
                else "out_of_scope_non_pocket"
                if candidate.region != "pocket"
                else "cautious_support"
            ),
        }
    )

columns = ["rank", "mutation", "pocket_support", "confidence_penalty", "decision"]
with (OUTPUT_DIR / "mutation_priorities.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns)
    writer.writeheader()
    writer.writerows(records)

assessment = {
    "pocket_reliable": False,
    "prioritize_mutations": False,
    "reason": (
        "Both nominated pocket intervals are below the reliability threshold: residues 210-230 have "
        "pLDDT 47 and PAE 14.2 A, and residues 231-270 have pLDDT 43 and PAE 17.8 A. "
        "All pocket mutations should be deferred pending structure validation."
    ),
}
(OUTPUT_DIR / "pocket_assessment.json").write_text(
    json.dumps(assessment, indent=2) + "\n", encoding="utf-8"
)

priority_rows = "\n".join(
    f"| {rank} | {candidate.mutation} | {candidate.region} | {candidate.plddt:.0f} | "
    f"{candidate.pae:.1f} | {candidate.ddg:.1f} | {candidate.activity_hypothesis} | "
    f"{records[rank - 1]['decision']} |"
    for rank, candidate in enumerate(ranked, start=1)
)
report = f"""# Low-confidence pocket assessment

## Conclusion

The nominated pocket is **not reliable enough to prioritize mutations**. Both pocket intervals fail the required pLDDT ≥ 70 and PAE ≤ 10 Å criteria. Every pocket candidate therefore receives `unsupported_low_confidence`, both confidence penalties, and `defer_structure_validation`.

## Ranked candidates

| Rank | Mutation | Region | Joined pLDDT | PAE to core (Å) | Predicted ΔΔG hypothesis (kcal/mol) | Predicted activity hypothesis | Decision |
|---:|---|---|---:|---:|---:|---|---|
{priority_rows}

Pocket candidates are ordered exactly by higher joined pLDDT, then lower PAE, then lower absolute predicted ΔΔG, then mutation name. Y218F ranks first because its interval has pLDDT 47 and PAE 14.2 Å, compared with pLDDT 43 and PAE 17.8 Å for R244A and L265W. R244A precedes L265W because their confidence values tie and |2.8| < |4.1| kcal/mol. S301A is non-pocket and is placed after all pocket candidates.

## Uncertainty propagation

pLDDT below 50 means the pocket coordinates should not be interpreted, while PAE above 10 Å indicates uncertain placement relative to the core. The simultaneous failures apply to all three pocket candidates and prevent mutation prioritization from this structure hypothesis. The supplied ΔΔG and activity fields are model hypotheses used only as prescribed tie-breaker/context; they are not measured binding, stability, or activity effects.
"""
(OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")

print(" > ".join(candidate.mutation for candidate in ranked))
