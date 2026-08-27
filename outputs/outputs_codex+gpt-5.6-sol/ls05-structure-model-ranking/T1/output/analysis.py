from __future__ import annotations

import csv
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path


OUTPUT_DIR = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[5]
INPUT_DIR = REPO / "inputs" / "ls05-structure-model-ranking"
CRITICAL_START, CRITICAL_END = 181, 240


@dataclass(frozen=True)
class Candidate:
    model_id: str
    tm_score: Decimal
    lddt: Decimal
    rmsd_a: Decimal
    aligned_residues: int
    mapping_complete: bool
    critical_risk: Decimal
    critical_overlap: int

    def ranking_tuple(self) -> tuple[object, ...]:
        return (
            not self.mapping_complete,
            -self.tm_score,
            -self.lddt,
            self.rmsd_a,
            -self.aligned_residues,
            self.critical_risk,
            self.model_id,
        )


def load_rows(path: Path, fields: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != fields:
            raise ValueError(f"Schema mismatch in {path.name}: {reader.fieldnames}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"No rows in {path.name}")
    missing = [(index + 2, key) for index, row in enumerate(rows) for key, value in row.items() if value == ""]
    if missing:
        raise ValueError(f"Missing values in {path.name}: {missing}")
    return rows


def weighted_critical_risk(error_rows: list[dict[str, str]], model_id: str) -> tuple[Decimal, int]:
    weighted_sum = Decimal("0")
    overlap_total = 0
    for row in error_rows:
        if row["model_id"] != model_id:
            continue
        row_start = int(row["residue_start"])
        row_end = int(row["residue_end"])
        if row_start > row_end:
            raise ValueError(f"Reversed interval for {model_id}")
        overlap_start = max(row_start, CRITICAL_START)
        overlap_end = min(row_end, CRITICAL_END)
        if overlap_start <= overlap_end:
            overlap = overlap_end - overlap_start + 1
            weighted_sum += Decimal(row["mean_error_a"]) * overlap
            overlap_total += overlap
    if overlap_total == 0:
        raise ValueError(f"Critical region is uncovered for {model_id}")
    return weighted_sum / overlap_total, overlap_total


metric_fields = (
    "model_id",
    "tm_score",
    "lddt",
    "rmsd_a",
    "aligned_residues",
    "chain_mapping_complete",
)
error_fields = ("model_id", "residue_start", "residue_end", "mean_error_a")
metric_rows = load_rows(INPUT_DIR / "model_metrics.csv", metric_fields)
error_rows = load_rows(INPUT_DIR / "residue_errors.csv", error_fields)

candidates: list[Candidate] = []
for row in metric_rows:
    mapping = row["chain_mapping_complete"].lower()
    if mapping not in {"true", "false"}:
        raise ValueError(f"Invalid mapping flag for {row['model_id']}")
    risk, overlap = weighted_critical_risk(error_rows, row["model_id"])
    candidates.append(
        Candidate(
            model_id=row["model_id"],
            tm_score=Decimal(row["tm_score"]),
            lddt=Decimal(row["lddt"]),
            rmsd_a=Decimal(row["rmsd_a"]),
            aligned_residues=int(row["aligned_residues"]),
            mapping_complete=mapping == "true",
            critical_risk=risk,
            critical_overlap=overlap,
        )
    )

if len({candidate.model_id for candidate in candidates}) != len(candidates):
    raise ValueError("model_id values must be unique")
if {row["model_id"] for row in error_rows} != {candidate.model_id for candidate in candidates}:
    raise ValueError("Model sets differ between input files")

ranked = sorted(candidates, key=Candidate.ranking_tuple)
records: list[dict[str, object]] = []
for rank, candidate in enumerate(ranked, 1):
    decision = (
        "preferred"
        if rank == 1
        else "alternate"
        if candidate.mapping_complete
        else "reject_incomplete_mapping"
    )
    records.append(
        {
            "rank": rank,
            "model_id": candidate.model_id,
            "global_score": str(candidate.tm_score),
            "interface_score": int(candidate.mapping_complete),
            "critical_residue_risk": f"{candidate.critical_risk:.1f}",
            "decision": decision,
        }
    )

columns = [
    "rank",
    "model_id",
    "global_score",
    "interface_score",
    "critical_residue_risk",
    "decision",
]
with (OUTPUT_DIR / "model_ranking.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns)
    writer.writeheader()
    writer.writerows(records)

rank_rows = "\n".join(
    f"| {row['rank']} | {row['model_id']} | {row['global_score']} | "
    f"{row['interface_score']} | {row['critical_residue_risk']} | {row['decision']} |"
    for row in records
)
audit_rows = "\n".join(
    f"| {candidate.model_id} | {'complete' if candidate.mapping_complete else 'incomplete'} | "
    f"{candidate.tm_score} | {candidate.lddt} | {candidate.rmsd_a} | "
    f"{candidate.aligned_residues} | {candidate.critical_risk:.1f} | {candidate.critical_overlap} |"
    for candidate in ranked
)

report = f"""# Structural model ranking

## Final ranking

| Rank | Model | Global score | Interface score | Critical-residue risk (Å) | Decision |
|---:|---|---:|---:|---:|---|
{rank_rows}

## Audited input metrics

| Model | Chain mapping | TM-score | lDDT | RMSD (Å) | Aligned residues | Critical risk (Å) | Critical overlap (residues) |
|---|---|---:|---:|---:|---:|---:|---:|
{audit_rows}

Both input schemas match exactly, all fields are populated, model identifiers are unique and consistent across files, and each model contributes 60 residue positions to the length-weighted mean over residues 181–240.

## Decision trace

1. **model_A versus model_C:** both have complete chain mappings, so the first tuple element ties. TM-score is the first differing element; 0.81 exceeds 0.62, making model_A preferred. model_C's lower critical-region risk (4.1 Å versus 4.8 Å) occurs only at the sixth discriminator and cannot reverse that decision.
2. **model_C versus model_B:** chain-mapping completeness is the first discriminator. model_C is therefore ahead despite model_B's higher TM-score (0.74 versus 0.62), higher lDDT, and lower RMSD. Incomplete mapping also fixes model_B's interface score at 0 and its decision at `reject_incomplete_mapping`.

Critical-region uncertainty is highest for model_B (8.5 Å), then model_A (4.8 Å), then model_C (4.1 Å). These risks can resolve a tie only after mapping completeness, TM-score, lDDT, RMSD, and aligned-residue count all tie.

## Evidence boundary

The supplied files contain comparison metrics but no coordinates or PDB identifier. Accordingly, no coordinate geometry, atom-level quality, physical interface, composition, visualization, or experimental property was calculated or claimed. `interface_score` is strictly the rule-prescribed binary chain-mapping indicator. With a deterministic frozen tuple and only three candidates, inferential hypothesis tests, p-values, effect sizes, and fitted composite scores would be unsupported and were not introduced.
"""
(OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")

print([candidate.model_id for candidate in ranked])
