from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


OUTPUT_DIR = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[5]
INPUT_DIR = REPO / "inputs" / "ls05-structure-model-ranking"
CRITICAL = range(181, 241)

MODEL_COLUMNS = {
    "model_id",
    "tm_score",
    "lddt",
    "rmsd_a",
    "aligned_residues",
    "chain_mapping_complete",
}
ERROR_COLUMNS = {"model_id", "residue_start", "residue_end", "mean_error_a"}


@dataclass(frozen=True)
class Model:
    model_id: str
    tm_score: float
    lddt: float
    rmsd_a: float
    aligned_residues: int
    mapping_complete: bool
    critical_risk: float

    @property
    def sort_key(self) -> tuple[object, ...]:
        return (
            not self.mapping_complete,
            -self.tm_score,
            -self.lddt,
            self.rmsd_a,
            -self.aligned_residues,
            self.critical_risk,
            self.model_id,
        )


def read_csv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if set(reader.fieldnames or []) != required:
            raise ValueError(f"Unexpected columns in {path.name}: {reader.fieldnames}")
        rows = list(reader)
    if not rows or any(value == "" for row in rows for value in row.values()):
        raise ValueError(f"Empty or missing values in {path.name}")
    return rows


def risk_for_model(error_rows: list[dict[str, str]], model_id: str) -> float:
    # Accumulate by residue so overlap length is explicit and coverage gaps are caught.
    weighted: dict[int, list[float]] = {residue: [] for residue in CRITICAL}
    for row in error_rows:
        if row["model_id"] != model_id:
            continue
        start, end = int(row["residue_start"]), int(row["residue_end"])
        if start > end:
            raise ValueError(f"Reversed interval for {model_id}: {start}-{end}")
        for residue in CRITICAL:
            if start <= residue <= end:
                weighted[residue].append(float(row["mean_error_a"]))
    uncovered = [residue for residue, values in weighted.items() if not values]
    duplicated = [residue for residue, values in weighted.items() if len(values) > 1]
    if uncovered or duplicated:
        raise ValueError(
            f"Critical-region coverage invalid for {model_id}; "
            f"uncovered={uncovered}, duplicated={duplicated}"
        )
    return sum(values[0] for values in weighted.values()) / len(weighted)


metric_rows = read_csv(INPUT_DIR / "model_metrics.csv", MODEL_COLUMNS)
error_rows = read_csv(INPUT_DIR / "residue_errors.csv", ERROR_COLUMNS)

models: list[Model] = []
for row in metric_rows:
    mapping_text = row["chain_mapping_complete"].lower()
    if mapping_text not in {"true", "false"}:
        raise ValueError(f"Invalid mapping flag for {row['model_id']}")
    models.append(
        Model(
            model_id=row["model_id"],
            tm_score=float(row["tm_score"]),
            lddt=float(row["lddt"]),
            rmsd_a=float(row["rmsd_a"]),
            aligned_residues=int(row["aligned_residues"]),
            mapping_complete=mapping_text == "true",
            critical_risk=risk_for_model(error_rows, row["model_id"]),
        )
    )

if len({model.model_id for model in models}) != len(models):
    raise ValueError("Duplicate model_id")

ranked = sorted(models, key=lambda model: model.sort_key)
records: list[dict[str, str | int | float]] = []
for rank, model in enumerate(ranked, start=1):
    decision = (
        "preferred"
        if rank == 1
        else "alternate"
        if model.mapping_complete
        else "reject_incomplete_mapping"
    )
    records.append(
        {
            "rank": rank,
            "model_id": model.model_id,
            "global_score": f"{model.tm_score:.2f}",
            "interface_score": int(model.mapping_complete),
            "critical_residue_risk": f"{model.critical_risk:.1f}",
            "decision": decision,
        }
    )

headers = [
    "rank",
    "model_id",
    "global_score",
    "interface_score",
    "critical_residue_risk",
    "decision",
]
with (OUTPUT_DIR / "model_ranking.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=headers)
    writer.writeheader()
    writer.writerows(records)

result_rows = "\n".join(
    f"| {row['rank']} | {row['model_id']} | {row['global_score']} | "
    f"{row['interface_score']} | {row['critical_residue_risk']} | {row['decision']} |"
    for row in records
)
diagnostic_rows = "\n".join(
    f"| {model.model_id} | {'complete' if model.mapping_complete else 'incomplete'} | "
    f"{model.tm_score:.2f} | {model.lddt:.2f} | {model.rmsd_a:.1f} | "
    f"{model.aligned_residues} | {model.critical_risk:.1f} |"
    for model in ranked
)
report = f"""# Structural model ranking

## Frozen-rule result

| Rank | Model | Global score | Interface score | Critical-residue risk (Å) | Decision |
|---:|---|---:|---:|---:|---|
{result_rows}

## Input and rule audit

| Model | Chain mapping | TM-score | lDDT | RMSD (Å) | Aligned residues | Residues 181–240 risk (Å) |
|---|---|---:|---:|---:|---:|---:|
{diagnostic_rows}

The analysis found all required columns, no missing values, unique model identifiers, and exactly one error value covering every residue in the critical interval 181–240 for each model. Because the supplied rule is deterministic, no hypothesis test, p-value, or fitted composite score is warranted for these three models.

## Ranking interpretation

Chain-mapping completeness is the first and therefore decisive discriminator. It places model_A and model_C ahead of model_B, even though model_B's TM-score (0.74) is higher than model_C's (0.62). The prescribed `interface_score` records only this completeness state; it is not evidence of physical interface quality.

Among complete mappings, model_A is preferred because its higher TM-score (0.81 versus 0.62) is encountered before all remaining tie-breakers. Critical-region uncertainty is evaluated over residues 181–240. Model_C's risk is lower than model_A's (4.1 versus 4.8 Å), but risk is sixth in the tuple and cannot reverse the earlier TM-score ordering. Model_B combines incomplete mapping with the highest risk (8.5 Å) and is rejected for incomplete mapping.

Only supplied comparison metrics are interpreted. No coordinate-level geometry, physical interface, or experimental property is claimed.
"""
(OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")

print(" > ".join(model.model_id for model in ranked))
