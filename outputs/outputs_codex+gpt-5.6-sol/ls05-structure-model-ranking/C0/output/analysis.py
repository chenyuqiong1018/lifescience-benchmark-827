from __future__ import annotations

import csv
from pathlib import Path


OUTPUT_DIR = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[5]
INPUT_DIR = REPO / "inputs" / "ls05-structure-model-ranking"
CRITICAL_START = 181
CRITICAL_END = 240


def as_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def critical_risk(rows: list[dict[str, str]], model_id: str) -> float:
    weighted_error = 0.0
    overlap_length = 0
    for row in rows:
        if row["model_id"] != model_id:
            continue
        start = max(int(row["residue_start"]), CRITICAL_START)
        end = min(int(row["residue_end"]), CRITICAL_END)
        if start <= end:
            length = end - start + 1
            weighted_error += length * float(row["mean_error_a"])
            overlap_length += length
    if overlap_length == 0:
        raise ValueError(f"No critical-region coverage for {model_id}")
    return weighted_error / overlap_length


with (INPUT_DIR / "model_metrics.csv").open(newline="", encoding="utf-8") as handle:
    models = list(csv.DictReader(handle))
with (INPUT_DIR / "residue_errors.csv").open(newline="", encoding="utf-8") as handle:
    errors = list(csv.DictReader(handle))

for model in models:
    model["mapping_complete"] = as_bool(model["chain_mapping_complete"])
    model["risk"] = critical_risk(errors, model["model_id"])

# Frozen rule, in priority order: complete mapping; higher TM-score; higher lDDT;
# lower RMSD; more aligned residues; lower critical-region error; lexical ID.
models.sort(
    key=lambda row: (
        not row["mapping_complete"],
        -float(row["tm_score"]),
        -float(row["lddt"]),
        float(row["rmsd_a"]),
        -int(row["aligned_residues"]),
        float(row["risk"]),
        row["model_id"],
    )
)

ranked: list[dict[str, str | int | float]] = []
for rank, model in enumerate(models, start=1):
    if rank == 1:
        decision = "preferred"
    elif model["mapping_complete"]:
        decision = "alternate"
    else:
        decision = "reject_incomplete_mapping"
    ranked.append(
        {
            "rank": rank,
            "model_id": model["model_id"],
            "global_score": model["tm_score"],
            "interface_score": 1 if model["mapping_complete"] else 0,
            "critical_residue_risk": f'{float(model["risk"]):.1f}',
            "decision": decision,
        }
    )

fieldnames = [
    "rank",
    "model_id",
    "global_score",
    "interface_score",
    "critical_residue_risk",
    "decision",
]
with (OUTPUT_DIR / "model_ranking.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(ranked)

table_rows = "\n".join(
    f'| {row["rank"]} | {row["model_id"]} | {row["global_score"]} | '
    f'{row["interface_score"]} | {row["critical_residue_risk"]} | {row["decision"]} |'
    for row in ranked
)
report = f"""# Structural model ranking

## Result

| Rank | Model | Global score | Interface score | Critical-residue risk (Å) | Decision |
|---:|---|---:|---:|---:|---|
{table_rows}

## Interpretation

The ranking applies the supplied stable tuple exactly. Complete chain mapping is the first discriminator, so model_A and model_C rank ahead of model_B even though model_B has a higher TM-score than model_C. The interface score is only the prescribed mapping-completeness indicator: 1 for a complete mapping and 0 for an incomplete mapping. It is not a coordinate-derived interface-quality claim.

Within the mapping-complete group, model_A ranks first because its TM-score (0.81) exceeds model_C's (0.62). Model_C has lower critical-region uncertainty (4.1 Å versus 4.8 Å for model_A), but critical-region risk is the sixth tie-breaker and therefore does not override the earlier TM-score comparison. Model_B has both incomplete mapping and the highest critical-region uncertainty (8.5 Å), so its decision is `reject_incomplete_mapping`.

Critical-region risk is the length-weighted mean error over residues 181–240. Each supplied model has one row spanning that full interval, so the result equals that row's mean error. These values are already-computed comparison metrics; no coordinate-level, physical-interface, or experimental property is inferred.
"""
(OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")

print(f"Wrote {len(ranked)} ranked models to {OUTPUT_DIR}")
