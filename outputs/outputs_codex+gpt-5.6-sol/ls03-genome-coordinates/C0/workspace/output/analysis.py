#!/usr/bin/env python3
"""Analyze enhancer-promoter distance, contact, and lagged transcription."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path


CONTACT_THRESHOLD_NM = 260.0
EXPECTED_TIMEPOINTS = 250
MIN_LAG = -20
MAX_LAG = 20


def load_data(path: Path):
    cells = defaultdict(list)
    required = [
        "cell_id", "time", "enh_x", "enh_y", "enh_z",
        "prom_x", "prom_y", "prom_z", "transcription",
    ]
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != required:
            raise ValueError(f"Unexpected columns: {reader.fieldnames}")
        for row_number, row in enumerate(reader, 2):
            cell_id = int(row["cell_id"])
            time = int(row["time"])
            enh = tuple(float(row[name]) for name in ("enh_x", "enh_y", "enh_z"))
            prom = tuple(float(row[name]) for name in ("prom_x", "prom_y", "prom_z"))
            transcription = int(row["transcription"])
            values = (*enh, *prom)
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"Non-finite coordinate on row {row_number}")
            if transcription not in (0, 1):
                raise ValueError(f"Non-binary transcription on row {row_number}")
            distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(enh, prom)))
            cells[cell_id].append((time, distance, int(distance <= CONTACT_THRESHOLD_NM), transcription))

    for cell_id, rows in cells.items():
        rows.sort(key=lambda item: item[0])
        times = [item[0] for item in rows]
        if len(rows) != EXPECTED_TIMEPOINTS or times != list(range(EXPECTED_TIMEPOINTS)):
            raise ValueError(f"Cell {cell_id} does not contain exactly times 0..249")
    if not cells:
        raise ValueError("No cells found")
    return dict(sorted(cells.items()))


def pearson_binary(pairs):
    n = len(pairs)
    if n < 2:
        return float("nan")
    sx = sum(x for x, _ in pairs)
    sy = sum(y for _, y in pairs)
    sxx = sum(x * x for x, _ in pairs)
    syy = sum(y * y for _, y in pairs)
    sxy = sum(x * y for x, y in pairs)
    numerator = n * sxy - sx * sy
    denominator = math.sqrt((n * sxx - sx * sx) * (n * syy - sy * sy))
    return numerator / denominator if denominator else float("nan")


def cell_metrics(cells):
    results = []
    for cell_id, rows in cells.items():
        n = len(rows)
        results.append({
            "cell_id": cell_id,
            "n_timepoints": n,
            "mean_distance_nm": sum(row[1] for row in rows) / n,
            "contact_fraction": sum(row[2] for row in rows) / n,
            "transcription_fraction": sum(row[3] for row in rows) / n,
        })
    return results


def lag_metrics(cells):
    results = []
    for lag in range(MIN_LAG, MAX_LAG + 1):
        pairs = []
        for rows in cells.values():
            if lag >= 0:
                indices = range(0, len(rows) - lag)
            else:
                indices = range(-lag, len(rows))
            pairs.extend((rows[t][2], rows[t + lag][3]) for t in indices)
        association = pearson_binary(pairs)
        if not math.isfinite(association):
            raise ValueError(f"Undefined association at lag {lag}")
        results.append({"lag": lag, "association": association, "n_observations": len(pairs)})
    return results


def write_outputs(out_dir: Path, cells, per_cell, lagged):
    with (out_dir / "cell_metrics.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "cell_id", "n_timepoints", "mean_distance_nm",
                "contact_fraction", "transcription_fraction",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        for row in per_cell:
            writer.writerow({
                **row,
                "mean_distance_nm": f'{row["mean_distance_nm"]:.12f}',
                "contact_fraction": f'{row["contact_fraction"]:.12f}',
                "transcription_fraction": f'{row["transcription_fraction"]:.12f}',
            })

    with (out_dir / "lag_analysis.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["lag", "association", "n_observations"],
            lineterminator="\n",
        )
        writer.writeheader()
        for row in lagged:
            writer.writerow({**row, "association": f'{row["association"]:.12f}'})

    strongest = max(lagged, key=lambda row: (abs(row["association"]), -abs(row["lag"]), -row["lag"]))
    zero = next(row for row in lagged if row["lag"] == 0)
    mean_contact = sum(row["contact_fraction"] for row in per_cell) / len(per_cell)
    mean_tx = sum(row["transcription_fraction"] for row in per_cell) / len(per_cell)
    mean_distance = sum(row["mean_distance_nm"] for row in per_cell) / len(per_cell)
    direction = (
        "contact leads the later transcription measurement"
        if strongest["lag"] > 0 else
        "transcription leads the later contact measurement"
        if strongest["lag"] < 0 else
        "contact and transcription are contemporaneous"
    )
    report = f"""# Enhancer-promoter distance and transcription dynamics

## Summary

The dataset contains {len(cells)} cells and {sum(len(rows) for rows in cells.values()):,}
observations, with exactly {EXPECTED_TIMEPOINTS} timepoints per cell. Across cells,
the mean of cell-level mean enhancer-promoter distances is {mean_distance:.3f} nm,
the mean contact fraction at the supplied `≤{CONTACT_THRESHOLD_NM:.0f} nm` threshold is
{mean_contact:.4f}, and the mean transcription fraction is {mean_tx:.4f}.

The strongest absolute pooled lag association occurs at **lag {strongest['lag']:+d}**:
Pearson `r = {strongest['association']:.6f}` from
{strongest['n_observations']:,} within-cell pairs. Under the supplied sign convention,
{direction}. At lag 0, `r = {zero['association']:.6f}` from
{zero['n_observations']:,} pairs.

## Method

For every row, enhancer-promoter distance was calculated as the three-dimensional
Euclidean distance in nanometres. Contact equals 1 when distance is at most 260 nm.
Cell summaries use all 250 rows for that cell. For each integer lag from -20 through
+20, contact at time `t` was paired with transcription at `t+lag` only when both
indices existed within the same cell. Pearson correlation was then calculated after
pooling those valid within-cell pairs. No pairs cross cell boundaries.

## Interpretation and limitations

The lag peak is a temporal association, not evidence that physical contact causes
transcription. These observational trajectories cannot establish causal direction,
exclude a shared upstream regulator, remove time-dependent or cell-state confounding,
or show that experimentally changing distance would change transcription. Binary
contact also discards variation above and below the 260-nm threshold. A causal claim
would require an intervention (for example targeted perturbation of contact), suitable
controls, and a design that addresses temporal autocorrelation and cell-level effects.

## Reproduction

Run `python analysis.py` from any directory. The script locates the repository from
its own path, validates the 250-point within-cell time grids, and deterministically
rewrites both CSV files and this report.
"""
    with (out_dir / "report.md").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(report)

    return strongest


def main():
    script = Path(__file__).resolve()
    repo = next(parent for parent in script.parents if parent.name == "lifescience-benchmark-827")
    input_path = repo / "inputs" / "ls03-genome-coordinates" / "single_cell_dynamics_question.csv"
    cells = load_data(input_path)
    per_cell = cell_metrics(cells)
    lagged = lag_metrics(cells)
    strongest = write_outputs(script.parent, cells, per_cell, lagged)
    print({
        "cells": len(cells),
        "rows": sum(len(rows) for rows in cells.values()),
        "strongest_lag": strongest["lag"],
        "association": strongest["association"],
        "n_observations": strongest["n_observations"],
    })


if __name__ == "__main__":
    main()
