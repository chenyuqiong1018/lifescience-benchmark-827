#!/usr/bin/env python3
"""Analyze retinal cell composition under the supplied frozen marker rule."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path

import numpy as np


SCALE_FACTOR = 10_000.0
ELIGIBILITY_FRACTION = 0.01


def load_panel(path: Path) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["cell_type", "markers"]:
            raise ValueError("MARKER_PANEL.tsv must have cell_type and markers columns")
        panel = [
            {"cell_type": row["cell_type"], "markers": [x.strip() for x in row["markers"].split(",")]}
            for row in reader
        ]
    if not panel or any(not row["markers"] for row in panel):
        raise ValueError("The marker panel must contain nonempty marker lists")
    return panel


def load_genes(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames is None or "gene_symbols" not in reader.fieldnames:
            raise ValueError("Gene file lacks gene_symbols")
        return [row["gene_symbols"] for row in reader]


def read_dimensions(stream) -> tuple[int, int, int]:
    if stream.readline().strip() != "%%MatrixMarket matrix coordinate integer general":
        raise ValueError("Expected an integer, coordinate, general Matrix Market file")
    line = stream.readline()
    while line.startswith("%"):
        line = stream.readline()
    return tuple(int(value) for value in line.split())


def extract_matrix_statistics(
    matrix_path: Path, genes: list[str], marker_names: list[str]
) -> tuple[np.ndarray, np.ndarray, dict[str, int | float]]:
    marker_to_position = {name: i for i, name in enumerate(marker_names)}
    marker_gene_rows: dict[int, int] = {}
    occurrences = {name: 0 for name in marker_names}
    for gene_row, symbol in enumerate(genes, start=1):
        position = marker_to_position.get(symbol)
        if position is not None:
            occurrences[symbol] += 1
            marker_gene_rows[gene_row] = position
    bad = {name: count for name, count in occurrences.items() if count != 1}
    if bad:
        raise ValueError(f"Every marker must occur exactly once: {bad}")

    with gzip.open(matrix_path, "rt", encoding="ascii", newline="") as stream:
        n_genes, n_cells, expected_entries = read_dimensions(stream)
        if n_genes != len(genes):
            raise ValueError("Matrix gene dimension and gene table disagree")
        libraries = np.zeros(n_cells, dtype=np.float64)
        marker_counts = np.zeros((len(marker_names), n_cells), dtype=np.float64)
        parsed = 0
        for line in stream:
            if not line.strip() or line.startswith("%"):
                continue
            gene_row, cell_column, count = map(int, line.split())
            if gene_row < 1 or gene_row > n_genes or cell_column < 1 or cell_column > n_cells:
                raise ValueError("Stored coordinate is outside matrix dimensions")
            if count < 0:
                raise ValueError("Negative counts are invalid")
            cell = cell_column - 1
            libraries[cell] += count
            marker_position = marker_gene_rows.get(gene_row)
            if marker_position is not None:
                marker_counts[marker_position, cell] += count
            parsed += 1
    if parsed != expected_entries:
        raise ValueError(f"Header says {expected_entries} entries, parsed {parsed}")
    if np.any(libraries == 0):
        raise ValueError("Zero-library cells cannot be normalized")
    qc = {
        "n_genes": n_genes,
        "n_cells": n_cells,
        "n_nonzero_entries": expected_entries,
        "empty_cells": int(np.count_nonzero(libraries == 0)),
        "library_min": int(libraries.min()),
        "library_median": float(np.median(libraries)),
        "library_max": int(libraries.max()),
    }
    return libraries, marker_counts, qc


def assign_types(
    libraries: np.ndarray,
    marker_counts: np.ndarray,
    marker_names: list[str],
    panel: list[dict[str, object]],
) -> tuple[np.ndarray, dict[str, int | float]]:
    marker_index = {name: i for i, name in enumerate(marker_names)}
    log_normalized = np.log1p(marker_counts * (SCALE_FACTOR / libraries)[None, :])
    score_rows = []
    for row in panel:
        positions = [marker_index[name] for name in row["markers"]]
        score_rows.append(log_normalized[positions].mean(axis=0))
    scores = np.vstack(score_rows)
    assignments = np.argmax(scores, axis=0)  # np.argmax returns first row in a tie
    ordered = np.sort(scores, axis=0)
    margin = ordered[-1] - ordered[-2]
    return assignments, {
        "exact_top_ties": int(np.sum(margin == 0)),
        "median_top_margin": float(np.median(margin)),
        "fraction_top_margin_below_0_05": float(np.mean(margin < 0.05)),
    }


def make_composition(
    sample: str, assignments: np.ndarray, panel: list[dict[str, object]]
) -> list[dict[str, object]]:
    counts = np.bincount(assignments, minlength=len(panel))
    return [
        {
            "sample": sample,
            "cell_type": row["cell_type"],
            "n_cells": int(counts[i]),
            "fraction": float(counts[i] / assignments.size),
        }
        for i, row in enumerate(panel)
    ]


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    p = successes / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return center - half, center + half


def approximate_ratio_interval(x1: int, n1: int, x2: int, n2: int) -> tuple[float, float]:
    ratio = (x2 / n2) / (x1 / n1)
    se = math.sqrt(1 / x2 - 1 / n2 + 1 / x1 - 1 / n1)
    return math.exp(math.log(ratio) - 1.959963984540054 * se), math.exp(math.log(ratio) + 1.959963984540054 * se)


def write_report(
    path: Path,
    call: dict[str, object],
    qc1: dict[str, int | float],
    qc2: dict[str, int | float],
    uncertainty1: dict[str, int | float],
    uncertainty2: dict[str, int | float],
) -> None:
    x1, n1 = int(call["sample_1_n_cells"]), int(qc1["n_cells"])
    x2, n2 = int(call["sample_2_n_cells"]), int(qc2["n_cells"])
    p1, p2 = float(call["sample_1_fraction"]), float(call["sample_2_fraction"])
    ratio = float(call["sample_2_to_sample_1_fraction_ratio"])
    p1_ci, p2_ci = wilson_interval(x1, n1), wilson_interval(x2, n2)
    ratio_ci = approximate_ratio_interval(x1, n1, x2, n2)
    text = f"""# Differential retinal cell composition

## Reproducible method

The matrices were used exactly as supplied: {qc1['n_genes']:,} genes x {n1:,} cells in sample 1 and {qc2['n_genes']:,} genes x {n2:,} cells in sample 2. No cells or genes were filtered. Each cell was divided by its full library size, scaled to 10,000, and transformed with `log1p`. For each marker-panel cell type, the arithmetic mean across its listed transformed markers was computed. The maximum score defined the label, and panel row order broke ties. Fractions use all matrix columns.

The frozen depleted-call rule was then applied without post hoc changes: retain types with sample-1 fraction >= 0.01 and choose the smallest sample-2/sample-1 fraction ratio.

## QC

- Sample 1: {qc1['n_nonzero_entries']:,} stored entries; library size {qc1['library_min']:,}–{qc1['library_max']:,}, median {qc1['library_median']:.1f}; {qc1['empty_cells']} empty cells.
- Sample 2: {qc2['n_nonzero_entries']:,} stored entries; library size {qc2['library_min']:,}–{qc2['library_max']:,}, median {qc2['library_median']:.1f}; {qc2['empty_cells']} empty cells.
- Every listed marker appeared exactly once, and parsed stored-entry counts matched both headers.

## Result

The depleted call is **{call['cell_type']}**: {x1:,}/{n1:,} cells in sample 1 (fraction {p1:.6f}) versus {x2:,}/{n2:,} in sample 2 (fraction {p2:.6f}). The sample-2/sample-1 fraction ratio is {ratio:.6f}, meaning the sample-2 fraction is about {100*ratio:.1f}% of the sample-1 fraction under the frozen annotation.

For scale only, cell-level binomial approximations give 95% Wilson intervals [{p1_ci[0]:.6f}, {p1_ci[1]:.6f}] and [{p2_ci[0]:.6f}, {p2_ci[1]:.6f}], and a log-ratio approximation [{ratio_ci[0]:.6f}, {ratio_ci[1]:.6f}]. These are descriptive conditional calculations, not valid sample-level biological inference when cells share donors or preparations.

## Annotation evidence and uncertainty

Labels are supported only by the supplied marker panels after the fixed transformation. Sample 1 has {uncertainty1['exact_top_ties']} exact top-score ties and a median winning-score margin of {uncertainty1['median_top_margin']:.6f}; sample 2 has {uncertainty2['exact_top_ties']} ties and median margin {uncertainty2['median_top_margin']:.6f}. Required tie handling was retained.

The marker rule is deterministic but biologically limited. It does not account for doublets, ambient RNA, marker overlap, continuous states, batch effects, or donor-level replication. With one aggregate matrix per sample, uncertainty between biological replicates and causality cannot be established. The call should therefore be interpreted as composition depletion under the supplied scoring rule, not proof that the population was biologically eliminated.
"""
    path.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample1", type=Path, required=True)
    parser.add_argument("--sample2", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--markers", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    panel = load_panel(args.markers)
    genes = load_genes(args.genes)
    markers = list(dict.fromkeys(name for row in panel for name in row["markers"]))
    library1, counts1, qc1 = extract_matrix_statistics(args.sample1, genes, markers)
    library2, counts2, qc2 = extract_matrix_statistics(args.sample2, genes, markers)
    labels1, uncertainty1 = assign_types(library1, counts1, markers, panel)
    labels2, uncertainty2 = assign_types(library2, counts2, markers, panel)
    rows1 = make_composition("sample_1", labels1, panel)
    rows2 = make_composition("sample_2", labels2, panel)

    with (args.output_dir / "composition.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, ["sample", "cell_type", "n_cells", "fraction"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows1 + rows2)

    sample1 = {row["cell_type"]: row for row in rows1}
    sample2 = {row["cell_type"]: row for row in rows2}
    eligible = [row["cell_type"] for row in rows1 if row["fraction"] >= ELIGIBILITY_FRACTION]
    called = min(eligible, key=lambda name: sample2[name]["fraction"] / sample1[name]["fraction"])
    call = {
        "cell_type": called,
        "sample_1_n_cells": sample1[called]["n_cells"],
        "sample_1_fraction": sample1[called]["fraction"],
        "sample_2_n_cells": sample2[called]["n_cells"],
        "sample_2_fraction": sample2[called]["fraction"],
        "sample_2_to_sample_1_fraction_ratio": sample2[called]["fraction"] / sample1[called]["fraction"],
        "sample_1_eligibility_threshold": ELIGIBILITY_FRACTION,
        "eligible_cell_types": eligible,
    }
    (args.output_dir / "depleted_call.json").write_text(json.dumps(call, indent=2) + "\n", encoding="utf-8")
    write_report(args.output_dir / "report.md", call, qc1, qc2, uncertainty1, uncertainty2)


if __name__ == "__main__":
    main()
