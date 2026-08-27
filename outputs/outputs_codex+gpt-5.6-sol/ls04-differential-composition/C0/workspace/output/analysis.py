#!/usr/bin/env python3
"""Frozen marker-score composition analysis for two Matrix Market samples."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np


TARGET_SUM = 10_000.0
MIN_SAMPLE1_FRACTION = 0.01


def read_marker_panel(path: Path) -> list[tuple[str, list[str]]]:
    panel: list[tuple[str, list[str]]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["cell_type", "markers"]:
            raise ValueError(f"Unexpected marker panel columns: {reader.fieldnames}")
        for row in reader:
            markers = [item.strip() for item in row["markers"].split(",") if item.strip()]
            if not markers:
                raise ValueError(f"No markers for {row['cell_type']}")
            panel.append((row["cell_type"], markers))
    if not panel:
        raise ValueError("Marker panel is empty")
    return panel


def read_gene_symbols(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "gene_symbols" not in reader.fieldnames:
            raise ValueError("Gene table must contain gene_symbols")
        return [row["gene_symbols"] for row in reader]


def parse_matrix_header(handle) -> tuple[int, int, int]:
    banner = handle.readline().strip()
    if banner != "%%MatrixMarket matrix coordinate integer general":
        raise ValueError(f"Unsupported Matrix Market banner: {banner}")
    line = handle.readline()
    while line.startswith("%"):
        line = handle.readline()
    n_genes, n_cells, n_entries = (int(item) for item in line.split())
    return n_genes, n_cells, n_entries


def load_sufficient_data(
    matrix_path: Path,
    gene_symbols: list[str],
    required_markers: list[str],
) -> tuple[np.ndarray, np.ndarray, dict[str, int | float]]:
    symbol_rows: dict[str, list[int]] = {marker: [] for marker in required_markers}
    for row_index, symbol in enumerate(gene_symbols, start=1):
        if symbol in symbol_rows:
            symbol_rows[symbol].append(row_index)
    missing = [marker for marker, rows in symbol_rows.items() if not rows]
    if missing:
        raise ValueError(f"Marker genes absent from gene table: {missing}")
    duplicates = {marker: rows for marker, rows in symbol_rows.items() if len(rows) != 1}
    if duplicates:
        raise ValueError(f"Marker gene symbols must map uniquely: {duplicates}")

    row_to_marker = {rows[0]: i for i, rows in enumerate(symbol_rows.values())}
    with gzip.open(matrix_path, "rt", encoding="ascii", newline="") as handle:
        n_genes, n_cells, n_entries = parse_matrix_header(handle)
        if n_genes != len(gene_symbols):
            raise ValueError(f"Matrix has {n_genes} genes but gene table has {len(gene_symbols)}")
        library_sizes = np.zeros(n_cells, dtype=np.float64)
        marker_counts = np.zeros((len(required_markers), n_cells), dtype=np.float64)
        seen_entries = 0
        for line in handle:
            if not line.strip() or line.startswith("%"):
                continue
            gene_index, cell_index, value = (int(item) for item in line.split())
            if not (1 <= gene_index <= n_genes and 1 <= cell_index <= n_cells):
                raise ValueError("Matrix coordinate outside declared dimensions")
            if value < 0:
                raise ValueError("Counts must be nonnegative")
            j = cell_index - 1
            library_sizes[j] += value
            marker_index = row_to_marker.get(gene_index)
            if marker_index is not None:
                marker_counts[marker_index, j] += value
            seen_entries += 1
    if seen_entries != n_entries:
        raise ValueError(f"Expected {n_entries} entries, parsed {seen_entries}")
    if np.any(library_sizes <= 0):
        raise ValueError("Every cell must have a positive library size for normalization")
    qc = {
        "n_genes": n_genes,
        "n_cells": n_cells,
        "n_nonzero_entries": n_entries,
        "empty_cells": int(np.count_nonzero(library_sizes == 0)),
        "min_library_size": float(library_sizes.min()),
        "median_library_size": float(np.median(library_sizes)),
        "max_library_size": float(library_sizes.max()),
    }
    return library_sizes, marker_counts, qc


def annotate(
    library_sizes: np.ndarray,
    marker_counts: np.ndarray,
    panel: list[tuple[str, list[str]]],
    required_markers: list[str],
) -> tuple[np.ndarray, np.ndarray, dict[str, float | int]]:
    marker_lookup = {marker: i for i, marker in enumerate(required_markers)}
    normalized = np.log1p(marker_counts / library_sizes[np.newaxis, :] * TARGET_SUM)
    scores = np.vstack(
        [normalized[[marker_lookup[m] for m in markers], :].mean(axis=0) for _, markers in panel]
    )
    labels = np.argmax(scores, axis=0)  # first maximum implements panel-order tie breaking
    ordered = np.sort(scores, axis=0)
    margins = ordered[-1] - ordered[-2]
    uncertainty = {
        "n_exact_top_score_ties": int(np.count_nonzero(margins == 0)),
        "median_top_score_margin": float(np.median(margins)),
        "fraction_margin_below_0_05": float(np.mean(margins < 0.05)),
    }
    return labels, scores, uncertainty


def composition_rows(
    sample_name: str, labels: np.ndarray, panel: list[tuple[str, list[str]]]
) -> list[dict[str, str | int | float]]:
    counts = np.bincount(labels, minlength=len(panel))
    return [
        {
            "sample": sample_name,
            "cell_type": cell_type,
            "n_cells": int(counts[i]),
            "fraction": float(counts[i] / labels.size),
        }
        for i, (cell_type, _) in enumerate(panel)
    ]


def write_report(
    path: Path,
    panel: list[tuple[str, list[str]]],
    rows: list[dict[str, str | int | float]],
    depleted: dict[str, object],
    qc: dict[str, dict[str, int | float]],
    uncertainty: dict[str, dict[str, float | int]],
) -> None:
    call = str(depleted["cell_type"])
    ratio = float(depleted["sample_2_to_sample_1_fraction_ratio"])
    f1 = float(depleted["sample_1_fraction"])
    f2 = float(depleted["sample_2_fraction"])
    n1 = int(depleted["sample_1_n_cells"])
    n2 = int(depleted["sample_2_n_cells"])
    report = f"""# Retinal differential composition

## Frozen workflow

Both input files are integer, coordinate-format Matrix Market count matrices with genes in rows and cells in columns. No cell or gene was removed. For every cell, all gene counts contributed to its library size; counts were divided by that library size, multiplied by 10,000, and transformed with `log1p`. For each of the {len(panel)} marker-panel rows, the score is the arithmetic mean of the transformed values of exactly the listed markers. Each cell was assigned to the largest score, with marker-panel row order resolving exact ties.

Fractions use the full declared matrix column count. The depleted population is selected exactly as specified: among types with sample-1 fraction >= 1%, choose the smallest sample-2/sample-1 fraction ratio.

## Quality control

- Sample 1: {qc['sample_1']['n_genes']:,} genes x {qc['sample_1']['n_cells']:,} cells; {qc['sample_1']['n_nonzero_entries']:,} stored entries; 0 empty libraries. Library sizes range from {qc['sample_1']['min_library_size']:.0f} to {qc['sample_1']['max_library_size']:.0f} counts (median {qc['sample_1']['median_library_size']:.1f}).
- Sample 2: {qc['sample_2']['n_genes']:,} genes x {qc['sample_2']['n_cells']:,} cells; {qc['sample_2']['n_nonzero_entries']:,} stored entries; 0 empty libraries. Library sizes range from {qc['sample_2']['min_library_size']:.0f} to {qc['sample_2']['max_library_size']:.0f} counts (median {qc['sample_2']['median_library_size']:.1f}).
- All listed marker symbols were present exactly once in the shared gene table. Parsed entry counts matched both Matrix Market headers.

## Depleted population

**{call}** is the frozen-rule depleted call. It changes from {n1:,}/{qc['sample_1']['n_cells']:,} cells ({f1:.6f}) in sample 1 to {n2:,}/{qc['sample_2']['n_cells']:,} cells ({f2:.6f}) in sample 2, for a sample-2/sample-1 fraction ratio of {ratio:.6f}. Full counts and fractions, including zero-count types, are in `composition.csv`.

## Annotation evidence and uncertainty

The evidence for every label is limited to relative expression of the marker sets in `MARKER_PANEL.tsv` after the frozen normalization. In sample 1, {uncertainty['sample_1']['n_exact_top_score_ties']} cells had an exact top-score tie and the median top-versus-runner-up score margin was {uncertainty['sample_1']['median_top_score_margin']:.6f}; in sample 2, the corresponding values were {uncertainty['sample_2']['n_exact_top_score_ties']} and {uncertainty['sample_2']['median_top_score_margin']:.6f}. Ties were retained and resolved by panel order, as required.

This deterministic marker rule is not a full biological annotation workflow. It does not model batch effects, ambient RNA, doublets, donor variability, uncertainty in marker specificity, or sampling uncertainty in the composition ratio. The depleted call is therefore a reproducible description under the supplied rule, not proof of biological loss or a causal effect. Confirmatory work would ordinarily inspect broader expression programs, technical covariates, replicate structure, and independent retinal annotations.
"""
    path.write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample1", type=Path, required=True)
    parser.add_argument("--sample2", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--markers", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    panel = read_marker_panel(args.markers)
    genes = read_gene_symbols(args.genes)
    required_markers = list(dict.fromkeys(marker for _, markers in panel for marker in markers))

    lib1, marker1, qc1 = load_sufficient_data(args.sample1, genes, required_markers)
    lib2, marker2, qc2 = load_sufficient_data(args.sample2, genes, required_markers)
    labels1, _, uncertain1 = annotate(lib1, marker1, panel, required_markers)
    labels2, _, uncertain2 = annotate(lib2, marker2, panel, required_markers)

    rows1 = composition_rows("sample_1", labels1, panel)
    rows2 = composition_rows("sample_2", labels2, panel)
    all_rows = rows1 + rows2
    with (args.output_dir / "composition.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "cell_type", "n_cells", "fraction"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_rows)

    by_type_1 = {str(row["cell_type"]): row for row in rows1}
    by_type_2 = {str(row["cell_type"]): row for row in rows2}
    eligible = [cell_type for cell_type, _ in panel if float(by_type_1[cell_type]["fraction"]) >= MIN_SAMPLE1_FRACTION]
    if not eligible:
        raise ValueError("No cell type meets the sample-1 1% threshold")
    called = min(
        eligible,
        key=lambda cell_type: float(by_type_2[cell_type]["fraction"]) / float(by_type_1[cell_type]["fraction"]),
    )
    depleted = {
        "cell_type": called,
        "sample_1_n_cells": int(by_type_1[called]["n_cells"]),
        "sample_1_fraction": float(by_type_1[called]["fraction"]),
        "sample_2_n_cells": int(by_type_2[called]["n_cells"]),
        "sample_2_fraction": float(by_type_2[called]["fraction"]),
        "sample_2_to_sample_1_fraction_ratio": float(by_type_2[called]["fraction"]) / float(by_type_1[called]["fraction"]),
        "sample_1_eligibility_threshold": MIN_SAMPLE1_FRACTION,
        "eligible_cell_types": eligible,
    }
    (args.output_dir / "depleted_call.json").write_text(
        json.dumps(depleted, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    write_report(
        args.output_dir / "report.md", panel, all_rows, depleted,
        {"sample_1": qc1, "sample_2": qc2},
        {"sample_1": uncertain1, "sample_2": uncertain2},
    )


if __name__ == "__main__":
    main()
