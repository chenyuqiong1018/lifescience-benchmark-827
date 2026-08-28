#!/usr/bin/env python3
"""Frozen retinal composition analysis with QC, uncertainty, and vector summary."""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import math
from pathlib import Path

import numpy as np


NORMALIZATION_TOTAL = 10_000.0
SAMPLE1_THRESHOLD = 0.01
Z_95 = 1.959963984540054


def marker_panel(path: Path) -> list[tuple[str, tuple[str, ...]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["cell_type", "markers"]:
            raise ValueError("Unexpected marker-panel schema")
        rows = [(r["cell_type"], tuple(x.strip() for x in r["markers"].split(",") if x.strip())) for r in reader]
    if not rows or any(not markers for _, markers in rows):
        raise ValueError("Every panel row requires markers")
    return rows


def gene_symbols(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "gene_symbols" not in reader.fieldnames:
            raise ValueError("Gene table lacks gene_symbols")
        return [row["gene_symbols"] for row in reader]


def matrix_header(handle) -> tuple[int, int, int]:
    banner = handle.readline().strip()
    if banner != "%%MatrixMarket matrix coordinate integer general":
        raise ValueError(f"Unsupported Matrix Market header: {banner}")
    line = handle.readline()
    while line.startswith("%"):
        line = handle.readline()
    dimensions = tuple(int(x) for x in line.split())
    if len(dimensions) != 3:
        raise ValueError("Invalid Matrix Market dimension line")
    return dimensions


def scan_counts(
    path: Path, genes: list[str], ordered_markers: list[str]
) -> tuple[np.ndarray, np.ndarray, dict[str, int | float]]:
    marker_position = {name: i for i, name in enumerate(ordered_markers)}
    occurrences = {name: 0 for name in ordered_markers}
    matrix_row_to_marker: dict[int, int] = {}
    for matrix_row, symbol in enumerate(genes, start=1):
        if symbol in marker_position:
            occurrences[symbol] += 1
            matrix_row_to_marker[matrix_row] = marker_position[symbol]
    invalid = {name: n for name, n in occurrences.items() if n != 1}
    if invalid:
        raise ValueError(f"Marker symbols must be present exactly once: {invalid}")

    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        n_genes, n_cells, declared_entries = matrix_header(handle)
        if n_genes != len(genes):
            raise ValueError("Gene dimension does not match gene table")
        library = np.zeros(n_cells, dtype=np.float64)
        selected = np.zeros((len(ordered_markers), n_cells), dtype=np.float64)
        parsed_entries = 0
        for line in handle:
            if not line.strip() or line.startswith("%"):
                continue
            row, column, count = (int(x) for x in line.split())
            if not (1 <= row <= n_genes and 1 <= column <= n_cells) or count < 0:
                raise ValueError("Invalid Matrix Market coordinate/count")
            cell = column - 1
            library[cell] += count
            selected_row = matrix_row_to_marker.get(row)
            if selected_row is not None:
                selected[selected_row, cell] += count
            parsed_entries += 1
    if parsed_entries != declared_entries:
        raise ValueError(f"Declared {declared_entries} entries, parsed {parsed_entries}")
    if np.any(library <= 0):
        raise ValueError("All cell libraries must be positive")
    qc = {
        "n_genes": n_genes,
        "n_cells": n_cells,
        "n_nonzero_entries": declared_entries,
        "empty_cells": int(np.sum(library == 0)),
        "min_library_size": int(library.min()),
        "median_library_size": float(np.median(library)),
        "max_library_size": int(library.max()),
    }
    return library, selected, qc


def frozen_annotation(
    library: np.ndarray,
    marker_counts: np.ndarray,
    ordered_markers: list[str],
    panel: list[tuple[str, tuple[str, ...]]],
) -> tuple[np.ndarray, dict[str, int | float]]:
    index = {name: i for i, name in enumerate(ordered_markers)}
    transformed = np.log1p(marker_counts / library[None, :] * NORMALIZATION_TOTAL)
    scores = np.vstack(
        [transformed[[index[name] for name in markers]].mean(axis=0) for _, markers in panel]
    )
    labels = np.argmax(scores, axis=0)  # documented first-row tie behavior
    sorted_scores = np.sort(scores, axis=0)
    margin = sorted_scores[-1] - sorted_scores[-2]
    diagnostics = {
        "exact_top_score_ties": int(np.count_nonzero(margin == 0)),
        "median_top_score_margin": float(np.median(margin)),
        "fraction_margin_below_0_05": float(np.mean(margin < 0.05)),
    }
    return labels, diagnostics


def rows_for_sample(
    sample: str, labels: np.ndarray, panel: list[tuple[str, tuple[str, ...]]]
) -> list[dict[str, object]]:
    counts = np.bincount(labels, minlength=len(panel))
    return [
        {
            "sample": sample,
            "cell_type": name,
            "n_cells": int(counts[i]),
            "fraction": float(counts[i] / labels.size),
        }
        for i, (name, _) in enumerate(panel)
    ]


def wilson(count: int, total: int) -> tuple[float, float]:
    p = count / total
    d = 1 + Z_95 * Z_95 / total
    center = (p + Z_95 * Z_95 / (2 * total)) / d
    half = Z_95 * math.sqrt(p * (1 - p) / total + Z_95 * Z_95 / (4 * total * total)) / d
    return center - half, center + half


def ratio_interval(x1: int, n1: int, x2: int, n2: int) -> tuple[float, float]:
    estimate = (x2 / n2) / (x1 / n1)
    se = math.sqrt(1 / x2 - 1 / n2 + 1 / x1 - 1 / n1)
    return math.exp(math.log(estimate) - Z_95 * se), math.exp(math.log(estimate) + Z_95 * se)


def composition_svg(
    path: Path,
    rows1: list[dict[str, object]],
    rows2: list[dict[str, object]],
    called: str,
) -> None:
    width, height = 880, 680
    left, right, top, bottom = 190, 35, 75, 70
    plot_w, plot_h = width - left - right, height - top - bottom
    maximum = max(float(row["fraction"]) for row in rows1 + rows2) * 1.08
    if maximum <= 0:
        maximum = 1.0

    def x(value: float) -> float:
        return left + value / maximum * plot_w

    spacing = plot_h / len(rows1)
    plot_elements: list[str] = []
    for i, (a, b) in enumerate(zip(rows1, rows2)):
        y = top + (i + 0.5) * spacing
        name = str(a["cell_type"])
        if name != str(b["cell_type"]):
            raise ValueError("Composition rows are not aligned")
        p1, p2 = float(a["fraction"]), float(b["fraction"])
        lo1, hi1 = wilson(int(a["n_cells"]), sum(int(r["n_cells"]) for r in rows1))
        lo2, hi2 = wilson(int(b["n_cells"]), sum(int(r["n_cells"]) for r in rows2))
        if name == called:
            plot_elements.append(f'<rect x="0" y="{y-spacing/2:.2f}" width="{width}" height="{spacing:.2f}" fill="#FFF3E6"/>')
        weight = "bold" if name == called else "normal"
        plot_elements.append(f'<text x="{left-12}" y="{y+4:.2f}" text-anchor="end" font-weight="{weight}">{html.escape(name)}</text>')
        plot_elements.append(f'<line x1="{x(p1):.2f}" y1="{y:.2f}" x2="{x(p2):.2f}" y2="{y:.2f}" stroke="#B0B0B0" stroke-width="1"/>')
        plot_elements.append(f'<line x1="{x(lo1):.2f}" y1="{y-3:.2f}" x2="{x(hi1):.2f}" y2="{y-3:.2f}" stroke="#0072B2" stroke-width="2"/>')
        plot_elements.append(f'<circle cx="{x(p1):.2f}" cy="{y-3:.2f}" r="4" fill="#0072B2"/>')
        plot_elements.append(f'<line x1="{x(lo2):.2f}" y1="{y+3:.2f}" x2="{x(hi2):.2f}" y2="{y+3:.2f}" stroke="#D55E00" stroke-width="2"/>')
        plot_elements.append(f'<rect x="{x(p2)-3.5:.2f}" y="{y-0.5:.2f}" width="7" height="7" fill="#D55E00"/>')

    tick_values = np.linspace(0, maximum, 6)
    ticks = "".join(
        f'<line x1="{x(float(v)):.2f}" y1="{top+plot_h}" x2="{x(float(v)):.2f}" y2="{top+plot_h+5}" stroke="#222"/>'
        f'<text x="{x(float(v)):.2f}" y="{top+plot_h+22}" text-anchor="middle">{100*float(v):.1f}%</text>'
        for v in tick_values
    )
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">
<title id="title">Retinal cell composition by sample</title>
<desc id="desc">Paired point estimates and descriptive cell-level Wilson intervals for sixteen marker-defined cell types. The frozen-rule depleted call, {html.escape(called)}, is highlighted.</desc>
<rect width="100%" height="100%" fill="white"/>
<g font-family="Arial, Helvetica, sans-serif" font-size="12" fill="#222">
<text x="{left}" y="28" font-size="18" font-weight="bold">Retinal cell composition</text>
<circle cx="{left}" cy="51" r="4" fill="#0072B2"/><text x="{left+10}" y="55">Sample 1</text>
<rect x="{left+105}" y="47" width="8" height="8" fill="#D55E00"/><text x="{left+120}" y="55">Sample 2</text>
{''.join(plot_elements)}
<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#222"/>
{ticks}
<text x="{left+plot_w/2:.2f}" y="{height-18}" text-anchor="middle" font-size="14">Fraction of all cells</text>
</g></svg>'''
    path.write_text(svg, encoding="utf-8")


def write_report(
    path: Path,
    call: dict[str, object],
    qc1: dict[str, int | float],
    qc2: dict[str, int | float],
    d1: dict[str, int | float],
    d2: dict[str, int | float],
) -> None:
    x1, n1 = int(call["sample_1_n_cells"]), int(qc1["n_cells"])
    x2, n2 = int(call["sample_2_n_cells"]), int(qc2["n_cells"])
    p1, p2 = float(call["sample_1_fraction"]), float(call["sample_2_fraction"])
    rr = float(call["sample_2_to_sample_1_fraction_ratio"])
    ci1, ci2, cir = wilson(x1, n1), wilson(x2, n2), ratio_interval(x1, n1, x2, n2)
    text = f"""# Retinal composition comparison

## QC and frozen analysis

Sample 1 contains {qc1['n_genes']:,} genes x {n1:,} cells and {qc1['n_nonzero_entries']:,} stored count entries; sample 2 contains {qc2['n_genes']:,} genes x {n2:,} cells and {qc2['n_nonzero_entries']:,} entries. No cell or gene was removed. All libraries were positive, header entry totals matched the parsed entries, and every listed marker symbol occurred exactly once.

Each cell's raw integer counts were divided by its full library size, multiplied by 10,000, and transformed with `log1p`. For each panel row, the arithmetic mean of its listed markers was computed; the largest score assigned the type, with panel order resolving ties. Fractions use the complete matrix column count. The call was restricted to types with sample-1 fraction >= 1% and selected by the minimum sample-2/sample-1 fraction ratio.

Library-size QC: sample 1 range {qc1['min_library_size']:,}–{qc1['max_library_size']:,}, median {qc1['median_library_size']:.1f}; sample 2 range {qc2['min_library_size']:,}–{qc2['max_library_size']:,}, median {qc2['median_library_size']:.1f}.

## Depleted population

The depleted call is **{call['cell_type']}**. It decreases from {x1:,}/{n1:,} cells (fraction {p1:.6f}) to {x2:,}/{n2:,} (fraction {p2:.6f}); the sample-2/sample-1 fraction ratio is {rr:.6f}. Thus, under the frozen labels, the sample-2 proportion is about {100*rr:.1f}% of the sample-1 proportion.

Conditional on treating cells as independent observations, descriptive 95% Wilson intervals are [{ci1[0]:.6f}, {ci1[1]:.6f}] and [{ci2[0]:.6f}, {ci2[1]:.6f}], and the log-ratio approximation is [{cir[0]:.6f}, {cir[1]:.6f}]. These intervals are shown in `composition_profile.svg`; they do not represent donor-level uncertainty.

## Annotation evidence, model scope, and uncertainty

Annotation evidence is restricted to the supplied marker sets. Sample 1 has {d1['exact_top_score_ties']} exact winning-score ties (median winning margin {d1['median_top_score_margin']:.6f}); sample 2 has {d2['exact_top_score_ties']} (median margin {d2['median_top_score_margin']:.6f}). Ties were retained and resolved as required.

The matrices are raw integer counts, which would be appropriate input to probabilistic single-cell models. However, scVI/scANVI training, scGPT embeddings, vocabulary-based gene dropping, learned clustering, or label transfer would change the explicitly frozen normalization and annotation rule. They were therefore not substituted for the required calculation. No checkpoint, reference labels, donor/batch covariates, or GPU modeling is needed to reproduce this artifact.

Biological uncertainty remains: the rule does not handle doublets, ambient RNA, marker overlap, continuous states, batch effects, or donor replication. With one aggregate matrix per sample, cells may be pseudoreplicates; the analysis cannot establish population-level variance, a causal loss mechanism, or generalization beyond these matrices. The result is a deterministic composition call under the supplied rule.
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

    panel = marker_panel(args.markers)
    genes = gene_symbols(args.genes)
    ordered_markers = list(dict.fromkeys(marker for _, markers in panel for marker in markers))
    lib1, selected1, qc1 = scan_counts(args.sample1, genes, ordered_markers)
    lib2, selected2, qc2 = scan_counts(args.sample2, genes, ordered_markers)
    labels1, diag1 = frozen_annotation(lib1, selected1, ordered_markers, panel)
    labels2, diag2 = frozen_annotation(lib2, selected2, ordered_markers, panel)
    rows1 = rows_for_sample("sample_1", labels1, panel)
    rows2 = rows_for_sample("sample_2", labels2, panel)

    with (args.output_dir / "composition.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, ["sample", "cell_type", "n_cells", "fraction"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows1 + rows2)

    s1 = {str(row["cell_type"]): row for row in rows1}
    s2 = {str(row["cell_type"]): row for row in rows2}
    eligible = [str(row["cell_type"]) for row in rows1 if float(row["fraction"]) >= SAMPLE1_THRESHOLD]
    called = min(eligible, key=lambda name: float(s2[name]["fraction"]) / float(s1[name]["fraction"]))
    call = {
        "cell_type": called,
        "sample_1_n_cells": int(s1[called]["n_cells"]),
        "sample_1_fraction": float(s1[called]["fraction"]),
        "sample_2_n_cells": int(s2[called]["n_cells"]),
        "sample_2_fraction": float(s2[called]["fraction"]),
        "sample_2_to_sample_1_fraction_ratio": float(s2[called]["fraction"]) / float(s1[called]["fraction"]),
        "sample_1_eligibility_threshold": SAMPLE1_THRESHOLD,
        "eligible_cell_types": eligible,
    }
    (args.output_dir / "depleted_call.json").write_text(json.dumps(call, indent=2) + "\n", encoding="utf-8")
    composition_svg(args.output_dir / "composition_profile.svg", rows1, rows2, called)
    write_report(args.output_dir / "report.md", call, qc1, qc2, diag1, diag2)


if __name__ == "__main__":
    main()
