#!/usr/bin/env python3
"""Differential composition analysis of two retinal single-cell matrices.

Implements the frozen composition rule from inputs/ANNOTATION_RULE.md:
  * Transpose each MTX matrix to cells-by-genes.
  * Divide each cell by its library size, multiply by 10,000, log1p.
  * For each cell type in inputs/MARKER_PANEL.tsv, compute the arithmetic
    mean of its listed marker genes per cell.
  * Assign each cell to the type with the largest mean; ties follow the
    row order of the marker panel.
  * No cells or genes are removed. Fractions use the matrix column count
    (total cells) as denominator.
  * Depleted call: among types with >=1% fraction in sample 1, the type
    with the smallest sample2/sample1 fraction ratio.

Outputs: output/composition.csv, output/depleted_call.json,
         output/qc.json (supporting metrics used by report.md).
"""
from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.io import mmread

ROOT = Path(__file__).resolve().parent.parent
INPUTS = ROOT / "inputs"
OUT = ROOT / "output"
SAMPLE_FILES = {
    "sample1": INPUTS / "differential.composition.q1.1.mtx.gz",
    "sample2": INPUTS / "differential.composition.q1.2.mtx.gz",
}
GENES_FILE = INPUTS / "differential.composition.q1.genes.txt.gz"
PANEL_FILE = INPUTS / "MARKER_PANEL.tsv"
SEED = 0


def load_panel():
    types, markers = [], []
    with open(PANEL_FILE, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            types.append(row["cell_type"].strip())
            markers.append([m.strip() for m in row["markers"].split(",")])
    return types, markers


def load_gene_symbols():
    with gzip.open(GENES_FILE, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        syms = [line.rstrip("\n").split("\t")[1] for line in fh]
    return syms


def annotate(X_gc: sp.csr_matrix, sym_to_idx, types, markers):
    """Return (labels, scores, margins) for a cells-by-genes CSR matrix."""
    n_cells = X_gc.shape[0]
    scores = np.full((n_cells, len(types)), -np.inf, dtype=np.float64)
    for j, (ct, ms) in enumerate(zip(types, markers)):
        idx = [sym_to_idx[m] for m in ms]  # all markers verified present
        sub = X_gc[:, idx]
        scores[:, j] = np.asarray(sub.mean(axis=1)).ravel()
    labels = scores.argmax(axis=1)  # argmax -> first max = panel row order on ties
    # margin: gap between best and second-best score (ambiguity diagnostic)
    top2 = np.sort(scores, axis=1)[:, -2:]
    margins = top2[:, 1] - top2[:, 0]
    return labels, scores, margins


def main():
    OUT.mkdir(exist_ok=True)
    types, markers = load_panel()
    syms = load_gene_symbols()
    sym_to_idx = {s: i for i, s in enumerate(syms)}
    missing = [m for ms in markers for m in ms if m not in sym_to_idx]
    assert not missing, f"missing marker genes: {missing}"

    composition_rows = []
    qc = {"marker_genes_missing": [], "samples": {}}
    fractions = {}  # sample -> {cell_type: fraction}
    counts = {}

    for sname, path in SAMPLE_FILES.items():
        X = mmread(str(path)).tocsr()  # genes x cells
        n_genes, n_cells = X.shape
        lib = np.asarray(X.sum(axis=0)).ravel().astype(np.float64)
        genes_detected = np.asarray((X > 0).sum(axis=0)).ravel()
        # cells-by-genes, CP10K + log1p (frozen rule)
        Xcg = X.transpose().tocsr()
        Xcg = sp.diags(10000.0 / lib) @ Xcg  # CP10K, keeps CSR
        np.log1p(Xcg.data, out=Xcg.data)
        Xcg.eliminate_zeros()

        labels, scores, margins = annotate(Xcg, sym_to_idx, types, markers)
        cnt = np.bincount(labels, minlength=len(types))

        fractions[sname] = {ct: cnt[j] / n_cells for j, ct in enumerate(types)}
        counts[sname] = {ct: int(cnt[j]) for j, ct in enumerate(types)}
        for j, ct in enumerate(types):
            composition_rows.append(
                {
                    "sample": sname,
                    "cell_type": ct,
                    "n_cells": int(cnt[j]),
                    "fraction": round(cnt[j] / n_cells, 6),
                }
            )

        # per-type mean marker score among cells assigned to the type
        assigned_mean = {}
        for j, ct in enumerate(types):
            sel = labels == j
            assigned_mean[ct] = (
                float(scores[sel, j].mean()) if sel.any() else None
            )
        qc["samples"][sname] = {
            "file": path.name,
            "n_cells": int(n_cells),
            "n_genes_in_matrix": int(n_genes),
            "library_size_mean": float(lib.mean()),
            "library_size_median": float(np.median(lib)),
            "library_size_min": float(lib.min()),
            "library_size_max": float(lib.max()),
            "genes_per_cell_median": float(np.median(genes_detected)),
            "genes_per_cell_mean": float(genes_detected.mean()),
            "median_assigned_marker_score_by_type": assigned_mean,
            "cells_with_top2_margin_lt_0.1": int((margins < 0.1).sum()),
            "cells_with_top2_margin_lt_0.01": int((margins < 0.01).sum()),
            "max_score_per_cell_median": float(scores.max(axis=1).mean()),
        }

    # depleted call per frozen rule
    eligible = [ct for ct in types if fractions["sample1"][ct] >= 0.01]
    ratios = {
        ct: fractions["sample2"][ct] / fractions["sample1"][ct] for ct in eligible
    }
    depleted = min(ratios, key=lambda ct: ratios[ct])

    ranking = sorted(
        (
            {
                "cell_type": ct,
                "fraction_sample1": fractions["sample1"][ct],
                "fraction_sample2": fractions["sample2"][ct],
                "ratio_sample2_over_sample1": ratios[ct],
                "eligible_sample1_fraction_ge_1pct": True,
            }
            for ct in eligible
        ),
        key=lambda d: d["ratio_sample2_over_sample1"],
    )
    ineligible = [
        {
            "cell_type": ct,
            "fraction_sample1": fractions["sample1"][ct],
            "fraction_sample2": fractions["sample2"][ct],
            "eligible_sample1_fraction_ge_1pct": False,
        }
        for ct in types
        if ct not in eligible
    ]

    call = {
        "depleted_cell_type": depleted,
        "criterion": (
            "Among cell types with >=1% fraction in sample 1, the one with the "
            "smallest sample2/sample1 fraction ratio (frozen rule in "
            "inputs/ANNOTATION_RULE.md)."
        ),
        "fraction_sample1": fractions["sample1"][depleted],
        "fraction_sample2": fractions["sample2"][depleted],
        "ratio_sample2_over_sample1": ratios[depleted],
        "ranking_by_ratio_among_eligible": ranking,
        "ineligible_types": ineligible,
    }

    with open(OUT / "composition.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "cell_type", "n_cells", "fraction"])
        w.writeheader()
        w.writerows(composition_rows)
    with open(OUT / "depleted_call.json", "w") as fh:
        json.dump(call, fh, indent=2)
    with open(OUT / "qc.json", "w") as fh:
        json.dump(qc, fh, indent=2)

    print(json.dumps(call, indent=2))
    for s in ("sample1", "sample2"):
        print(s, qc["samples"][s]["n_cells"], "cells; library median",
              int(qc["samples"][s]["library_size_median"]))


if __name__ == "__main__":
    main()
