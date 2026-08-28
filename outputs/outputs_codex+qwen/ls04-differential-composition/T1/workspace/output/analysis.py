#!/usr/bin/env python3
"""Differential composition analysis of two retinal single-cell count matrices.

Implements the frozen annotation rule from inputs/ANNOTATION_RULE.md:
  * Each genes-x-cells MTX matrix is transposed to cells x genes.
  * Each cell is divided by its library size (total UMI counts), multiplied by
    10,000 (CP10K) and log1p-transformed.
  * Each cell is scored for every cell type in inputs/MARKER_PANEL.tsv as the
    arithmetic mean of the normalized expression of the type's marker genes.
  * Each cell is assigned to the type with the largest mean; ties follow the
    row order of the marker panel (np.argmax returns the first maximum).
  * No cells or genes are removed. Fractions use the matrix column count
    (total cells) as denominator.
  * The depleted call is the listed type with the smallest sample-2/sample-1
    fraction ratio among types with at least 1% fraction in sample 1.

Outputs (written to output/):
  composition.csv      sample,cell_type,n_cells,fraction for both samples
  depleted_call.json   depleted population + full ratio evidence
  analysis_summary.json  QC + scoring diagnostics used by report.md
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

BASE = Path(__file__).resolve().parent.parent  # workspace root
INPUT = BASE / "inputs"
OUT = BASE / "output"

SAMPLES = {"sample1": "1", "sample2": "2"}  # label -> MTX suffix
FRACTION_MIN_S1 = 0.01  # depletion candidates need >=1% in sample 1


def load_counts(tag: str) -> sparse.csr_matrix:
    """Load one MTX file and return genes x cells CSR float32 counts."""
    path = INPUT / f"differential.composition.q1.{tag}.mtx.gz"
    mat = mmread(str(path)).tocsr()
    data = mat.data
    checks = {
        "shape_genes_x_cells": list(mat.shape),
        "nnz": int(mat.nnz),
        "non_negative": bool(data.min() >= 0),
        "integer_counts": bool(np.all(np.modf(data)[0] == 0)),
    }
    mat = mat.astype(np.float32)
    mat.check_format()
    return mat, checks


def qc_stats(Xcg: sparse.csr_matrix, gene_symbols: np.ndarray) -> dict:
    """Documentation-only QC (frozen rule forbids removing cells/genes)."""
    lib = np.asarray(Xcg.sum(axis=1), dtype=np.float64).ravel()
    n_detected = np.diff(Xcg.indptr).astype(np.int64)
    mt_cols = np.flatnonzero([str(g).startswith("MT-") for g in gene_symbols])
    mt_lib = np.asarray(Xcg[:, mt_cols].sum(axis=1), dtype=np.float64).ravel()
    mt_pct = np.where(lib > 0, 100.0 * mt_lib / lib, 0.0)

    def pct(x, q):
        return float(np.percentile(x, q))

    return {
        "n_cells": int(Xcg.shape[0]),
        "zero_library_cells": int((lib == 0).sum()),
        "library_size": {
            "min": float(lib.min()), "p5": pct(lib, 5), "median": pct(lib, 50),
            "p95": pct(lib, 95), "max": float(lib.max()),
            "mean": float(lib.mean()),
        },
        "genes_detected": {
            "min": int(n_detected.min()), "median": float(np.median(n_detected)),
            "max": int(n_detected.max()),
        },
        "mito_pct": {"median": pct(mt_pct, 50), "p95": pct(mt_pct, 95)},
        "cells_lib_lt_500": int((lib < 500).sum()),
        "cells_genes_lt_200": int((n_detected < 200).sum()),
        "cells_mito_gt_20pct": int((mt_pct > 20).sum()),
    }


def cp10k_log1p(Xcg: sparse.csr_matrix, lib: np.ndarray) -> sparse.csr_matrix:
    """Per-cell library-size normalization to CP10K followed by log1p."""
    scale = np.zeros_like(lib, dtype=np.float64)
    nz = lib > 0
    scale[nz] = 10000.0 / lib[nz]
    Xn = Xcg.multiply(scale[:, None]).tocsr()
    Xn.data = np.log1p(Xn.data)
    return Xn


def score_and_assign(Xn: sparse.csr_matrix, panel: pd.DataFrame,
                     sym2idx: dict) -> tuple[np.ndarray, np.ndarray]:
    """Return (labels, scores) under the frozen rule."""
    n_cells = Xn.shape[0]
    n_types = len(panel)
    scores = np.empty((n_cells, n_types), dtype=np.float64)
    for t, (_, row) in enumerate(panel.iterrows()):
        idx = [sym2idx[g] for g in row["markers"].split(",")]
        scores[:, t] = np.asarray(Xn[:, idx].mean(axis=1), dtype=np.float64).ravel()
    labels = scores.argmax(axis=1)  # first max wins => panel row-order ties
    return labels, scores


def main() -> None:
    OUT.mkdir(exist_ok=True)
    genes = pd.read_csv(INPUT / "differential.composition.q1.genes.txt.gz", sep="\t")
    gene_symbols = genes["gene_symbols"].to_numpy()
    assert genes["gene_symbols"].is_unique
    sym2idx = {s: i for i, s in enumerate(gene_symbols)}

    panel = pd.read_csv(INPUT / "MARKER_PANEL.tsv", sep="\t")
    panel_markers = panel["markers"].str.split(",")
    missing = sorted({g for ms in panel_markers for g in ms} - set(sym2idx))
    assert not missing, f"marker genes missing from gene list: {missing}"

    types = panel["cell_type"].to_numpy()
    summary = {"rule": "frozen marker-mean rule (see inputs/ANNOTATION_RULE.md)",
               "n_cell_types": len(types),
               "n_marker_entries": int(sum(len(m) for m in panel_markers)),
               "missing_marker_genes": missing,
               "shared_markers": sorted({g for m in panel_markers for g in m
                                         if sum(g in n for n in panel_markers) > 1}),
               "samples": {}}

    comp_rows = []
    per_sample = {}
    for label, tag in SAMPLES.items():
        mat, checks = load_counts(tag)           # genes x cells
        Xcg = mat.T.tocsr()                      # cells x genes
        lib = np.asarray(Xcg.sum(axis=1), dtype=np.float64).ravel()
        qc = qc_stats(Xcg, gene_symbols)
        qc["mtx_integrity"] = checks
        Xn = cp10k_log1p(Xcg, lib)
        labels, scores = score_and_assign(Xn, panel, sym2idx)

        n_cells_total = Xcg.shape[0]
        counts = np.bincount(labels, minlength=len(types))
        fracs = counts / n_cells_total           # denominator = matrix column count
        for ct, n, f in zip(types, counts, fracs):
            comp_rows.append({"sample": label, "cell_type": ct,
                              "n_cells": int(n), "fraction": round(float(f), 6)})

        # annotation evidence / uncertainty diagnostics
        order = np.sort(scores, axis=1)
        margin = order[:, -1] - order[:, -2]     # best minus runner-up score
        winning = scores[np.arange(n_cells_total), labels]
        per_sample[label] = {"lib": lib, "labels": labels, "fracs": fracs}
        summary["samples"][label] = {
            "qc": qc,
            "n_assigned": int(n_cells_total),
            "counts": {ct: int(n) for ct, n in zip(types, counts)},
            "fractions": {ct: round(float(f), 6) for ct, f in zip(types, fracs)},
            "mean_winning_marker_score": {
                ct: round(float(winning[labels == t].mean()), 4)
                if counts[t] else None for t, ct in enumerate(types)},
            "assignment_margin": {
                "median": round(float(np.median(margin)), 4),
                "p10": round(float(np.percentile(margin, 10)), 4),
                "frac_below_0.1": round(float((margin < 0.1).mean()), 4),
            },
        }

    comp = pd.DataFrame(comp_rows)
    comp.to_csv(OUT / "composition.csv", index=False)

    # depletion call: smallest s2/s1 ratio among types >=1% in sample 1
    f1 = per_sample["sample1"]["fracs"]
    f2 = per_sample["sample2"]["fracs"]
    ratios = {ct: (float(f2[t] / f1[t]) if f1[t] > 0 else None)
              for t, ct in enumerate(types)}
    cand = [t for t in range(len(types)) if f1[t] >= FRACTION_MIN_S1]
    dep = min(cand, key=lambda t: ratios[types[t]])
    call = {
        "depleted_cell_type": types[dep],
        "sample1": "sample1", "sample2": "sample2",
        "sample1_fraction": round(float(f1[dep]), 6),
        "sample2_fraction": round(float(f2[dep]), 6),
        "sample1_n_cells": int(per_sample["sample1"]["labels"].tolist().count(dep)),
        "sample2_n_cells": int(per_sample["sample2"]["labels"].tolist().count(dep)),
        "ratio_sample2_over_sample1": round(ratios[types[dep]], 6),
        "criterion": ("smallest sample-2/sample-1 fraction ratio among cell types "
                      f"with fraction >= {FRACTION_MIN_S1:.0%} in sample 1 "
                      "(frozen rule in inputs/ANNOTATION_RULE.md)"),
        "candidate_types": [types[t] for t in cand],
        "all_fraction_ratios": {ct: (round(r, 6) if r is not None else None)
                                for ct, r in ratios.items()},
    }
    (OUT / "depleted_call.json").write_text(json.dumps(call, indent=2) + "\n")
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(json.dumps({"composition_rows": len(comp),
                      "depleted": call["depleted_cell_type"],
                      "ratio": call["ratio_sample2_over_sample1"],
                      "s1_frac": call["sample1_fraction"],
                      "s2_frac": call["sample2_fraction"]}, indent=2))


if __name__ == "__main__":
    main()
