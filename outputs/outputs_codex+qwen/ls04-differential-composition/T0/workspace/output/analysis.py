#!/usr/bin/env python3
"""
Differential composition analysis of two retinal single-cell expression matrices.

Task
----
Compare two retinal single-cell count matrices (sample 1 and sample 2) and identify
the cell population that is severely depleted in sample 2.

Inputs (under ./inputs)
-----------------------
- differential.composition.q1.1.mtx.gz   : sample 1 counts, genes x cells (MTX, integer)
- differential.composition.q1.2.mtx.gz   : sample 2 counts, genes x cells (MTX, integer)
- differential.composition.q1.genes.txt.gz : gene_ids / gene_symbols (one per matrix row)
- MARKER_PANEL.tsv                       : cell_type -> comma-separated marker gene symbols
- ANNOTATION_RULE.md                     : the frozen, reproducible annotation rule

Frozen annotation rule (from ANNOTATION_RULE.md)
------------------------------------------------
For each matrix: transpose to cells x genes; divide every cell by its library size;
multiply by 10,000; apply log1p. For each cell type in MARKER_PANEL.tsv, compute the
arithmetic mean across its listed marker genes. Assign each cell to the type with the
largest mean; ties follow the row order in the marker panel. Do NOT remove cells or
genes. Fractions use the matrix column count (number of cells) as denominator. The
depleted call is the listed type with the smallest sample-2/sample-1 fraction ratio
among types having at least 1% fraction in sample 1.

Outputs (under ./output)
------------------------
- composition.csv            : sample, cell_type, n_cells, fraction
- depleted_call.json         : the depleted population + supporting evidence
- qc.json                    : per-sample QC metrics (documentation only; no filtering)
- annotation_evidence.json   : marker-score evidence + uncertainty per cell type/sample
- analysis.py                : this script

Notes
-----
- This marker-score rule is intentionally simple and is frozen for reproducibility of
  the benchmark artifact; it is NOT a substitute for a full biological annotation
  workflow. QC below is descriptive only -- no cells or genes are removed.
- Deterministic: numpy argmax returns the first maximum, matching the panel-order tie
  rule because cell types are kept in marker-panel row order throughout.
"""

from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
INPUT_DIR = "inputs"
OUTPUT_DIR = "output"

MATRIX_FILES = {
    "sample_1": os.path.join(INPUT_DIR, "differential.composition.q1.1.mtx.gz"),
    "sample_2": os.path.join(INPUT_DIR, "differential.composition.q1.2.mtx.gz"),
}
GENES_FILE = os.path.join(INPUT_DIR, "differential.composition.q1.genes.txt.gz")
PANEL_FILE = os.path.join(INPUT_DIR, "MARKER_PANEL.tsv")

MIN_FRACTION_S1 = 0.01  # depleted call restricted to types >= 1% in sample 1
RNG_SEED = 0            # not strictly needed (deterministic), kept for transparency

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ----------------------------------------------------------------------------
# Loading helpers
# ----------------------------------------------------------------------------
def load_matrix(path: str) -> sparse.csc_matrix:
    """Load an MTX count matrix as genes x cells CSC sparse."""
    mat = mmread(path, spmatrix=True)
    return mat.tocsc()


def load_genes(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def load_panel(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def build_symbol_index(genes_df: pd.DataFrame) -> dict:
    """Map gene symbol -> matrix row index (no duplicate symbols verified upstream)."""
    return {sym: i for i, sym in enumerate(genes_df["gene_symbols"].to_numpy())}


# ----------------------------------------------------------------------------
# Normalization (frozen rule)
# ----------------------------------------------------------------------------
def normalize_genes_x_cells(mat: sparse.csc_matrix):
    """CP10k + log1p normalization, keeping genes x cells orientation.

    Divide every cell (column) by its library size, multiply by 10,000, log1p.
    Returns (normalized sparse matrix, per-cell library sizes).
    """
    lib = np.asarray(mat.sum(axis=0)).ravel().astype(np.float64)  # per-cell totals
    safe = np.where(lib > 0, lib, 1.0)  # guard against empty droplets (none expected)
    inv = sparse.diags(1.0 / safe)
    normed = (mat @ inv) * 10000.0
    normed = normed.log1p()  # log1p(0)=0 so sparsity structure is preserved
    return normed, lib


# ----------------------------------------------------------------------------
# Marker scoring
# ----------------------------------------------------------------------------
def build_marker_matrix(panel_df: pd.DataFrame, sym_to_idx: dict, n_genes: int):
    """Build a (genes x types) weight matrix so that X_norm @ M yields, for each
    cell and each type, the arithmetic mean of that type's marker genes.

    M[g, t] = 1 / (# markers of type t)  if gene g is a listed marker of type t.
    Missing markers would be skipped, but all markers are verified present upstream.
    """
    cell_types = panel_df["cell_type"].tolist()
    rows, cols, vals = [], [], []
    missing = []
    for t, markers_str in enumerate(panel_df["markers"].tolist()):
        markers = [m.strip() for m in str(markers_str).split(",") if m.strip()]
        k = len(markers)
        for m in markers:
            if m in sym_to_idx:
                rows.append(sym_to_idx[m])
                cols.append(t)
                vals.append(1.0 / k)
            else:
                missing.append(m)
    M = sparse.csr_matrix((vals, (rows, cols)), shape=(n_genes, len(cell_types)))
    # Return as dense (n_genes x n_types is small) so sparse @ dense -> dense scores.
    M = M.toarray()
    return M, cell_types, missing


def qc_metrics(mat: sparse.csc_matrix, sym_to_idx: dict, genes_df: pd.DataFrame) -> dict:
    """Descriptive QC metrics for a genes x cells matrix (no filtering applied)."""
    lib = np.asarray(mat.sum(axis=0)).ravel().astype(np.float64)      # counts per cell
    n_genes_per_cell = np.asarray((mat > 0).sum(axis=0)).ravel()        # detected genes

    # Mitochondrial fraction (MT- genes)
    mt_genes = [g for g in genes_df["gene_symbols"] if str(g).startswith("MT-")]
    mt_idx = [sym_to_idx[g] for g in mt_genes if g in sym_to_idx]
    if mt_idx:
        mt_counts = np.asarray(mat[mt_idx, :].sum(axis=0)).ravel().astype(np.float64)
        mt_frac = np.divide(mt_counts, lib, out=np.zeros_like(lib), where=lib > 0)
    else:
        mt_frac = np.zeros_like(lib)

    def _stats(x):
        x = np.asarray(x, dtype=np.float64)
        return {
            "min": float(np.min(x)),
            "median": float(np.median(x)),
            "mean": float(np.mean(x)),
            "max": float(np.max(x)),
        }

    nnz = mat.nnz
    total_possible = mat.shape[0] * mat.shape[1]
    return {
        "n_genes": int(mat.shape[0]),
        "n_cells": int(mat.shape[1]),
        "nnz": int(nnz),
        "sparsity": float(1.0 - nnz / total_possible),
        "total_counts": float(lib.sum()),
        "library_size_per_cell": _stats(lib),
        "genes_detected_per_cell": _stats(n_genes_per_cell),
        "mito_fraction_per_cell": _stats(mt_frac),
        "n_mt_genes_used": len(mt_idx),
        "cells_with_zero_counts": int(np.sum(lib == 0)),
    }


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main() -> None:
    genes_df = load_genes(GENES_FILE)
    panel_df = load_panel(PANEL_FILE)
    sym_to_idx = build_symbol_index(genes_df)
    n_genes = genes_df.shape[0]

    M, cell_types, missing_markers = build_marker_matrix(panel_df, sym_to_idx, n_genes)
    if missing_markers:
        raise RuntimeError(f"Marker genes missing from gene table: {missing_markers}")

    composition_rows = []
    qc_all = {}
    evidence_all = {}
    fractions = {}  # sample -> dict(cell_type -> fraction)

    for sample, path in MATRIX_FILES.items():
        mat = load_matrix(path)  # genes x cells
        assert mat.shape[0] == n_genes, f"{sample}: gene dimension mismatch"

        # QC (descriptive only; frozen rule forbids removing cells/genes)
        qc_all[sample] = qc_metrics(mat, sym_to_idx, genes_df)

        # Normalization + marker scoring + assignment
        normed, lib = normalize_genes_x_cells(mat)
        scores = normed.transpose() @ M           # cells x types (mean marker expr)
        scores = np.asarray(scores, dtype=np.float64)
        labels = np.argmax(scores, axis=1)        # ties -> first (panel order)

        n_cells_total = mat.shape[1]
        fractions[sample] = {}
        counts = np.bincount(labels, minlength=len(cell_types))

        # Annotation evidence + uncertainty for this sample
        # margin = score(winner) - score(runner-up); small margin => ambiguous call
        sorted_scores = np.sort(scores, axis=1)
        margin = sorted_scores[:, -1] - sorted_scores[:, -2]
        ambiguous_frac = float(np.mean(margin < 1e-6))

        per_type_evidence = {}
        for t, ct in enumerate(cell_types):
            sel = labels == t
            n_assigned = int(counts[t])
            fraction = n_assigned / n_cells_total
            fractions[sample][ct] = fraction
            composition_rows.append({
                "sample": sample,
                "cell_type": ct,
                "n_cells": n_assigned,
                "fraction": fraction,
            })
            if n_assigned > 0:
                type_margins = margin[sel]
                type_scores = scores[sel, t]
                per_type_evidence[ct] = {
                    "n_cells": n_assigned,
                    "fraction": fraction,
                    "mean_marker_score_assigned": float(np.mean(type_scores)),
                    "median_marker_score_assigned": float(np.median(type_scores)),
                    "mean_margin_over_runner_up": float(np.mean(type_margins)),
                    "frac_assigned_with_zero_margin": float(np.mean(type_margins < 1e-6)),
                }
            else:
                per_type_evidence[ct] = {
                    "n_cells": 0, "fraction": 0.0,
                    "mean_marker_score_assigned": None,
                    "median_marker_score_assigned": None,
                    "mean_margin_over_runner_up": None,
                    "frac_assigned_with_zero_margin": None,
                }
        evidence_all[sample] = {
            "n_cells_total": int(n_cells_total),
            "global_ambiguous_fraction_zero_margin": ambiguous_frac,
            "per_type": per_type_evidence,
        }

    # ------------------------------------------------------------------
    # Composition CSV (all listed types for both samples)
    # ------------------------------------------------------------------
    comp_df = pd.DataFrame(composition_rows, columns=["sample", "cell_type", "n_cells", "fraction"])
    comp_path = os.path.join(OUTPUT_DIR, "composition.csv")
    comp_df.to_csv(comp_path, index=False, float_format="%.10g")

    # ------------------------------------------------------------------
    # Depleted call: smallest s2/s1 fraction ratio among types >=1% in sample 1
    # ------------------------------------------------------------------
    ratios = {}
    eligible = []
    for ct in cell_types:
        f1 = fractions["sample_1"][ct]
        f2 = fractions["sample_2"][ct]
        if f1 >= MIN_FRACTION_S1:
            eligible.append(ct)
            ratios[ct] = (f2 / f1) if f1 > 0 else float("inf")

    depleted = min(eligible, key=lambda c: ratios[c]) if eligible else None

    depleted_call = {
        "depleted_cell_type": depleted,
        "rule": (
            "Listed type with the smallest sample_2/sample_1 fraction ratio "
            "among types with >=1% fraction in sample_1."
        ),
        "min_fraction_sample_1_threshold": MIN_FRACTION_S1,
        "fraction_sample_1": fractions["sample_1"].get(depleted) if depleted else None,
        "fraction_sample_2": fractions["sample_2"].get(depleted) if depleted else None,
        "ratio_sample_2_over_sample_1": ratios.get(depleted) if depleted else None,
        "n_cells_sample_1": int(comp_df[(comp_df["sample"] == "sample_1") & (comp_df["cell_type"] == depleted)].n_cells.iloc[0]) if depleted else None,
        "n_cells_sample_2": int(comp_df[(comp_df["sample"] == "sample_2") & (comp_df["cell_type"] == depleted)].n_cells.iloc[0]) if depleted else None,
        "eligible_cell_types": eligible,
        "all_ratios_sample_2_over_sample_1": ratios,
        "annotation_method": "frozen marker-mean rule (see ANNOTATION_RULE.md)",
        "note": (
            "Simple marker-score rule frozen for reproducibility; not a full "
            "biological annotation workflow."
        ),
    }
    with open(os.path.join(OUTPUT_DIR, "depleted_call.json"), "w") as fh:
        json.dump(depleted_call, fh, indent=2)

    # ------------------------------------------------------------------
    # QC + evidence JSON
    # ------------------------------------------------------------------
    with open(os.path.join(OUTPUT_DIR, "qc.json"), "w") as fh:
        json.dump(qc_all, fh, indent=2)
    with open(os.path.join(OUTPUT_DIR, "annotation_evidence.json"), "w") as fh:
        json.dump(evidence_all, fh, indent=2)

    # ------------------------------------------------------------------
    # Console summary
    # ------------------------------------------------------------------
    print("Composition (fraction):")
    pivot = comp_df.pivot(index="cell_type", columns="sample", values="fraction")
    pivot = pivot.reindex(cell_types)
    pivot["ratio_s2_s1"] = pivot["sample_2"] / pivot["sample_1"]
    print(pivot.to_string(float_format=lambda x: f"{x:0.4f}"))
    print("\nEligible (>=1% in sample1):", eligible)
    print("DEPLETED in sample 2:", depleted)
    if depleted:
        print(f"  sample_1 fraction={fractions['sample_1'][depleted]:.4f}, "
              f"sample_2 fraction={fractions['sample_2'][depleted]:.4f}, "
              f"ratio={ratios[depleted]:.4f}")


if __name__ == "__main__":
    np.random.seed(RNG_SEED)
    main()


