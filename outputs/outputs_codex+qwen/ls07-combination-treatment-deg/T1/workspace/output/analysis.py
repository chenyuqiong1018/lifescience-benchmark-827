"""Differential expression analysis: Cisplatin_IC50_CBD_IC50 vs DMSO.

Uses PyDESeq2 0.5.0 on the six frozen samples (3 replicates per group) from
inputs/. Genes are pre-filtered to those with raw count > 10 in at least one
of the six selected samples. Group is the only design term; the combination
treatment is the numerator and DMSO the denominator.

Outputs (written to output/):
  - differential_expression.csv : full DE table keyed by Ensembl ID
  - summary.json                : run metadata and counts
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# --- anndata compatibility shim ---------------------------------------------
# PyDESeq2 0.5.0 stores 1-D arrays in AnnData.varm/obsm and later indexes with
# them. anndata >= 0.11 reshapes 1-D values to (n, 1) on assignment, which
# breaks PyDESeq2 0.5.0 (IndexError in fit_genewise_dispersions). Restore the
# pre-0.11 behavior (keep 1-D entries 1-D) for this process.
import anndata._core.aligned_mapping as _aligned_mapping

_orig_validate_value = _aligned_mapping.AlignedMappingBase._validate_value


def _validate_value_keep_1d(self, val, key):
    was_1d_ndarray = isinstance(val, np.ndarray) and val.ndim == 1
    out = _orig_validate_value(self, val, key)
    if (
        was_1d_ndarray
        and isinstance(out, np.ndarray)
        and out.shape == (val.shape[0], 1)
    ):
        out = out.reshape(-1)
    return out


_aligned_mapping.AlignedMappingBase._validate_value = _validate_value_keep_1d
# ---------------------------------------------------------------------------

import pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT = Path(__file__).resolve().parent.parent  # workspace root
INPUTS = ROOT / "inputs"
OUTDIR = ROOT / "output"

NUMERATOR = "Cisplatin_IC50_CBD_IC50"
DENOMINATOR = "DMSO"
PADJ_CUTOFF = 0.05
LFC_CUTOFF = 0.5
BASEMEAN_CUTOFF = 10.0


def main() -> None:
    # --- Load inputs -------------------------------------------------------
    layout = pd.read_csv(INPUTS / "sample_layout.csv")
    counts = pd.read_csv(INPUTS / "counts_raw_unfiltered.csv", index_col=0)
    gene_map = pd.read_csv(INPUTS / "ensg_to_gene_name.tsv", sep="\t")

    sel = layout[layout["Group"].isin([NUMERATOR, DENOMINATOR])].copy()
    assert sorted(sel["Group"].value_counts()) == [3, 3], sel
    sel["counts_col"] = sel["SampleID"].str.replace("-", "_")
    missing = [c for c in sel["counts_col"] if c not in counts.columns]
    assert not missing, missing
    samples = list(sel["counts_col"])

    # Six-sample count matrix (samples x genes) and Group metadata.
    six = counts.loc[:, samples].T
    six.index.name = "sample"
    metadata = (
        sel.set_index("counts_col")
        .loc[samples, ["SampleID", "Group"]]
        .rename_axis("sample")
    )
    assert list(metadata["Group"].value_counts().sort_index().values) == [3, 3]

    # --- Pre-fit filter: raw count > 10 in at least one of the six samples -
    keep = six.max(axis=0) > 10
    six = six.loc[:, keep].astype(int)
    n_genes_total = counts.shape[0]
    n_genes_tested = int(six.shape[1])

    # --- DESeq2 (Group as only design term; refit Cook's outliers) ---------
    dds = DeseqDataSet(
        counts=six,
        metadata=metadata,
        design="~Group",
        refit_cooks=True,
    )
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=["Group", NUMERATOR, DENOMINATOR])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        stat_res.summary()
    res = stat_res.results_df.copy()
    res.index.name = "ensembl_id"

    # --- Passing genes -----------------------------------------------------
    passed = (
        (res["padj"] < PADJ_CUTOFF)
        & (res["log2FoldChange"].abs() > LFC_CUTOFF)
        & (res["baseMean"] > BASEMEAN_CUTOFF)
    )
    passed = passed.fillna(False).astype(bool)
    res["passed"] = passed

    res = res.join(
        gene_map.set_index("ENSG")["gene_name"], on="ensembl_id", how="left"
    )
    cols = ["gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat",
            "pvalue", "padj", "passed"]
    res = res[cols]

    csv_path = OUTDIR / "differential_expression.csv"
    res.to_csv(csv_path)  # NaN -> empty cell (null preserved, never zero)

    # --- Summary -------------------------------------------------------------
    sig = res.loc[res["passed"]]
    summary = {
        "contrast": f"{NUMERATOR} vs {DENOMINATOR}",
        "numerator": NUMERATOR,
        "denominator": DENOMINATOR,
        "design": "~Group",
        "samples": {
            NUMERATOR: sel.loc[sel["Group"] == NUMERATOR, "SampleID"].tolist(),
            DENOMINATOR: sel.loc[sel["Group"] == DENOMINATOR, "SampleID"].tolist(),
        },
        "prefilter": "raw count > 10 in at least one of the six selected samples",
        "n_genes_in_raw_matrix": n_genes_total,
        "n_genes_tested": n_genes_tested,
        "n_genes_dropped_by_prefilter": n_genes_total - n_genes_tested,
        "n_padj_null": int(res["padj"].isna().sum()),
        "thresholds": {
            "padj": PADJ_CUTOFF,
            "abs_log2FoldChange": LFC_CUTOFF,
            "baseMean": BASEMEAN_CUTOFF,
        },
        "n_passed": int(passed.sum()),
        "n_up": int((passed & (res["log2FoldChange"] > 0)).sum()),
        "n_down": int((passed & (res["log2FoldChange"] < 0)).sum()),
        "software": {
            "pydeseq2": pydeseq2.__version__,
            "pandas": pd.__version__,
            "refit_cooks": True,
        },
    }
    with open(OUTDIR / "summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))
    print("top passing genes by padj:")
    print(sig.sort_values("padj").head(10).to_string())


if __name__ == "__main__":
    main()

