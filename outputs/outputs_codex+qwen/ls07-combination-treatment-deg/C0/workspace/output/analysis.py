#!/usr/bin/env python3
"""Differential expression: Cisplatin_IC50_CBD_IC50 (numerator) vs DMSO (denominator).

Uses PyDESeq2 0.5.0 with refit_cooks=True and a single design term (Group,
encoded as `condition`, formula "~condition"). Only the 3 replicates of each
of the two groups are used. Genes are retained pre-fit when at least one of
the six selected samples has a raw count > 10. Unavailable adjusted p-values
are preserved as null (empty CSV cells / JSON null), never converted to zero.

Outputs (written next to this script, i.e. output/):
  differential_expression.csv, summary.json
"""
import json
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------
# Compatibility shim: anndata >= 0.13 reshapes 1-D aligned arrays (varm/obsm)
# to shape (n, 1) on assignment, which breaks PyDESeq2 0.5.0 (it relies on
# the anndata <= 0.12 behaviour of round-tripping 1-D arrays). Restore the
# old behaviour for 1-D ndarrays only; all other values go through the
# original validation unchanged.
# --------------------------------------------------------------------------
from anndata._core import aligned_mapping as _aligned_mapping

_orig_validate_value = _aligned_mapping.AlignedMappingBase._validate_value


def _validate_value_keep_1d(self, val, key):
    if isinstance(val, np.ndarray) and val.ndim == 1:
        for i, axis in enumerate(self.axes):
            if self.parent.shape[axis] == len(val):
                continue
            raise ValueError(
                f"Value passed to {type(self).__name__} key {key!r} has length "
                f"{len(val)} but parent shape[{axis}] is {self.parent.shape[axis]}."
            )
        return val
    return _orig_validate_value(self, val, key)


_aligned_mapping.AlignedMappingBase._validate_value = _validate_value_keep_1d
# --------------------------------------------------------------------------

import pandas as pd
import pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

HERE = Path(__file__).resolve().parent          # output/
INPUTS = HERE.parent / "inputs"

NUMERATOR = "Cisplatin_IC50_CBD_IC50"  # combination treatment
DENOMINATOR = "DMSO"

# ---------------------------------------------------------------- load inputs
counts = pd.read_csv(INPUTS / "counts_raw_unfiltered.csv", index_col=0)
layout = pd.read_csv(INPUTS / "sample_layout.csv")
gene_map = pd.read_csv(INPUTS / "ensg_to_gene_name.tsv", sep="\t")

# ----------------------------------------------------------- sample selection
layout = layout[layout["Group"].isin([NUMERATOR, DENOMINATOR])].copy()
assert sorted(layout["Group"].value_counts()) == [3, 3], "expected 3 replicates per group"
layout["count_col"] = layout["SampleID"].str.replace("-", "_")
assert layout["count_col"].isin(counts.columns).all(), "sample mismatch"

sample_order = layout.sort_values("count_col")["count_col"].tolist()
counts6 = counts.loc[:, sample_order].astype(int)

# --------------------------------------------------------- pre-fit gene filter
keep = (counts6 > 10).any(axis=1)
counts_filt = counts6.loc[keep]
n_genes_raw, n_genes_kept = counts.shape[0], counts_filt.shape[0]

# -------------------------------------------------------------- DESeq2 model
clinical = layout.set_index("count_col").loc[sample_order, ["SampleID", "Group"]].copy()
clinical.columns = ["sample_id", "condition"]
clinical["condition"] = clinical["condition"].astype(str)

dds = DeseqDataSet(
    counts=counts_filt.T,                 # samples x genes, as required
    metadata=clinical,
    design="~condition",                  # Group is the only design term
    refit_cooks=True,
)
dds.deseq2()

stat_res = DeseqStats(
    dds,
    contrast=["condition", NUMERATOR, DENOMINATOR],  # combination / DMSO
)
stat_res.summary()
res = stat_res.results_df.copy()          # baseMean, log2FoldChange, lfcSE, stat, pvalue, padj

# ------------------------------------------------------------- gene name join
gene_map = gene_map.drop_duplicates("ENSG").set_index("ENSG")["gene_name"]
res.insert(0, "gene_name", res.index.map(gene_map))   # unmapped -> NaN (empty cell)

# ------------------------------------------------------------ pass/fail calls
padj = res["padj"]
lfc = res["log2FoldChange"]
base_mean = res["baseMean"]
res["significant"] = (
    padj.notna() & (padj < 0.05) & (lfc.abs() > 0.5) & (base_mean > 10)
)

# ------------------------------------------------------------------ CSV write
res.index.name = "gene_id"                # Ensembl IDs key the DE table
res.to_csv(HERE / "differential_expression.csv")  # NaN -> empty cell (null)

# --------------------------------------------------------------- summary.json
sig = res["significant"]
summary = {
    "comparison": f"{NUMERATOR} vs {DENOMINATOR}",
    "numerator": NUMERATOR,
    "denominator": DENOMINATOR,
    "design": "~ condition (Group as only design term)",
    "software": {"pydeseq2": pydeseq2.__version__, "pandas": pd.__version__},
    "refit_cooks": True,
    "samples": {
        NUMERATOR: layout.loc[layout["Group"] == NUMERATOR, "SampleID"].tolist(),
        DENOMINATOR: layout.loc[layout["Group"] == DENOMINATOR, "SampleID"].tolist(),
    },
    "prefilter_rule": "keep gene if raw count > 10 in at least one of the 6 selected samples",
    "n_genes_input": int(n_genes_raw),
    "n_genes_after_prefilter": int(n_genes_kept),
    "n_genes_tested": int(res.shape[0]),
    "n_padj_unavailable": int(res["padj"].isna().sum()),
    "thresholds": {"padj": 0.05, "abs_log2FoldChange": 0.5, "baseMean": 10.0},
    "n_significant": int(sig.sum()),
    "n_significant_up": int((sig & (lfc > 0)).sum()),
    "n_significant_down": int((sig & (lfc < 0)).sum()),
}
with open(HERE / "summary.json", "w", encoding="utf-8") as fh:
    json.dump(summary, fh, indent=2, allow_nan=False)

print(json.dumps(summary, indent=2))
print("Top genes by padj:")
print(res.dropna(subset=["padj"]).sort_values("padj").head(10).to_string())
