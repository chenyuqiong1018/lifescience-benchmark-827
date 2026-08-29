#!/usr/bin/env python
"""Differential expression: Cisplatin_IC50_CBD_IC50 vs DMSO with PyDESeq2 0.5.0.

Uses only the 3 replicates of Cisplatin_IC50_CBD_IC50 (numerator) and the 3
replicates of DMSO (denominator). Design: Group as the only term.
Prefilter: keep genes with raw count > 10 in at least one of the six samples.
Pass rule: padj < 0.05, abs(log2FoldChange) > 0.5, baseMean > 10.
Unavailable adjusted p-values are preserved as null (empty CSV cells / JSON null).
"""

import inspect
import json
import os
import textwrap

import numpy as np
import pandas as pd


def _apply_anndata_compat_patch() -> None:
    """Keep 1-D arrays 1-D when stored in AnnData obsm/varm.

    anndata >= 0.13 reshapes every 1-D ndarray to shape (n, 1) on storage, but
    PyDESeq2 0.5.0 stores per-gene/per-sample vectors (e.g. varm["non_zero"])
    and indexes them as 1-D arrays. This shim neutralizes only that reshape
    branch, restoring the pre-0.13 behavior PyDESeq2 was built against.
    """
    from anndata._core import aligned_mapping as am

    src = textwrap.dedent(inspect.getsource(am.AlignedMappingBase._validate_value))
    old_body = "val = val.reshape((val.shape[0], 1))"
    if old_body not in src:
        raise RuntimeError(
            "anndata internals changed; 1-D reshape compatibility patch not applicable"
        )
    src = src.replace(old_body, "pass  # keep 1-D arrays 1-D (PyDESeq2 compat)")
    src = src.replace("def _validate_value", "def _validate_value_keep_1d", 1)
    ns = dict(am.__dict__)
    exec(compile(src, "<anndata_1d_compat>", "exec"), ns)
    am.AlignedMappingBase._validate_value = ns["_validate_value_keep_1d"]


_apply_anndata_compat_patch()

import pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.normpath(os.path.join(HERE, "..", "inputs"))
NUMERATOR = "Cisplatin_IC50_CBD_IC50"   # combination treatment
DENOMINATOR = "DMSO"


def main() -> None:
    # ------------------------------------------------------------------ load
    counts = pd.read_csv(os.path.join(INPUT_DIR, "counts_raw_unfiltered.csv"), index_col=0)
    layout = pd.read_csv(os.path.join(INPUT_DIR, "sample_layout.csv"))
    gene_map = pd.read_csv(os.path.join(INPUT_DIR, "ensg_to_gene_name.tsv"), sep="\t")
    gene_map = (
        gene_map.dropna(subset=["gene_name"])
        .drop_duplicates(subset=["ENSG"], keep="first")
        .set_index("ENSG")["gene_name"]
    )

    # layout SampleID uses dashes ("9-1"); counts columns use underscores ("9_1")
    layout["count_col"] = layout["SampleID"].str.replace("-", "_", regex=False)

    # ------------------------------------------- select the six frozen samples
    sel = layout[layout["Group"].isin([NUMERATOR, DENOMINATOR])].copy()
    assert len(sel) == 6, f"expected 6 samples, got {len(sel)}"
    assert set(sel["Group"]) == {NUMERATOR, DENOMINATOR}
    sample_cols = sel["count_col"].tolist()
    missing = [c for c in sample_cols if c not in counts.columns]
    assert not missing, f"count columns missing: {missing}"

    counts6 = counts.loc[:, sample_cols].astype(np.int64)
    n_genes_input = int(counts6.shape[0])

    # ----------------- prefilter: raw count > 10 in >=1 of the six samples
    keep_mask = (counts6 > 10).any(axis=1)
    counts6 = counts6.loc[keep_mask]
    n_genes_prefilter = int(counts6.shape[0])
    print(f"genes input={n_genes_input} kept_by_prefilter={n_genes_prefilter}")

    # ------------------------------------------------------------- metadata
    metadata = sel.set_index("count_col").loc[sample_cols, ["Group"]].copy()
    metadata["Group"] = pd.Categorical(
        metadata["Group"], categories=[DENOMINATOR, NUMERATOR]
    )

    # -------------------------------------------------------------- PyDESeq2
    dds = DeseqDataSet(
        counts=counts6.T,  # PyDESeq2 expects samples (rows) x genes (columns)
        metadata=metadata,
        design="~Group",   # Group is the only design term
        refit_cooks=True,
    )
    print(f"genes in DeseqDataSet={dds.n_vars} samples={dds.n_obs}")
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=["Group", NUMERATOR, DENOMINATOR])
    stat_res.summary()
    res = stat_res.results_df.copy()
    res.index.name = "gene_id"
    n_genes_tested = int(res.shape[0])
    print(f"genes tested={n_genes_tested}")

    # -------------------------------------------------------- annotate names
    res.insert(0, "gene_name", res.index.map(gene_map))

    # ------------------------------------------------------------- pass rule
    padj = res["padj"]
    passed = (
        padj.notna()
        & (padj < 0.05)
        & (res["log2FoldChange"].abs() > 0.5)
        & (res["baseMean"] > 10)
    )
    res["significant"] = passed
    res = res[["gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "significant"]]
    sig = res.loc[passed]
    n_up = int((sig["log2FoldChange"] > 0).sum())
    n_down = int((sig["log2FoldChange"] < 0).sum())
    n_null_padj = int(res["padj"].isna().sum())
    print(f"significant={len(sig)} up={n_up} down={n_down} padj_null={n_null_padj}")

    # ------------------------------------- write CSV (missing values stay empty)
    csv_path = os.path.join(HERE, "differential_expression.csv")
    res.to_csv(csv_path, index=True, index_label="gene_id")

    # ------------------------------------------------------------- summary.json
    def _f(x):
        return None if pd.isna(x) else float(x)

    summary = {
        "contrast": f"{NUMERATOR} vs {DENOMINATOR}",
        "numerator": NUMERATOR,
        "denominator": DENOMINATOR,
        "design": "Group (single term)",
        "samples": {
            NUMERATOR: sel.loc[sel["Group"] == NUMERATOR, "SampleID"].tolist(),
            DENOMINATOR: sel.loc[sel["Group"] == DENOMINATOR, "SampleID"].tolist(),
        },
        "software": {"pydeseq2": pydeseq2.__version__, "refit_cooks": True},
        "prefilter": {
            "rule": "keep gene if raw count > 10 in at least one of the six selected samples",
            "genes_input": n_genes_input,
            "genes_kept": n_genes_prefilter,
        },
        "n_genes_tested": n_genes_tested,
        "n_padj_unavailable": n_null_padj,
        "pass_criteria": {"padj": 0.05, "abs_log2FoldChange": 0.5, "baseMean": 10},
        "n_significant": int(passed.sum()),
        "n_up": n_up,
        "n_down": n_down,
        "significant_genes": [
            {
                "gene_id": idx,
                "gene_name": (None if pd.isna(row["gene_name"]) else row["gene_name"]),
                "baseMean": _f(row["baseMean"]),
                "log2FoldChange": _f(row["log2FoldChange"]),
                "padj": _f(row["padj"]),
            }
            for idx, row in sig.sort_values("padj").iterrows()
        ],
    }
    with open(os.path.join(HERE, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, allow_nan=False)

    print("wrote", csv_path)
    print("wrote", os.path.join(HERE, "summary.json"))


if __name__ == "__main__":
    main()
