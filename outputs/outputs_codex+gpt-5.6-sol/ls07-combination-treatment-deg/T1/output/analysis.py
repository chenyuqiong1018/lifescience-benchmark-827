"""Differential expression for combination treatment versus DMSO."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


TREATMENT = "Cisplatin_IC50_CBD_IC50"
CONTROL = "DMSO"
CONTROL_SAMPLES = ["3-1", "3-2", "3-3"]
TREATMENT_SAMPLES = ["9-1", "9-2", "9-3"]
SAMPLES = CONTROL_SAMPLES + TREATMENT_SAMPLES

OUTPUT_DIR = Path(__file__).resolve().parent
REPO_DIR = Path(__file__).resolve().parents[5]
INPUT_DIR = REPO_DIR / "inputs" / "ls07-combination-treatment-deg"


def json_number(value: object) -> float | None:
    """Return finite numeric values and encode unavailable values as JSON null."""
    if pd.isna(value):
        return None
    number = float(value)
    return number if np.isfinite(number) else None


def main() -> None:
    layout = pd.read_csv(INPUT_DIR / "sample_layout.csv")
    chosen = layout.loc[layout["SampleID"].isin(SAMPLES), ["SampleID", "Group"]].copy()
    chosen = chosen.set_index("SampleID").loc[SAMPLES]

    expected_groups = [CONTROL] * 3 + [TREATMENT] * 3
    if chosen.index.tolist() != SAMPLES or chosen["Group"].tolist() != expected_groups:
        raise ValueError("Selected samples or groups do not match the prespecified six-sample contrast")

    raw = pd.read_csv(INPUT_DIR / "counts_raw_unfiltered.csv", index_col=0)
    count_columns = [sample.replace("-", "_") for sample in SAMPLES]
    missing_columns = sorted(set(count_columns).difference(raw.columns))
    if missing_columns:
        raise ValueError(f"Missing count columns: {missing_columns}")

    selected = raw.loc[:, count_columns].copy()
    selected.columns = SAMPLES
    if selected.isna().any().any() or (selected < 0).any().any():
        raise ValueError("Counts must be complete and non-negative")
    if not all(pd.api.types.is_integer_dtype(dtype) for dtype in selected.dtypes):
        raise ValueError("Counts must be integer-valued")

    keep = selected.gt(10).any(axis=1)
    filtered = selected.loc[keep].copy()
    counts = filtered.T.astype(int)
    metadata = chosen.copy()
    metadata["Group"] = pd.Categorical(
        metadata["Group"], categories=[CONTROL, TREATMENT], ordered=False
    )

    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design="~Group",
        refit_cooks=True,
        n_cpus=1,
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(
        dds,
        contrast=["Group", TREATMENT, CONTROL],
        alpha=0.05,
        quiet=True,
    )
    stats.summary()

    results = stats.results_df.reindex(filtered.index).copy()
    required = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    if any(column not in results.columns for column in required):
        raise ValueError("PyDESeq2 result table is missing required columns")

    mapping = pd.read_csv(INPUT_DIR / "ensg_to_gene_name.tsv", sep="\t", dtype=str)
    gene_name = mapping.drop_duplicates("ENSG").set_index("ENSG")["gene_name"]

    passes = (
        results["padj"].notna()
        & results["padj"].lt(0.05)
        & results["log2FoldChange"].abs().gt(0.5)
        & results["baseMean"].gt(10)
    )
    direction = np.where(
        passes & results["log2FoldChange"].gt(0),
        "up",
        np.where(passes & results["log2FoldChange"].lt(0), "down", "not_passing"),
    )

    output = results.loc[:, required].copy()
    output.insert(0, "gene_name", output.index.to_series().map(gene_name).fillna(""))
    output.insert(0, "gene_id", output.index)
    output["passes"] = passes.to_numpy(dtype=bool)
    output["direction"] = direction
    output.to_csv(OUTPUT_DIR / "differential_expression.csv", index=False, na_rep="")

    pass_count = int(passes.sum())
    up_count = int((passes & results["log2FoldChange"].gt(0)).sum())
    down_count = int((passes & results["log2FoldChange"].lt(0)).sum())
    unavailable_padj = int(results["padj"].isna().sum())
    summary = {
        "comparison": f"{TREATMENT} versus {CONTROL}",
        "numerator": TREATMENT,
        "denominator": CONTROL,
        "design": "~Group",
        "design_terms": ["Group"],
        "control_samples": CONTROL_SAMPLES,
        "treatment_samples": TREATMENT_SAMPLES,
        "sample_count": 6,
        "input_gene_count": int(raw.shape[0]),
        "filter": "at least one selected sample has raw count > 10",
        "retained_gene_count": int(filtered.shape[0]),
        "pydeseq2_version": pydeseq2.__version__,
        "refit_cooks": True,
        "contrast": ["Group", TREATMENT, CONTROL],
        "threshold_rule": "padj < 0.05 and abs(log2FoldChange) > 0.5 and baseMean > 10",
        "passing_gene_count": pass_count,
        "upregulated_passing_gene_count": up_count,
        "downregulated_passing_gene_count": down_count,
        "unavailable_adjusted_pvalue_count": unavailable_padj,
        "unavailable_adjusted_pvalues_encoded_as": None,
        "external_biomarker_queries_used": False,
        "external_transcript_queries_used": False,
    }
    with (OUTPUT_DIR / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")

    report = f"""# Differential-expression report

PyDESeq2 {pydeseq2.__version__} was used to compare `{TREATMENT}` (numerator; samples {', '.join(TREATMENT_SAMPLES)}) with `{CONTROL}` (denominator; samples {', '.join(CONTROL_SAMPLES)}). No other groups were used. The sole design term was `Group` (`~Group`), `refit_cooks=True`, and the standard `DeseqStats` contrast was `['Group', '{TREATMENT}', '{CONTROL}']`.

Before model fitting, genes were retained when at least one of the six selected samples had a raw count greater than 10. This retained {filtered.shape[0]:,} of {raw.shape[0]:,} input genes. A gene passed only when adjusted p-value < 0.05, absolute log2 fold change > 0.5, and base mean > 10. {pass_count:,} genes passed: {up_count:,} had positive and {down_count:,} had negative log2 fold change. Adjusted p-values were unavailable for {unavailable_padj:,} retained genes and remain empty in the CSV and null in the JSON summary; they were not treated as zero.

The supplied Ensembl-to-gene-name table was used only to add labels; statistical rows remain keyed by Ensembl gene ID. No external biomarker or transcript database was queried. These results identify expression changes statistically associated with the treatment contrast under this experiment. They do not establish that the combination treatment causally regulates any individual gene, nor do they isolate interaction effects between the two agents.
"""
    (OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
