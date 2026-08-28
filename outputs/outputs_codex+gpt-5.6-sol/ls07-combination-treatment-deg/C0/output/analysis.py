"""PyDESeq2 analysis for combination treatment versus DMSO."""

from __future__ import annotations

import importlib.metadata
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


TREATMENT = "Cisplatin_IC50_CBD_IC50"
CONTROL = "DMSO"
FACTOR = "Group"
EXPECTED_VERSION = "0.5.0"
PADJ_THRESHOLD = 0.05
ABS_LFC_THRESHOLD = 0.5
BASE_MEAN_THRESHOLD = 10.0


def load_selected_inputs(input_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[str]], int]:
    layout = pd.read_csv(input_dir / "sample_layout.csv", dtype=str)
    required_layout = {"SampleID", FACTOR}
    if not required_layout.issubset(layout.columns):
        raise ValueError(f"sample_layout.csv lacks {sorted(required_layout - set(layout.columns))}")

    selected = layout.loc[layout[FACTOR].isin([CONTROL, TREATMENT]), ["SampleID", FACTOR]].copy()
    group_sizes = selected.groupby(FACTOR, observed=True).size().to_dict()
    if group_sizes != {CONTROL: 3, TREATMENT: 3}:
        raise ValueError(f"Expected three replicates per target group; found {group_sizes}")
    if selected["SampleID"].duplicated().any():
        raise ValueError("Selected SampleID values are not unique")

    selected["count_column"] = selected["SampleID"].str.replace("-", "_", regex=False)
    if selected["count_column"].duplicated().any():
        raise ValueError("Selected count-column identifiers are not unique")

    all_counts = pd.read_csv(input_dir / "counts_raw_unfiltered.csv", index_col=0)
    input_gene_count = int(all_counts.shape[0])
    wanted_columns = selected["count_column"].tolist()
    missing_columns = [name for name in wanted_columns if name not in all_counts.columns]
    if missing_columns:
        raise ValueError(f"Selected samples missing from count matrix: {missing_columns}")

    counts = all_counts.loc[:, wanted_columns].apply(pd.to_numeric, errors="raise")
    values = counts.to_numpy()
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise ValueError("Selected raw counts must be finite non-negative integers")
    counts = counts.astype(np.int64)

    keep = (counts > 10).any(axis=1)
    filtered = counts.loc[keep]
    if filtered.empty:
        raise ValueError("Raw-count filter removed every gene")

    metadata = selected.set_index("count_column")[[FACTOR]].copy()
    metadata[FACTOR] = pd.Categorical(metadata[FACTOR], categories=[CONTROL, TREATMENT])
    filtered = filtered.loc[:, metadata.index]
    samples = {
        CONTROL: selected.loc[selected[FACTOR] == CONTROL, "SampleID"].tolist(),
        TREATMENT: selected.loc[selected[FACTOR] == TREATMENT, "SampleID"].tolist(),
    }
    return filtered.T, metadata, samples, input_gene_count


def run_deseq(counts: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    version = importlib.metadata.version("pydeseq2")
    if version != EXPECTED_VERSION:
        raise RuntimeError(f"Expected pydeseq2 {EXPECTED_VERSION}, found {version}")

    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design=f"~{FACTOR}",
        refit_cooks=True,
        n_cpus=1,
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(
        dds,
        contrast=[FACTOR, TREATMENT, CONTROL],
        alpha=PADJ_THRESHOLD,
        quiet=True,
    )
    stats.summary()
    results = stats.results_df.copy()
    expected = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    missing = [name for name in expected if name not in results.columns]
    if missing:
        raise ValueError(f"PyDESeq2 result lacks columns: {missing}")
    return results.loc[:, expected]


def add_mapping_and_calls(results: pd.DataFrame, input_dir: Path) -> pd.DataFrame:
    mapping = pd.read_csv(input_dir / "ensg_to_gene_name.tsv", sep="\t", dtype=str)
    if not {"ENSG", "gene_name"}.issubset(mapping.columns):
        raise ValueError("Mapping file must contain ENSG and gene_name")
    mapping = mapping.drop_duplicates("ENSG", keep="first").set_index("ENSG")["gene_name"]

    table = results.copy()
    table.insert(0, "gene_name", table.index.to_series().map(mapping))
    table.insert(0, "gene_id", table.index.astype(str))
    table["passes"] = (
        table["padj"].notna()
        & (table["padj"] < PADJ_THRESHOLD)
        & (table["log2FoldChange"].abs() > ABS_LFC_THRESHOLD)
        & (table["baseMean"] > BASE_MEAN_THRESHOLD)
    )
    table["direction"] = np.where(
        table["passes"] & (table["log2FoldChange"] > 0),
        "up",
        np.where(table["passes"] & (table["log2FoldChange"] < 0), "down", "not_passing"),
    )
    table = table.sort_values(
        ["passes", "padj", "gene_id"], ascending=[False, True, True], na_position="last"
    ).reset_index(drop=True)
    return table


def make_summary(table: pd.DataFrame, samples: dict[str, list[str]], input_genes: int) -> dict[str, object]:
    passing = table["passes"]
    return {
        "software": {"name": "PyDESeq2", "version": EXPECTED_VERSION, "refit_cooks": True},
        "design": {"formula": "~Group", "only_design_term": FACTOR},
        "contrast": {"factor": FACTOR, "numerator": TREATMENT, "denominator": CONTROL},
        "samples": samples,
        "replicates_per_group": {CONTROL: 3, TREATMENT: 3},
        "prefilter": "retain genes with raw count > 10 in at least one selected sample",
        "thresholds": {
            "padj_lt": PADJ_THRESHOLD,
            "abs_log2FoldChange_gt": ABS_LFC_THRESHOLD,
            "baseMean_gt": BASE_MEAN_THRESHOLD,
        },
        "input_gene_count": input_genes,
        "retained_gene_count": int(len(table)),
        "genes_with_unavailable_padj": int(table["padj"].isna().sum()),
        "passing_gene_count": int(passing.sum()),
        "upregulated_passing_count": int((table["direction"] == "up").sum()),
        "downregulated_passing_count": int((table["direction"] == "down").sum()),
    }


def make_report(summary: dict[str, object]) -> str:
    counts = summary
    report = f"""# Combination-treatment differential expression

## Analysis

PyDESeq2 {EXPECTED_VERSION} compared `{TREATMENT}` (numerator) with `{CONTROL}` (denominator). The model used only `Group` (`~Group`) and exactly three replicates from each group: DMSO 3-1, 3-2, 3-3 and combination treatment 9-1, 9-2, 9-3. No other layout samples entered the count matrix or metadata.

Before fitting, genes were retained when at least one of these six samples had a raw count greater than 10. This retained {counts['retained_gene_count']:,} of {counts['input_gene_count']:,} input genes. The fit used `refit_cooks=True` and the standard `DeseqStats` contrast `Group, {TREATMENT}, {CONTROL}`.

## Results

A gene passes only when all three strict criteria hold: `padj < 0.05`, `abs(log2FoldChange) > 0.5`, and `baseMean > 10`. There are {counts['passing_gene_count']:,} passing genes: {counts['upregulated_passing_count']:,} higher and {counts['downregulated_passing_count']:,} lower in combination treatment relative to DMSO. Adjusted p-values are unavailable for {counts['genes_with_unavailable_padj']:,} retained genes; these remain empty in the CSV and are counted explicitly rather than converted to zero.

`differential_expression.csv` contains all filtered-gene statistics keyed by Ensembl ID, with capsule-supplied display names where available, plus pass and direction fields.

These differential-expression results are statistical associations for this experiment. They do not by themselves establish that the combination treatment caused a particular molecular mechanism or downstream phenotype.
"""
    if len(re.findall(r"\b\w+[\w-]*\b", report)) > 500:
        raise ValueError("Report exceeds 500 words")
    return report


def main() -> None:
    output_dir = Path(__file__).resolve().parent
    repo_root = Path(__file__).resolve().parents[5]
    input_dir = repo_root / "inputs" / "ls07-combination-treatment-deg"
    counts, metadata, samples, input_genes = load_selected_inputs(input_dir)
    results = run_deseq(counts, metadata)
    table = add_mapping_and_calls(results, input_dir)
    summary = make_summary(table, samples, input_genes)

    table.to_csv(output_dir / "differential_expression.csv", index=False, na_rep="")
    (output_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    (output_dir / "report.md").write_text(make_report(summary), encoding="utf-8")


if __name__ == "__main__":
    main()
