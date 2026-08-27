"""Differential expression for combination treatment versus DMSO."""

from __future__ import annotations

import importlib.metadata
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


NUMERATOR = "Cisplatin_IC50_CBD_IC50"
DENOMINATOR = "DMSO"
DESIGN_TERM = "Group"
PYDESEQ2_VERSION = "0.5.0"
FILTER_MIN_COUNT = 10
PADJ_CUTOFF = 0.05
LFC_CUTOFF = 0.5
BASE_MEAN_CUTOFF = 10.0


def prepare_data(input_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[str]], int]:
    layout = pd.read_csv(input_dir / "sample_layout.csv", dtype=str)
    if not {"SampleID", DESIGN_TERM}.issubset(layout.columns):
        raise ValueError("sample_layout.csv must contain SampleID and Group")
    chosen = layout[layout[DESIGN_TERM].isin([DENOMINATOR, NUMERATOR])][
        ["SampleID", DESIGN_TERM]
    ].copy()
    observed = chosen[DESIGN_TERM].value_counts().to_dict()
    if observed != {DENOMINATOR: 3, NUMERATOR: 3}:
        raise ValueError(f"Target groups must have exactly three replicates: {observed}")
    if len(chosen) != 6 or chosen["SampleID"].nunique() != 6:
        raise ValueError("The selected model must contain six unique samples")

    chosen["matrix_id"] = chosen["SampleID"].str.replace("-", "_", regex=False)
    raw = pd.read_csv(input_dir / "counts_raw_unfiltered.csv", index_col=0)
    total_genes = int(raw.shape[0])
    matrix_ids = chosen["matrix_id"].tolist()
    absent = sorted(set(matrix_ids) - set(raw.columns))
    if absent:
        raise ValueError(f"Selected count columns are absent: {absent}")

    six_counts = raw[matrix_ids].apply(pd.to_numeric, errors="raise")
    array = six_counts.to_numpy()
    if (array < 0).any() or not np.isfinite(array).all() or not np.equal(array, np.floor(array)).all():
        raise ValueError("Selected counts must be finite non-negative raw integers")
    six_counts = six_counts.astype("int64")
    retained = six_counts[(six_counts > FILTER_MIN_COUNT).any(axis=1)]

    metadata = chosen.set_index("matrix_id")[[DESIGN_TERM]]
    metadata[DESIGN_TERM] = pd.Categorical(
        metadata[DESIGN_TERM], categories=[DENOMINATOR, NUMERATOR]
    )
    retained = retained.loc[:, metadata.index]
    sample_ids = {
        DENOMINATOR: chosen.loc[chosen[DESIGN_TERM] == DENOMINATOR, "SampleID"].tolist(),
        NUMERATOR: chosen.loc[chosen[DESIGN_TERM] == NUMERATOR, "SampleID"].tolist(),
    }
    return retained.T, metadata, sample_ids, total_genes


def fit(counts: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    actual_version = importlib.metadata.version("pydeseq2")
    if actual_version != PYDESEQ2_VERSION:
        raise RuntimeError(f"Expected pydeseq2 {PYDESEQ2_VERSION}, found {actual_version}")
    dataset = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design="~Group",
        refit_cooks=True,
        n_cpus=1,
        quiet=True,
    )
    dataset.deseq2()
    statistics = DeseqStats(
        dataset,
        contrast=[DESIGN_TERM, NUMERATOR, DENOMINATOR],
        alpha=PADJ_CUTOFF,
        quiet=True,
    )
    statistics.summary()
    columns = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    missing = [column for column in columns if column not in statistics.results_df.columns]
    if missing:
        raise ValueError(f"Missing DE result columns: {missing}")
    return statistics.results_df[columns].copy()


def annotate(results: pd.DataFrame, input_dir: Path) -> pd.DataFrame:
    names = pd.read_csv(input_dir / "ensg_to_gene_name.tsv", sep="\t", dtype=str)
    if not {"ENSG", "gene_name"}.issubset(names.columns):
        raise ValueError("Mapping input lacks ENSG/gene_name")
    lookup = names.drop_duplicates("ENSG").set_index("ENSG")["gene_name"]

    output = results.copy()
    output.insert(0, "gene_name", output.index.to_series().map(lookup))
    output.insert(0, "gene_id", output.index.astype(str))
    output["passes"] = (
        output["padj"].notna()
        & output["padj"].lt(PADJ_CUTOFF)
        & output["log2FoldChange"].abs().gt(LFC_CUTOFF)
        & output["baseMean"].gt(BASE_MEAN_CUTOFF)
    )
    output["direction"] = "not_passing"
    output.loc[output["passes"] & output["log2FoldChange"].gt(0), "direction"] = "up"
    output.loc[output["passes"] & output["log2FoldChange"].lt(0), "direction"] = "down"
    return output.sort_values(
        ["passes", "padj", "gene_id"],
        ascending=[False, True, True],
        na_position="last",
    ).reset_index(drop=True)


def summarize(table: pd.DataFrame, samples: dict[str, list[str]], total_genes: int) -> dict[str, object]:
    return {
        "software": {"name": "PyDESeq2", "version": PYDESEQ2_VERSION, "refit_cooks": True},
        "design": {"formula": "~Group", "only_design_term": DESIGN_TERM},
        "contrast": {"factor": DESIGN_TERM, "numerator": NUMERATOR, "denominator": DENOMINATOR},
        "samples": samples,
        "replicates_per_group": {DENOMINATOR: 3, NUMERATOR: 3},
        "prefilter": "retain genes with raw count > 10 in at least one selected sample",
        "thresholds": {
            "padj_lt": PADJ_CUTOFF,
            "abs_log2FoldChange_gt": LFC_CUTOFF,
            "baseMean_gt": BASE_MEAN_CUTOFF,
        },
        "input_gene_count": total_genes,
        "retained_gene_count": int(len(table)),
        "genes_with_unavailable_padj": int(table["padj"].isna().sum()),
        "passing_gene_count": int(table["passes"].sum()),
        "upregulated_passing_count": int((table["direction"] == "up").sum()),
        "downregulated_passing_count": int((table["direction"] == "down").sum()),
        "external_transcript_queries_used": False,
    }


def report(summary: dict[str, object]) -> str:
    text = f"""# Combination treatment versus DMSO: differential expression

## Scope and method

The analysis used only the six requested samples: DMSO 3-1, 3-2, 3-3 and `{NUMERATOR}` 9-1, 9-2, 9-3. PyDESeq2 {PYDESEQ2_VERSION} fitted `~Group` with `refit_cooks=True`; the standard contrast used `{NUMERATOR}` as numerator and DMSO as denominator. No other condition, covariate, transcript lookup, sequence, haplotype, or genome-track data entered the model.

Genes were retained before fitting when any selected sample had a raw count greater than 10. This retained {summary['retained_gene_count']:,} of {summary['input_gene_count']:,} genes.

## Thresholded findings

Passing genes simultaneously satisfy `padj < 0.05`, `abs(log2FoldChange) > 0.5`, and `baseMean > 10`. A total of {summary['passing_gene_count']:,} genes pass: {summary['upregulated_passing_count']:,} up and {summary['downregulated_passing_count']:,} down in combination treatment relative to DMSO. For {summary['genes_with_unavailable_padj']:,} retained genes, `padj` is unavailable and remains an empty CSV cell, never zero.

The complete filtered result table is keyed by Ensembl ID; the supplied mapping contributes display names without dropping unmapped genes.

These estimates identify statistical associations between treatment group and expression in this experiment. They do not establish that the combination treatment caused a specific mechanism or phenotype.
"""
    if len(re.findall(r"\b\w+[\w-]*\b", text)) > 500:
        raise ValueError("Report exceeds 500 words")
    return text


def main() -> None:
    output_dir = Path(__file__).resolve().parent
    repo_root = Path(__file__).resolve().parents[5]
    input_dir = repo_root / "inputs" / "ls07-combination-treatment-deg"
    counts, metadata, samples, total_genes = prepare_data(input_dir)
    table = annotate(fit(counts, metadata), input_dir)
    summary = summarize(table, samples, total_genes)
    table.to_csv(output_dir / "differential_expression.csv", index=False, na_rep="")
    (output_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    (output_dir / "report.md").write_text(report(summary), encoding="utf-8")


if __name__ == "__main__":
    main()
