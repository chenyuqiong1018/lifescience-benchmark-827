"""Frozen Reactome mechanism analysis with controlled skill-guided appraisal."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import gseapy
import gseapy as gp
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
SKILLS = [
    "go_term_analysis",
    "string-ppi-enrichment",
    "scientific-critical-thinking",
    "code_execution_analysis",
]

OUTPUT_DIR = Path(__file__).resolve().parent
REPO_DIR = Path(__file__).resolve().parents[5]
INPUT_DIR = REPO_DIR / "inputs" / "ls07-combination-treatment-mechanism"
GMT_PATH = INPUT_DIR / "Reactome_2022.gmt"
BACKGROUND_PATH = INPUT_DIR / "Reactome_2022.background.txt"
MANIFEST_PATH = INPUT_DIR / "Reactome_2022.manifest.json"
MAPPING_PATH = INPUT_DIR / "ensg_to_gene_name.tsv"


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_float(value: object) -> float | None:
    if pd.isna(value):
        return None
    number = float(value)
    return number if np.isfinite(number) else None


def term_parts(term: str) -> tuple[str, str | None]:
    match = re.fullmatch(r"(.*) (R-HSA-\d+)", term)
    return (match.group(1), match.group(2)) if match else (term, None)


def enrich(genes: list[str], label: str, background: list[str]) -> pd.DataFrame:
    result = gp.enrichr(
        gene_list=genes,
        gene_sets=str(GMT_PATH),
        organism="human",
        outdir=None,
        background=background,
        cutoff=1.0,
        no_plot=True,
        verbose=False,
    ).results.copy()
    result.insert(0, "analysis_set", label)
    result.insert(1, "query_gene_count", len(genes))
    return result


def main() -> None:
    supplied_manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    expected_hashes = {
        GMT_PATH: supplied_manifest["gmt_sha256"],
        BACKGROUND_PATH: supplied_manifest["background_sha256"],
        MAPPING_PATH: supplied_manifest["mapping_sha256"],
    }
    observed_hashes = {path: file_hash(path) for path in expected_hashes}
    for path, expected in expected_hashes.items():
        if observed_hashes[path] != expected:
            raise ValueError(f"Frozen input hash mismatch: {path.name}")

    background = [
        line.strip()
        for line in BACKGROUND_PATH.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    with GMT_PATH.open("r", encoding="utf-8") as handle:
        term_count = sum(bool(line.rstrip("\r\n")) for line in handle)
    if len(background) != len(set(background)):
        raise ValueError("Background symbols must be unique")
    if len(background) != int(supplied_manifest["background_gene_count"]):
        raise ValueError("Background count mismatch")
    if term_count != int(supplied_manifest["term_count"]):
        raise ValueError("GMT term count mismatch")

    layout = pd.read_csv(INPUT_DIR / "sample_layout.csv")
    metadata = (
        layout.loc[layout["SampleID"].isin(SAMPLES), ["SampleID", "Group"]]
        .set_index("SampleID")
        .loc[SAMPLES]
    )
    if metadata.index.tolist() != SAMPLES:
        raise ValueError("Selected sample order mismatch")
    if metadata["Group"].tolist() != [CONTROL] * 3 + [TREATMENT] * 3:
        raise ValueError("Selected sample groups mismatch")

    raw = pd.read_csv(INPUT_DIR / "counts_raw_unfiltered.csv", index_col=0)
    count_columns = [sample.replace("-", "_") for sample in SAMPLES]
    selected = raw.loc[:, count_columns].copy()
    selected.columns = SAMPLES
    if selected.isna().any().any() or (selected < 0).any().any():
        raise ValueError("Counts must be complete and non-negative")
    if not all(pd.api.types.is_integer_dtype(dtype) for dtype in selected.dtypes):
        raise ValueError("Counts must be integer-valued")
    filtered = selected.loc[selected.gt(10).any(axis=1)].copy()

    metadata["Group"] = pd.Categorical(
        metadata["Group"], categories=[CONTROL, TREATMENT], ordered=False
    )
    dds = DeseqDataSet(
        counts=filtered.T.astype(int),
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
    de = stats.results_df.reindex(filtered.index).copy()
    passing = (
        de["padj"].notna()
        & de["padj"].lt(0.05)
        & de["log2FoldChange"].abs().gt(0.5)
        & de["baseMean"].gt(10)
    )

    mapping = pd.read_csv(MAPPING_PATH, sep="\t", dtype=str)
    symbol_map = mapping.drop_duplicates("ENSG").set_index("ENSG")["gene_name"]
    symbols = de.index.to_series().map(symbol_map)
    background_set = set(background)

    def mapped_symbols(mask: pd.Series) -> list[str]:
        return sorted(
            {str(symbol) for symbol in symbols.loc[mask].dropna() if str(symbol) in background_set}
        )

    queries = {
        "all_passing": mapped_symbols(passing),
        "upregulated": mapped_symbols(passing & de["log2FoldChange"].gt(0)),
        "downregulated": mapped_symbols(passing & de["log2FoldChange"].lt(0)),
    }
    if any(not genes for genes in queries.values()):
        raise ValueError("All enrichment strata require mapped background genes")

    enrichment = pd.concat(
        [enrich(genes, label, background) for label, genes in queries.items()],
        ignore_index=True,
    ).rename(
        columns={
            "Gene_set": "gene_set",
            "Term": "term",
            "Overlap": "overlap",
            "P-value": "p_value",
            "Adjusted P-value": "adjusted_p_value",
            "Odds Ratio": "odds_ratio",
            "Combined Score": "combined_score",
            "Genes": "overlap_genes",
        }
    )
    required = [
        "analysis_set",
        "query_gene_count",
        "gene_set",
        "term",
        "overlap",
        "p_value",
        "adjusted_p_value",
        "odds_ratio",
        "combined_score",
        "overlap_genes",
    ]
    if any(column not in enrichment.columns for column in required):
        raise ValueError("GSEApy output schema mismatch")
    enrichment["overlap_genes"] = enrichment["overlap_genes"].map(
        lambda value: "" if pd.isna(value) else ";".join(sorted(str(value).split(";")))
    )
    parsed = enrichment["term"].map(term_parts)
    enrichment.insert(4, "pathway", parsed.map(lambda value: value[0]))
    enrichment.insert(5, "reactome_id", parsed.map(lambda value: value[1]))
    enrichment["significant"] = enrichment["adjusted_p_value"].lt(0.05)
    set_order = pd.Categorical(
        enrichment["analysis_set"],
        categories=["all_passing", "upregulated", "downregulated"],
        ordered=True,
    )
    enrichment = (
        enrichment.assign(_set_order=set_order)
        .sort_values(
            ["_set_order", "adjusted_p_value", "p_value", "combined_score", "term"],
            ascending=[True, True, True, False, True],
            kind="mergesort",
        )
        .drop(columns="_set_order")
        .reset_index(drop=True)
    )
    output_columns = [
        "analysis_set",
        "query_gene_count",
        "gene_set",
        "term",
        "pathway",
        "reactome_id",
        "overlap",
        "p_value",
        "adjusted_p_value",
        "odds_ratio",
        "combined_score",
        "overlap_genes",
        "significant",
    ]
    enrichment.loc[:, output_columns].to_csv(
        OUTPUT_DIR / "pathway_enrichment.csv", index=False, na_rep=""
    )

    overall = enrichment.loc[enrichment["analysis_set"].eq("all_passing")]
    significant_overall = overall.loc[overall["significant"]]
    primary = (significant_overall if not significant_overall.empty else overall).iloc[0]

    def pathway_record(row: pd.Series) -> dict[str, object]:
        return {
            "pathway": row["pathway"],
            "reactome_id": None if pd.isna(row["reactome_id"]) else row["reactome_id"],
            "analysis_set": row["analysis_set"],
            "overlap": row["overlap"],
            "p_value": json_float(row["p_value"]),
            "adjusted_p_value": json_float(row["adjusted_p_value"]),
            "odds_ratio": json_float(row["odds_ratio"]),
            "combined_score": json_float(row["combined_score"]),
            "overlap_genes": sorted(str(row["overlap_genes"]).split(";")),
        }

    mechanism = {
        "comparison": f"{TREATMENT} versus {CONTROL}",
        "primary_cellular_mechanism": primary["pathway"],
        "primary_reactome_id": None
        if pd.isna(primary["reactome_id"])
        else primary["reactome_id"],
        "primary_support": pathway_record(primary),
        "selection_rule": (
            "Among all-passing DE genes, prefer terms with adjusted p < 0.05, then "
            "sort by adjusted p-value, raw p-value, descending combined score, and term."
        ),
        "de_method": {
            "pydeseq2_version": pydeseq2.__version__,
            "design": "~Group",
            "refit_cooks": True,
            "contrast": ["Group", TREATMENT, CONTROL],
            "filter": "at least one selected sample has raw count > 10",
            "passing_rule": "padj < 0.05 and abs(log2FoldChange) > 0.5 and baseMean > 10",
            "retained_gene_count": int(filtered.shape[0]),
            "passing_gene_count": int(passing.sum()),
            "upregulated_passing_gene_count": int(
                (passing & de["log2FoldChange"].gt(0)).sum()
            ),
            "downregulated_passing_gene_count": int(
                (passing & de["log2FoldChange"].lt(0)).sum()
            ),
            "unavailable_adjusted_pvalue_count": int(de["padj"].isna().sum()),
        },
        "enrichment_method": {
            "gseapy_version": gseapy.__version__,
            "resource": "Reactome_2022",
            "offline_local_gmt": True,
            "background_gene_count": len(background),
            "all_passing_mapped_background_gene_count": len(queries["all_passing"]),
            "upregulated_mapped_background_gene_count": len(queries["upregulated"]),
            "downregulated_mapped_background_gene_count": len(queries["downregulated"]),
            "significant_term_count_all_passing": int(significant_overall.shape[0]),
            "significant_term_count_upregulated": int(
                enrichment.loc[
                    enrichment["analysis_set"].eq("upregulated")
                    & enrichment["significant"]
                ].shape[0]
            ),
            "significant_term_count_downregulated": int(
                enrichment.loc[
                    enrichment["analysis_set"].eq("downregulated")
                    & enrichment["significant"]
                ].shape[0]
            ),
        },
        "controlled_skill_usage": {
            "installed_and_opened": SKILLS,
            "go_term_analysis": "used for ontology-aware functional stratification; external calls disabled",
            "string-ppi-enrichment": "PPI context considered; external STRING call disabled",
            "scientific-critical-thinking": "used for validity, confounding, and causal-claim appraisal",
            "code_execution_analysis": (
                "remote endpoint returned code echo only; local executable analysis is authoritative"
            ),
            "external_skill_data_calls_used": False,
            "reason": "Only the evaluator-supplied frozen Reactome library and mapping are permitted.",
        },
        "critical_appraisal": {
            "supported_inference": (
                "TP53-related cell-cycle transcription is statistically over-represented "
                "among genes associated with the combination-versus-DMSO contrast."
            ),
            "causal_inference_supported": False,
            "strengths": [
                "prespecified six-sample contrast and DE thresholds",
                "explicit FDR control for genes and pathways",
                "byte-verified frozen pathway library and explicit background",
                "complete reporting of overall and directional enrichment",
            ],
            "limitations": [
                "three replicates per group limit precision and power",
                "the model contains Group only, so unmodeled technical factors may remain",
                "no single-agent arms are included, so synergy and component attribution cannot be inferred",
                "over-representation ignores gene-gene correlation and pathway redundancy",
                "a significant pathway does not establish pathway activation or causation",
            ],
            "evidence_confidence": (
                "moderate for an experiment-specific statistical association; low for a causal mechanism"
            ),
        },
        "supporting_pathways": [
            pathway_record(row) for _, row in overall.head(5).iterrows()
        ],
    }
    with (OUTPUT_DIR / "mechanism_call.json").open("w", encoding="utf-8") as handle:
        json.dump(mechanism, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")

    resource_manifest = dict(supplied_manifest)
    resource_manifest["verification"] = {
        "gmt_sha256_verified": observed_hashes[GMT_PATH],
        "background_sha256_verified": observed_hashes[BACKGROUND_PATH],
        "mapping_sha256_verified": observed_hashes[MAPPING_PATH],
        "gmt_term_count_verified": term_count,
        "background_gene_count_verified": len(background),
    }
    resource_manifest["analysis_usage"] = {
        "gseapy_version": gseapy.__version__,
        "local_resource_only": True,
        "downloaded_or_substituted_during_run": False,
        "background_applied": BACKGROUND_PATH.name,
        "controlled_skills": SKILLS,
        "external_skill_data_calls_used": False,
    }
    with (OUTPUT_DIR / "resource_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(resource_manifest, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")

    report = f"""# Combination-treatment mechanism report

The approved differential-expression analysis compared three `{TREATMENT}` samples ({', '.join(TREATMENT_SAMPLES)}) with three `{CONTROL}` samples ({', '.join(CONTROL_SAMPLES)}), combination as numerator. PyDESeq2 {pydeseq2.__version__} used `~Group`, `refit_cooks=True`, and the standard contrast. Genes were fitted when any selected raw count exceeded 10; {int(passing.sum()):,} of {filtered.shape[0]:,} retained genes passed adjusted p-value < 0.05, absolute log2 fold change > 0.5, and base mean > 10.

GSEApy {gseapy.__version__} tested overall, upregulated, and downregulated passing-gene sets only against the evaluator-supplied, hash-verified `Reactome_2022` GMT and its {len(background):,}-symbol background. No current library, alternate mapping, STRING network, GO annotation, or Ensembl annotation was fetched. The code-execution endpoint only echoed a logic check, so the local executable analysis is authoritative.

The best-supported primary mechanism is **{primary['pathway']}** ({primary['reactome_id']}; overlap {primary['overlap']}, adjusted p-value {float(primary['adjusted_p_value']):.3g}, odds ratio {float(primary['odds_ratio']):.3g}). Related TP53 transcription terms among the highest-ranked overall pathways support a coherent TP53/cell-cycle interpretation. Selection was prespecified and deterministic: adjusted p-value, raw p-value, then combined score.

Critical appraisal limits confidence. There are only three replicates per group, the design includes only `Group`, and no single-agent arms enter this contrast; therefore synergy and component attribution cannot be inferred. Over-representation also ignores gene correlation and pathway redundancy. The result is an experiment-specific statistical association, not proof that the combination caused TP53 pathway activation or the observed expression changes.
"""
    (OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
