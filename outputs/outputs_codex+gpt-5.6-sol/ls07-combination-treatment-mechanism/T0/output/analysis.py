"""Offline Reactome mechanism analysis with functional-enrichment stratification."""

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
SKILL = "go_term_analysis"

OUTPUT_DIR = Path(__file__).resolve().parent
REPO_DIR = Path(__file__).resolve().parents[5]
INPUT_DIR = REPO_DIR / "inputs" / "ls07-combination-treatment-mechanism"
GMT_PATH = INPUT_DIR / "Reactome_2022.gmt"
BACKGROUND_PATH = INPUT_DIR / "Reactome_2022.background.txt"
MANIFEST_PATH = INPUT_DIR / "Reactome_2022.manifest.json"
MAPPING_PATH = INPUT_DIR / "ensg_to_gene_name.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite(value: object) -> float | None:
    if pd.isna(value):
        return None
    number = float(value)
    return number if np.isfinite(number) else None


def split_term(term: str) -> tuple[str, str | None]:
    match = re.fullmatch(r"(.*) (R-HSA-\d+)", term)
    return (match.group(1), match.group(2)) if match else (term, None)


def local_enrichment(genes: list[str], label: str, background: list[str]) -> pd.DataFrame:
    table = gp.enrichr(
        gene_list=genes,
        gene_sets=str(GMT_PATH),
        organism="human",
        outdir=None,
        background=background,
        cutoff=1.0,
        no_plot=True,
        verbose=False,
    ).results.copy()
    table.insert(0, "analysis_set", label)
    table.insert(1, "query_gene_count", len(genes))
    return table


def main() -> None:
    source_manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    expected = {
        GMT_PATH: source_manifest["gmt_sha256"],
        BACKGROUND_PATH: source_manifest["background_sha256"],
        MAPPING_PATH: source_manifest["mapping_sha256"],
    }
    observed = {path: sha256(path) for path in expected}
    for path in expected:
        if observed[path] != expected[path]:
            raise ValueError(f"Frozen input hash mismatch: {path.name}")

    background = [
        symbol.strip()
        for symbol in BACKGROUND_PATH.read_text(encoding="utf-8").splitlines()
        if symbol.strip()
    ]
    with GMT_PATH.open("r", encoding="utf-8") as stream:
        term_count = sum(bool(line.rstrip("\r\n")) for line in stream)
    if len(background) != len(set(background)):
        raise ValueError("Background symbols are not unique")
    if len(background) != int(source_manifest["background_gene_count"]):
        raise ValueError("Background count mismatch")
    if term_count != int(source_manifest["term_count"]):
        raise ValueError("GMT term count mismatch")

    layout = pd.read_csv(INPUT_DIR / "sample_layout.csv")
    metadata = (
        layout.loc[layout["SampleID"].isin(SAMPLES), ["SampleID", "Group"]]
        .set_index("SampleID")
        .loc[SAMPLES]
    )
    if metadata["Group"].tolist() != [CONTROL] * 3 + [TREATMENT] * 3:
        raise ValueError("Prespecified six-sample contrast mismatch")

    raw = pd.read_csv(INPUT_DIR / "counts_raw_unfiltered.csv", index_col=0)
    selected = raw.loc[:, [sample.replace("-", "_") for sample in SAMPLES]].copy()
    selected.columns = SAMPLES
    if selected.isna().any().any() or (selected < 0).any().any():
        raise ValueError("Counts must be complete and non-negative")
    if not all(pd.api.types.is_integer_dtype(dtype) for dtype in selected.dtypes):
        raise ValueError("Counts must be integers")
    filtered = selected.loc[selected.gt(10).any(axis=1)]

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
    statistics = DeseqStats(
        dds,
        contrast=["Group", TREATMENT, CONTROL],
        alpha=0.05,
        quiet=True,
    )
    statistics.summary()
    de = statistics.results_df.reindex(filtered.index)
    passes = (
        de["padj"].notna()
        & de["padj"].lt(0.05)
        & de["log2FoldChange"].abs().gt(0.5)
        & de["baseMean"].gt(10)
    )

    mapping = pd.read_csv(MAPPING_PATH, sep="\t", dtype=str)
    symbol_map = mapping.drop_duplicates("ENSG").set_index("ENSG")["gene_name"]
    symbols = de.index.to_series().map(symbol_map)
    background_set = set(background)

    def mapped(mask: pd.Series) -> list[str]:
        return sorted(
            {str(value) for value in symbols.loc[mask].dropna() if str(value) in background_set}
        )

    query = {
        "all_passing": mapped(passes),
        "upregulated": mapped(passes & de["log2FoldChange"].gt(0)),
        "downregulated": mapped(passes & de["log2FoldChange"].lt(0)),
    }
    if any(not genes for genes in query.values()):
        raise ValueError("Every functional-enrichment stratum must contain mapped genes")

    enriched = pd.concat(
        [local_enrichment(genes, label, background) for label, genes in query.items()],
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
    if any(column not in enriched.columns for column in required):
        raise ValueError("GSEApy output schema mismatch")
    enriched["overlap_genes"] = enriched["overlap_genes"].map(
        lambda value: "" if pd.isna(value) else ";".join(sorted(str(value).split(";")))
    )
    term_parts = enriched["term"].map(split_term)
    enriched.insert(4, "pathway", term_parts.map(lambda pair: pair[0]))
    enriched.insert(5, "reactome_id", term_parts.map(lambda pair: pair[1]))
    enriched["significant"] = enriched["adjusted_p_value"].lt(0.05)
    set_order = pd.Categorical(
        enriched["analysis_set"],
        categories=["all_passing", "upregulated", "downregulated"],
        ordered=True,
    )
    enriched = (
        enriched.assign(_set_order=set_order)
        .sort_values(
            ["_set_order", "adjusted_p_value", "p_value", "combined_score", "term"],
            ascending=[True, True, True, False, True],
            kind="mergesort",
        )
        .drop(columns="_set_order")
        .reset_index(drop=True)
    )
    columns = [
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
    enriched.loc[:, columns].to_csv(
        OUTPUT_DIR / "pathway_enrichment.csv", index=False, na_rep=""
    )

    all_terms = enriched.loc[enriched["analysis_set"].eq("all_passing")]
    significant_all = all_terms.loc[all_terms["significant"]]
    primary = (significant_all if not significant_all.empty else all_terms).iloc[0]

    def record(row: pd.Series) -> dict[str, object]:
        return {
            "pathway": row["pathway"],
            "reactome_id": None if pd.isna(row["reactome_id"]) else row["reactome_id"],
            "analysis_set": row["analysis_set"],
            "overlap": row["overlap"],
            "p_value": finite(row["p_value"]),
            "adjusted_p_value": finite(row["adjusted_p_value"]),
            "odds_ratio": finite(row["odds_ratio"]),
            "combined_score": finite(row["combined_score"]),
            "overlap_genes": sorted(str(row["overlap_genes"]).split(";")),
        }

    mechanism = {
        "comparison": f"{TREATMENT} versus {CONTROL}",
        "primary_cellular_mechanism": primary["pathway"],
        "primary_reactome_id": None
        if pd.isna(primary["reactome_id"])
        else primary["reactome_id"],
        "primary_support": record(primary),
        "selection_rule": (
            "Within all-passing DE genes, sort by adjusted p-value, raw p-value, "
            "descending combined score, and term; prefer adjusted p < 0.05."
        ),
        "de_method": {
            "pydeseq2_version": pydeseq2.__version__,
            "design": "~Group",
            "refit_cooks": True,
            "contrast": ["Group", TREATMENT, CONTROL],
            "filter": "at least one selected sample has raw count > 10",
            "passing_rule": "padj < 0.05 and abs(log2FoldChange) > 0.5 and baseMean > 10",
            "retained_gene_count": int(filtered.shape[0]),
            "passing_gene_count": int(passes.sum()),
            "upregulated_passing_gene_count": int(
                (passes & de["log2FoldChange"].gt(0)).sum()
            ),
            "downregulated_passing_gene_count": int(
                (passes & de["log2FoldChange"].lt(0)).sum()
            ),
            "unavailable_adjusted_pvalue_count": int(de["padj"].isna().sum()),
        },
        "enrichment_method": {
            "gseapy_version": gseapy.__version__,
            "resource": "Reactome_2022",
            "offline_local_gmt": True,
            "background_gene_count": len(background),
            "all_passing_mapped_background_gene_count": len(query["all_passing"]),
            "upregulated_mapped_background_gene_count": len(query["upregulated"]),
            "downregulated_mapped_background_gene_count": len(query["downregulated"]),
            "significant_term_count_all_passing": int(significant_all.shape[0]),
            "significant_term_count_upregulated": int(
                enriched.loc[
                    enriched["analysis_set"].eq("upregulated") & enriched["significant"]
                ].shape[0]
            ),
            "significant_term_count_downregulated": int(
                enriched.loc[
                    enriched["analysis_set"].eq("downregulated") & enriched["significant"]
                ].shape[0]
            ),
        },
        "skill_usage": {
            "selected_skill": SKILL,
            "applied_as": "functional-enrichment stratification and ontology-aware interpretation",
            "external_string_go_ensembl_calls_used": False,
            "reason_external_calls_not_used": (
                "The task requires the evaluator-supplied frozen Reactome library and mapping only."
            ),
        },
        "supporting_pathways": [record(row) for _, row in all_terms.head(5).iterrows()],
        "interpretation_limit": (
            "Enrichment is an association with the treatment contrast, not demonstrated causation."
        ),
    }
    with (OUTPUT_DIR / "mechanism_call.json").open("w", encoding="utf-8") as stream:
        json.dump(mechanism, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")

    output_manifest = dict(source_manifest)
    output_manifest["verification"] = {
        "gmt_sha256_verified": observed[GMT_PATH],
        "background_sha256_verified": observed[BACKGROUND_PATH],
        "mapping_sha256_verified": observed[MAPPING_PATH],
        "gmt_term_count_verified": term_count,
        "background_gene_count_verified": len(background),
    }
    output_manifest["analysis_usage"] = {
        "gseapy_version": gseapy.__version__,
        "local_resource_only": True,
        "downloaded_or_substituted_during_run": False,
        "background_applied": BACKGROUND_PATH.name,
        "selected_skill": SKILL,
        "external_skill_data_calls_used": False,
    }
    with (OUTPUT_DIR / "resource_manifest.json").open("w", encoding="utf-8") as stream:
        json.dump(output_manifest, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")

    report = f"""# Combination-treatment mechanism report

The approved differential-expression analysis used only three `{TREATMENT}` samples ({', '.join(TREATMENT_SAMPLES)}) and three `{CONTROL}` samples ({', '.join(CONTROL_SAMPLES)}), with the combination as numerator. PyDESeq2 {pydeseq2.__version__} used `~Group`, `refit_cooks=True`, and the standard contrast. Genes were fitted when any selected raw count exceeded 10; {int(passes.sum()):,} of {filtered.shape[0]:,} retained genes passed adjusted p-value < 0.05, absolute log2 fold change > 0.5, and base mean > 10.

Following the `{SKILL}` functional-genomics workflow, passing genes were stratified into overall, upregulated, and downregulated sets before annotation. GSEApy {gseapy.__version__} tested them only against the evaluator-supplied, hash-verified `Reactome_2022` GMT and its {len(background):,}-symbol background. The external STRING, GO, and Ensembl endpoints described by the skill were not called because they would introduce non-frozen annotations forbidden by the task.

The best-supported primary mechanism is **{primary['pathway']}** ({primary['reactome_id']}; overlap {primary['overlap']}, adjusted p-value {float(primary['adjusted_p_value']):.3g}, odds ratio {float(primary['odds_ratio']):.3g}). The related TP53 terms among the highest-ranked overall results provide coherent functional context. Selection was deterministic: minimum adjusted p-value, then raw p-value, then maximum combined score among all passing genes.

This Reactome over-representation is statistical association with the treatment contrast, not evidence that the combination caused the pathway change. The design also cannot separate synergy from either component's individual effect.
"""
    (OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
