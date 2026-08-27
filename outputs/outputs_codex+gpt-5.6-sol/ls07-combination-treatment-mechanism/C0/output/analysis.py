"""Frozen differential expression and local Reactome enrichment."""

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

OUTPUT_DIR = Path(__file__).resolve().parent
REPO_DIR = Path(__file__).resolve().parents[5]
INPUT_DIR = REPO_DIR / "inputs" / "ls07-combination-treatment-mechanism"
GMT_PATH = INPUT_DIR / "Reactome_2022.gmt"
BACKGROUND_PATH = INPUT_DIR / "Reactome_2022.background.txt"
SOURCE_MANIFEST_PATH = INPUT_DIR / "Reactome_2022.manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_or_none(value: object) -> float | None:
    if pd.isna(value):
        return None
    number = float(value)
    return number if np.isfinite(number) else None


def parse_term(term: str) -> tuple[str, str | None]:
    match = re.match(r"^(.*) (R-HSA-\d+)$", term)
    if match:
        return match.group(1), match.group(2)
    return term, None


def enrich(gene_list: list[str], label: str, background: list[str]) -> pd.DataFrame:
    result = gp.enrichr(
        gene_list=gene_list,
        gene_sets=str(GMT_PATH),
        organism="human",
        outdir=None,
        background=background,
        cutoff=1.0,
        no_plot=True,
        verbose=False,
    ).results.copy()
    result.insert(0, "analysis_set", label)
    result.insert(1, "query_gene_count", len(gene_list))
    return result


def main() -> None:
    source_manifest = json.loads(SOURCE_MANIFEST_PATH.read_text(encoding="utf-8"))
    expected_hashes = {
        GMT_PATH: source_manifest["gmt_sha256"],
        BACKGROUND_PATH: source_manifest["background_sha256"],
        INPUT_DIR / "ensg_to_gene_name.tsv": source_manifest["mapping_sha256"],
    }
    observed_hashes = {path: sha256(path) for path in expected_hashes}
    for path, expected in expected_hashes.items():
        if observed_hashes[path] != expected:
            raise ValueError(f"Frozen resource hash mismatch: {path.name}")

    with GMT_PATH.open("r", encoding="utf-8") as handle:
        gmt_term_count = sum(1 for line in handle if line.rstrip("\n\r"))
    background = [
        line.strip()
        for line in BACKGROUND_PATH.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(background) != len(set(background)):
        raise ValueError("Background must contain unique gene symbols")
    if gmt_term_count != int(source_manifest["term_count"]):
        raise ValueError("GMT term count does not match manifest")
    if len(background) != int(source_manifest["background_gene_count"]):
        raise ValueError("Background count does not match manifest")

    layout = pd.read_csv(INPUT_DIR / "sample_layout.csv")
    metadata = (
        layout.loc[layout["SampleID"].isin(SAMPLES), ["SampleID", "Group"]]
        .set_index("SampleID")
        .loc[SAMPLES]
    )
    expected_groups = [CONTROL] * 3 + [TREATMENT] * 3
    if metadata.index.tolist() != SAMPLES or metadata["Group"].tolist() != expected_groups:
        raise ValueError("The six prespecified samples or groups do not match")

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
    passes = (
        de["padj"].notna()
        & de["padj"].lt(0.05)
        & de["log2FoldChange"].abs().gt(0.5)
        & de["baseMean"].gt(10)
    )

    mapping = pd.read_csv(INPUT_DIR / "ensg_to_gene_name.tsv", sep="\t", dtype=str)
    symbol_map = mapping.drop_duplicates("ENSG").set_index("ENSG")["gene_name"]
    symbols = de.index.to_series().map(symbol_map)
    background_set = set(background)

    def symbols_for(mask: pd.Series) -> list[str]:
        values = symbols.loc[mask]
        return sorted({str(value) for value in values.dropna() if str(value) in background_set})

    all_genes = symbols_for(passes)
    up_genes = symbols_for(passes & de["log2FoldChange"].gt(0))
    down_genes = symbols_for(passes & de["log2FoldChange"].lt(0))
    if not all_genes or not up_genes or not down_genes:
        raise ValueError("At least one mapped background gene is required in each enrichment set")

    enriched = pd.concat(
        [
            enrich(all_genes, "all_passing", background),
            enrich(up_genes, "upregulated", background),
            enrich(down_genes, "downregulated", background),
        ],
        ignore_index=True,
    )
    rename = {
        "Gene_set": "gene_set",
        "Term": "term",
        "Overlap": "overlap",
        "P-value": "p_value",
        "Adjusted P-value": "adjusted_p_value",
        "Odds Ratio": "odds_ratio",
        "Combined Score": "combined_score",
        "Genes": "overlap_genes",
    }
    enriched = enriched.rename(columns=rename)
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
        raise ValueError("GSEApy result table is missing required columns")
    enriched["overlap_genes"] = enriched["overlap_genes"].map(
        lambda value: "" if pd.isna(value) else ";".join(sorted(str(value).split(";")))
    )
    parsed = enriched["term"].map(parse_term)
    enriched.insert(4, "pathway", parsed.map(lambda item: item[0]))
    enriched.insert(5, "reactome_id", parsed.map(lambda item: item[1]))
    enriched["significant"] = enriched["adjusted_p_value"].lt(0.05)
    order = pd.Categorical(
        enriched["analysis_set"],
        categories=["all_passing", "upregulated", "downregulated"],
        ordered=True,
    )
    enriched = (
        enriched.assign(_set_order=order)
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
    enriched.loc[:, output_columns].to_csv(
        OUTPUT_DIR / "pathway_enrichment.csv", index=False, na_rep=""
    )

    all_results = enriched.loc[enriched["analysis_set"].eq("all_passing")]
    significant_all = all_results.loc[all_results["significant"]]
    candidates = significant_all if not significant_all.empty else all_results
    primary = candidates.iloc[0]

    def pathway_record(row: pd.Series) -> dict[str, object]:
        return {
            "pathway": row["pathway"],
            "reactome_id": None if pd.isna(row["reactome_id"]) else row["reactome_id"],
            "analysis_set": row["analysis_set"],
            "overlap": row["overlap"],
            "p_value": finite_or_none(row["p_value"]),
            "adjusted_p_value": finite_or_none(row["adjusted_p_value"]),
            "odds_ratio": finite_or_none(row["odds_ratio"]),
            "combined_score": finite_or_none(row["combined_score"]),
            "overlap_genes": sorted(str(row["overlap_genes"]).split(";")),
        }

    mechanism_call = {
        "comparison": f"{TREATMENT} versus {CONTROL}",
        "primary_cellular_mechanism": primary["pathway"],
        "primary_reactome_id": None
        if pd.isna(primary["reactome_id"])
        else primary["reactome_id"],
        "primary_support": pathway_record(primary),
        "selection_rule": (
            "Among all-passing DE genes, choose the Reactome term with the smallest "
            "GSEApy adjusted p-value, then p-value, then largest combined score; use "
            "the same ordering without the significance restriction if none has adjusted p < 0.05."
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
            "all_passing_mapped_background_gene_count": len(all_genes),
            "upregulated_mapped_background_gene_count": len(up_genes),
            "downregulated_mapped_background_gene_count": len(down_genes),
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
        "supporting_pathways": [
            pathway_record(row) for _, row in all_results.head(5).iterrows()
        ],
        "interpretation_limit": (
            "Pathway over-representation supports an association with the treatment "
            "contrast; it does not demonstrate that the treatment caused the pathway change."
        ),
    }
    with (OUTPUT_DIR / "mechanism_call.json").open("w", encoding="utf-8") as handle:
        json.dump(mechanism_call, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")

    resource_manifest = dict(source_manifest)
    resource_manifest["verification"] = {
        "gmt_sha256_verified": observed_hashes[GMT_PATH],
        "background_sha256_verified": observed_hashes[BACKGROUND_PATH],
        "mapping_sha256_verified": observed_hashes[INPUT_DIR / "ensg_to_gene_name.tsv"],
        "gmt_term_count_verified": gmt_term_count,
        "background_gene_count_verified": len(background),
    }
    resource_manifest["analysis_usage"] = {
        "gseapy_version": gseapy.__version__,
        "local_resource_only": True,
        "downloaded_or_substituted_during_run": False,
        "background_applied": BACKGROUND_PATH.name,
    }
    with (OUTPUT_DIR / "resource_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(resource_manifest, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")

    report = f"""# Combination-treatment mechanism report

The approved six-sample differential-expression analysis compared `{TREATMENT}` (numerator; {', '.join(TREATMENT_SAMPLES)}) with `{CONTROL}` (denominator; {', '.join(CONTROL_SAMPLES)}), using only `Group` in the design. PyDESeq2 {pydeseq2.__version__} used `refit_cooks=True`; genes were retained when any selected raw count exceeded 10. Of {filtered.shape[0]:,} retained genes, {int(passes.sum()):,} passed adjusted p-value < 0.05, absolute log2 fold change > 0.5, and base mean > 10.

GSEApy {gseapy.__version__} tested the mapped passing genes against the evaluator-supplied, byte-verified `Reactome_2022` GMT and its explicit {len(background):,}-gene background. No current pathway library or alternate mapping was downloaded. The all-passing query contained {len(all_genes):,} unique background-mapped symbols; the up- and downregulated queries contained {len(up_genes):,} and {len(down_genes):,}, respectively.

The best-supported primary cellular mechanism is **{primary['pathway']}** ({primary['reactome_id']}; overlap {primary['overlap']}, adjusted p-value {float(primary['adjusted_p_value']):.3g}, odds ratio {float(primary['odds_ratio']):.3g}). This call follows a prespecified deterministic rule: among all-passing genes, minimize adjusted p-value and then raw p-value, with combined score as the next discriminator. The complete overall and directional enrichment tables are in `pathway_enrichment.csv`; resource provenance and integrity checks are in `resource_manifest.json`.

This is pathway over-representation associated with the treatment contrast, not proof that the combination treatment caused the pathway change. The experiment also does not isolate pharmacologic interaction from the effects of either component alone.
"""
    (OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
