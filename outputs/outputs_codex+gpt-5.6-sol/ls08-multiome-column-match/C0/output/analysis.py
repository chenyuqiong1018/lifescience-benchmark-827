"""Recover the bijection between RNA populations and permuted ATAC columns."""

from __future__ import annotations

import csv
import gzip
import itertools
import math
from collections import Counter
from pathlib import Path


OUTPUT_DIR = Path(__file__).resolve().parent
REPO_DIR = Path(__file__).resolve().parents[5]
INPUT_DIR = REPO_DIR / "inputs" / "ls08-multiome-column-match"
RNA_PATH = INPUT_DIR / "multiome.match.atac.rna.q1.rna.tsv.gz"
ATAC_PATH = INPUT_DIR / "multiome.match.atac.rna.q1.atac.tsv.gz"
COORD_PATH = INPUT_DIR / "ensembl112_gene_coordinates.tsv"
AUTOSOMES_AND_SEX = {str(number) for number in range(1, 23)} | {"X", "Y"}
TOP_GENE_COUNT = 2000


def parse_numeric_row(row: list[str], label: str) -> list[float]:
    if any(value == "" for value in row):
        raise ValueError(f"Missing numeric value in {label}; no zero imputation is allowed")
    values = [float(value) for value in row]
    if any(not math.isfinite(value) or value < 0 for value in values):
        raise ValueError(f"Non-finite or negative value in {label}")
    return values


def variance(values: list[float]) -> float:
    center = sum(values) / len(values)
    return sum((value - center) ** 2 for value in values) / len(values)


def pearson(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or not left:
        raise ValueError("Pearson vectors must have equal positive length")
    left_mean = sum(left) / len(left)
    right_mean = sum(right) / len(right)
    numerator = sum(
        (x - left_mean) * (y - right_mean) for x, y in zip(left, right)
    )
    left_ss = sum((x - left_mean) ** 2 for x in left)
    right_ss = sum((y - right_mean) ** 2 for y in right)
    if left_ss <= 0 or right_ss <= 0:
        raise ValueError("Pearson correlation is undefined for a constant vector")
    return numerator / math.sqrt(left_ss * right_ss)


def main() -> None:
    with gzip.open(RNA_PATH, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[0] != "gene" or len(header) != 9:
            raise ValueError("Unexpected RNA matrix schema")
        rna_populations = header[1:]
        rna: dict[str, list[float]] = {}
        for row in reader:
            gene = row[0]
            if gene in rna:
                raise ValueError(f"Duplicate RNA gene symbol: {gene}")
            rna[gene] = parse_numeric_row(row[1:], f"RNA gene {gene}")

    with COORD_PATH.open("r", encoding="utf-8", newline="") as handle:
        coordinate_rows = list(csv.DictReader(handle, delimiter="\t"))
    symbol_counts = Counter(row["Gene name"] for row in coordinate_rows if row["Gene name"])
    gene_to_bin: dict[str, str] = {}
    for row in coordinate_rows:
        symbol = row["Gene name"]
        chromosome = row["Chromosome/scaffold name"]
        if (
            not symbol
            or symbol_counts[symbol] != 1
            or symbol not in rna
            or chromosome not in AUTOSOMES_AND_SEX
        ):
            continue
        strand = int(row["Strand"])
        if strand not in (-1, 1):
            raise ValueError(f"Unexpected strand for {symbol}: {strand}")
        tss_1_based = (
            int(row["Gene start (bp)"]) if strand == 1 else int(row["Gene end (bp)"])
        )
        if tss_1_based < 1:
            raise ValueError(f"Invalid TSS for {symbol}")
        tss_zero_based = tss_1_based - 1
        bin_start = (tss_zero_based // 10_000) * 10_000
        gene_to_bin[symbol] = f"chr{chromosome}_{bin_start}_{bin_start + 10_000}"

    needed_bins = set(gene_to_bin.values())
    with gzip.open(ATAC_PATH, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[0] != "peak" or len(header) != 9:
            raise ValueError("Unexpected ATAC matrix schema")
        atac_columns = header[1:]
        atac: dict[str, list[float]] = {}
        for row in reader:
            peak = row[0]
            if peak in needed_bins:
                if peak in atac:
                    raise ValueError(f"Duplicate ATAC bin: {peak}")
                atac[peak] = parse_numeric_row(row[1:], f"ATAC bin {peak}")

    mapped_genes = sorted(symbol for symbol, peak in gene_to_bin.items() if peak in atac)
    if len(mapped_genes) < TOP_GENE_COUNT:
        raise ValueError(f"Only {len(mapped_genes)} genes map to supplied ATAC bins")
    logged_rna = {
        gene: [math.log1p(value) for value in rna[gene]] for gene in mapped_genes
    }
    selected_genes = sorted(
        mapped_genes,
        key=lambda gene: (-variance(logged_rna[gene]), gene),
    )[:TOP_GENE_COUNT]

    rna_vectors = [
        [logged_rna[gene][index] for gene in selected_genes]
        for index in range(len(rna_populations))
    ]
    atac_vectors = [
        [math.log1p(atac[gene_to_bin[gene]][index]) for gene in selected_genes]
        for index in range(len(atac_columns))
    ]
    scores = [
        [pearson(rna_vector, atac_vector) for atac_vector in atac_vectors]
        for rna_vector in rna_vectors
    ]

    best_permutation: tuple[int, ...] | None = None
    best_total = -math.inf
    for permutation in itertools.permutations(range(len(atac_columns))):
        total = sum(scores[row][column] for row, column in enumerate(permutation))
        if total > best_total:
            best_total = total
            best_permutation = permutation
    if best_permutation is None:
        raise RuntimeError("No bijective assignment found")

    with (OUTPUT_DIR / "column_mapping.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        columns = ["rna_population", "atac_column", "match_score", "runner_up_score"]
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row_index, rna_population in enumerate(rna_populations):
            assigned = best_permutation[row_index]
            runner_up = sorted(scores[row_index], reverse=True)[1]
            writer.writerow(
                {
                    "rna_population": rna_population,
                    "atac_column": atac_columns[assigned],
                    "match_score": format(scores[row_index][assigned], ".12g"),
                    "runner_up_score": format(runner_up, ".12g"),
                }
            )

    with (OUTPUT_DIR / "score_matrix.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["rna_population", *atac_columns])
        for rna_population, row in zip(rna_populations, scores):
            writer.writerow([rna_population, *[format(value, ".12g") for value in row]])

    mapping_text = ", ".join(
        f"RNA {rna_populations[index]}→ATAC {atac_columns[column]}"
        for index, column in enumerate(best_permutation)
    )
    minimum_margin = min(
        scores[row][best_permutation[row]] - sorted(scores[row], reverse=True)[1]
        for row in range(len(rna_populations))
    )
    report = f"""# Multiome column matching

The frozen rule retained Ensembl 112 gene symbols occurring exactly once, present in RNA, and located on chromosomes 1-22, X, or Y. Strand-aware 1-based transcription start sites were converted to 0-based coordinates and assigned to their containing 10 kb ATAC bins. Both RNA TPM and mapped ATAC signal were transformed with `log1p`; the {TOP_GENE_COUNT:,} mapped genes with greatest variance across RNA populations were retained.

The shared biological signal is gene-proximal chromatin accessibility and expression covariation across the same populations: accessibility at a gene's TSS-containing bin should track that gene's RNA abundance. Pearson correlations across the selected genes produced the 8×8 score matrix. Exhaustive evaluation of all 8! assignments, equivalent to the maximum-weight Hungarian solution at this size, enforced a bijection and maximized total correlation ({best_total:.6f}).

Recovered mapping: {mapping_text}. `match_score` is the correlation selected by the global bijection; `runner_up_score` is the second-highest correlation within that RNA row. The smallest assigned-versus-row-runner-up margin is {minimum_margin:.6f}. Scores measure cross-modal concordance, not causal regulation, and the assignment is specific to the supplied normalization and gene-coordinate release.
"""
    (OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
