#!/usr/bin/env python3
"""Recover the RNA-to-ATAC population matching from the supplied matrices."""

from __future__ import annotations

import csv
import gzip
import math
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[5]
INPUT = ROOT / "inputs" / "ls08-multiome-column-match"
OUT = Path(__file__).resolve().parent
RNA_FILE = INPUT / "multiome.match.atac.rna.q1.rna.tsv.gz"
ATAC_FILE = INPUT / "multiome.match.atac.rna.q1.atac.tsv.gz"
COORD_FILE = INPUT / "ensembl112_gene_coordinates.tsv"
AUTOSOMES_SEX = {str(i) for i in range(1, 23)} | {"X", "Y"}


def finite_nonnegative(text: str, label: str) -> float:
    if text == "":
        raise ValueError(f"Missing value in {label}; zero imputation is not allowed")
    value = float(text)
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"Invalid value in {label}: {text!r}")
    return value


def read_rna() -> tuple[list[str], dict[str, list[float]]]:
    with gzip.open(RNA_FILE, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[0] != "gene" or len(header) != 9:
            raise ValueError("Expected gene plus eight RNA population columns")
        populations = header[1:]
        expression: dict[str, list[float]] = {}
        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"RNA width mismatch at line {row_number}")
            gene = row[0]
            if gene in expression:
                raise ValueError(f"Duplicate RNA gene: {gene}")
            expression[gene] = [finite_nonnegative(x, f"RNA line {row_number}") for x in row[1:]]
    return populations, expression


def unique_gene_tss(rna_genes: set[str]) -> dict[str, tuple[str, int]]:
    rows: list[dict[str, str]] = []
    counts: Counter[str] = Counter()
    with COORD_FILE.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "Gene name", "Chromosome/scaffold name", "Gene start (bp)",
            "Gene end (bp)", "Strand",
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Coordinate annotation has unexpected columns")
        for row in reader:
            symbol = row["Gene name"]
            if symbol:
                counts[symbol] += 1
            rows.append(row)

    result: dict[str, tuple[str, int]] = {}
    for row in rows:
        symbol = row["Gene name"]
        chrom = row["Chromosome/scaffold name"]
        if symbol not in rna_genes or counts[symbol] != 1 or chrom not in AUTOSOMES_SEX:
            continue
        strand = int(row["Strand"])
        if strand not in (-1, 1):
            raise ValueError(f"Unexpected strand for {symbol}: {strand}")
        start = int(row["Gene start (bp)"])
        end = int(row["Gene end (bp)"])
        tss_one_based = start if strand == 1 else end
        tss_zero_based = tss_one_based - 1
        bin_start = (tss_zero_based // 10_000) * 10_000
        result[symbol] = (f"chr{chrom}_{bin_start}_{bin_start + 10_000}", tss_one_based)
    return result


def read_needed_atac(needed: set[str]) -> tuple[list[str], dict[str, list[float]]]:
    retained: dict[str, list[float]] = {}
    with gzip.open(ATAC_FILE, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[0] != "peak" or len(header) != 9:
            raise ValueError("Expected peak plus eight ATAC columns")
        atac_columns = header[1:]
        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"ATAC width mismatch at line {row_number}")
            if row[0] in needed:
                if row[0] in retained:
                    raise ValueError(f"Duplicate ATAC bin: {row[0]}")
                retained[row[0]] = [finite_nonnegative(x, f"ATAC line {row_number}") for x in row[1:]]
    return atac_columns, retained


def population_variance(values: list[float]) -> float:
    mean = sum(values) / len(values)
    return sum((x - mean) ** 2 for x in values) / len(values)


def pearson(x: list[float], y: list[float]) -> float:
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("Pearson vectors must have equal non-trivial length")
    mx, my = sum(x) / len(x), sum(y) / len(y)
    numerator = sum((a - mx) * (b - my) for a, b in zip(x, y))
    dx = sum((a - mx) ** 2 for a in x)
    dy = sum((b - my) ** 2 for b in y)
    if dx == 0 or dy == 0:
        raise ValueError("Pearson correlation is undefined for a constant vector")
    return numerator / math.sqrt(dx * dy)


def hungarian_max(weights: list[list[float]]) -> list[int]:
    """Maximum-weight square assignment; return the assigned column per row."""
    n = len(weights)
    if n == 0 or any(len(row) != n for row in weights):
        raise ValueError("Hungarian algorithm requires a non-empty square matrix")
    # Standard shortest-augmenting-path Hungarian algorithm on costs=-weights.
    u = [0.0] * (n + 1)
    v = [0.0] * (n + 1)
    p = [0] * (n + 1)
    way = [0] * (n + 1)
    for i in range(1, n + 1):
        p[0] = i
        minv = [math.inf] * (n + 1)
        used = [False] * (n + 1)
        j0 = 0
        while True:
            used[j0] = True
            i0 = p[j0]
            delta = math.inf
            j1 = 0
            for j in range(1, n + 1):
                if used[j]:
                    continue
                cur = -weights[i0 - 1][j - 1] - u[i0] - v[j]
                if cur < minv[j]:
                    minv[j] = cur
                    way[j] = j0
                if minv[j] < delta:
                    delta = minv[j]
                    j1 = j
            for j in range(n + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break
    assignment = [-1] * n
    for j in range(1, n + 1):
        assignment[p[j] - 1] = j - 1
    if sorted(assignment) != list(range(n)):
        raise AssertionError("Assignment is not bijective")
    return assignment


def main() -> None:
    rna_populations, rna = read_rna()
    tss = unique_gene_tss(set(rna))
    atac_columns, atac = read_needed_atac({entry[0] for entry in tss.values()})

    mapped = [(gene, bin_name) for gene, (bin_name, _) in tss.items() if bin_name in atac]
    if len(mapped) < 2_000:
        raise ValueError(f"Only {len(mapped)} mapped genes; 2,000 are required")
    ranked = sorted(
        mapped,
        key=lambda pair: (-population_variance([math.log1p(x) for x in rna[pair[0]]]), pair[0]),
    )
    selected = ranked[:2_000]

    rna_vectors = [
        [math.log1p(rna[gene][column]) for gene, _ in selected]
        for column in range(8)
    ]
    atac_vectors = [
        [math.log1p(atac[bin_name][column]) for _, bin_name in selected]
        for column in range(8)
    ]
    scores = [[pearson(rv, av) for av in atac_vectors] for rv in rna_vectors]
    assignment = hungarian_max(scores)

    with (OUT / "score_matrix.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["rna_population", *atac_columns])
        for row_name, row in zip(rna_populations, scores):
            writer.writerow([row_name, *[f"{value:.12f}" for value in row]])

    with (OUT / "column_mapping.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["rna_population", "atac_column", "match_score", "runner_up_score"])
        for row_index, row_name in enumerate(rna_populations):
            assigned = assignment[row_index]
            second = sorted(scores[row_index], reverse=True)[1]
            writer.writerow([
                row_name,
                atac_columns[assigned],
                f"{scores[row_index][assigned]:.12f}",
                f"{second:.12f}",
            ])

    total = sum(scores[i][assignment[i]] for i in range(8))
    margins = [scores[i][assignment[i]] - sorted(scores[i], reverse=True)[1] for i in range(8)]
    mapping_text = ", ".join(
        f"RNA {rna_populations[i]} → ATAC {atac_columns[assignment[i]]}" for i in range(8)
    )
    report = f"""# Multiome column matching

## Result

The maximum-total-correlation bijection is: {mapping_text}. The assigned correlations sum to {total:.6f}; the smallest assigned-minus-row-runner-up margin is {min(margins):.6f}.

## Shared biological signal

RNA abundance and chromatin accessibility are different molecular layers, but active regulatory programs tend to make transcription start sites accessible while their linked genes are expressed. I therefore compared populations through matched gene-level profiles: each eligible gene's strand-aware TSS was mapped to its containing 10 kb ATAC bin, and the corresponding RNA TPM and ATAC value were log-transformed. This is the cross-modal biological signal used for matching; it does not imply a causal relationship for individual genes.

## Procedure and safeguards

The Ensembl release 112 annotation was filtered to symbols occurring exactly once in the complete annotation, present in RNA, and located on chromosomes 1–22, X, or Y. Coordinates were treated as one-based and converted to zero-based half-open ATAC bins. No missing value was replaced with zero. Among genes with an observed target bin, the 2,000 with highest variance across log1p RNA populations were selected. All 64 Pearson correlations were computed across the same ordered genes, then an O(n³) Hungarian maximum-weight assignment enforced the required bijection. `runner_up_score` is the second-highest correlation in that RNA row, irrespective of the global assignment.

The freshly installed `multiomics_integration` skill informed the cross-layer integration framing; its external TCGA, UniProt, STRING, and KEGG calls were not relevant to this self-contained column-matching input and were not used.
"""
    (OUT / "report.md").write_text(report, encoding="utf-8", newline="\n")
    print(f"Selected genes: {len(selected)}; total assignment score: {total:.6f}")
    print(mapping_text)


if __name__ == "__main__":
    main()
