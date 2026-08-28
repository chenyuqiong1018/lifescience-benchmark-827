#!/usr/bin/env python3
"""Audited, dependency-free RNA/ATAC population matching workflow."""

from __future__ import annotations

import csv
import gzip
import itertools
import math
from collections import Counter
from pathlib import Path


REPO = Path(__file__).resolve().parents[5]
INPUTS = REPO / "inputs" / "ls08-multiome-column-match"
OUTPUT = Path(__file__).resolve().parent
RNA_PATH = INPUTS / "multiome.match.atac.rna.q1.rna.tsv.gz"
ATAC_PATH = INPUTS / "multiome.match.atac.rna.q1.atac.tsv.gz"
ANNOTATION_PATH = INPUTS / "ensembl112_gene_coordinates.tsv"
VALID_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}


def parse_measurement(text: str, context: str) -> float:
    if not text.strip():
        raise ValueError(f"Missing measurement at {context}; no imputation is permitted")
    number = float(text)
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"Non-finite or negative measurement at {context}")
    return number


def load_rna() -> tuple[list[str], dict[str, list[float]], int]:
    values: dict[str, list[float]] = {}
    row_count = 0
    with gzip.open(RNA_PATH, "rt", encoding="utf-8", newline="") as stream:
        rows = csv.reader(stream, delimiter="\t")
        header = next(rows)
        if header[0] != "gene" or len(header) != 9 or len(set(header[1:])) != 8:
            raise ValueError("RNA header must contain gene and eight unique populations")
        for line, row in enumerate(rows, 2):
            row_count += 1
            if len(row) != 9:
                raise ValueError(f"RNA line {line} has {len(row)} fields, expected 9")
            gene = row[0]
            if not gene or gene in values:
                raise ValueError(f"Blank or duplicate RNA gene at line {line}: {gene!r}")
            values[gene] = [parse_measurement(x, f"RNA line {line}") for x in row[1:]]
    return header[1:], values, row_count


def load_annotation(rna_genes: set[str]) -> tuple[dict[str, str], dict[str, int]]:
    annotation_rows: list[dict[str, str]] = []
    symbol_counts: Counter[str] = Counter()
    with ANNOTATION_PATH.open("r", encoding="utf-8", newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        required = {
            "Gene name", "Chromosome/scaffold name", "Gene start (bp)",
            "Gene end (bp)", "Strand",
        }
        if rows.fieldnames is None or not required.issubset(rows.fieldnames):
            raise ValueError("Unexpected Ensembl annotation schema")
        for row in rows:
            symbol = row["Gene name"]
            if symbol:
                symbol_counts[symbol] += 1
            annotation_rows.append(row)

    gene_to_bin: dict[str, str] = {}
    eligible = 0
    for row in annotation_rows:
        gene = row["Gene name"]
        chrom = row["Chromosome/scaffold name"]
        if gene not in rna_genes or symbol_counts[gene] != 1 or chrom not in VALID_CHROMS:
            continue
        eligible += 1
        strand = int(row["Strand"])
        if strand == 1:
            tss_one_based = int(row["Gene start (bp)"])
        elif strand == -1:
            tss_one_based = int(row["Gene end (bp)"])
        else:
            raise ValueError(f"Unsupported strand {strand} for {gene}")
        # Ensembl coordinates are one-based; ATAC bins are zero-based, half-open.
        bin_start = ((tss_one_based - 1) // 10_000) * 10_000
        gene_to_bin[gene] = f"chr{chrom}_{bin_start}_{bin_start + 10_000}"
    audit = {
        "annotation_rows": len(annotation_rows),
        "nonblank_unique_symbols": sum(count == 1 for count in symbol_counts.values()),
        "eligible_genes": eligible,
    }
    return gene_to_bin, audit


def load_atac(needed_bins: set[str]) -> tuple[list[str], dict[str, list[float]], int]:
    retained: dict[str, list[float]] = {}
    seen_bins: set[str] = set()
    row_count = 0
    with gzip.open(ATAC_PATH, "rt", encoding="utf-8", newline="") as stream:
        rows = csv.reader(stream, delimiter="\t")
        header = next(rows)
        if header[0] != "peak" or len(header) != 9 or len(set(header[1:])) != 8:
            raise ValueError("ATAC header must contain peak and eight unique columns")
        for line, row in enumerate(rows, 2):
            row_count += 1
            if len(row) != 9:
                raise ValueError(f"ATAC line {line} has {len(row)} fields, expected 9")
            peak = row[0]
            if not peak or peak in seen_bins:
                raise ValueError(f"Blank or duplicate ATAC bin at line {line}: {peak!r}")
            seen_bins.add(peak)
            parsed = [parse_measurement(x, f"ATAC line {line}") for x in row[1:]]
            if peak in needed_bins:
                retained[peak] = parsed
    return header[1:], retained, row_count


def variance(vector: list[float]) -> float:
    center = sum(vector) / len(vector)
    return sum((value - center) ** 2 for value in vector) / len(vector)


def correlation(x: list[float], y: list[float]) -> float:
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("Pearson correlation needs aligned vectors of length at least two")
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)
    centered_x = [value - mean_x for value in x]
    centered_y = [value - mean_y for value in y]
    scale_x = sum(value * value for value in centered_x)
    scale_y = sum(value * value for value in centered_y)
    if scale_x <= 0 or scale_y <= 0:
        raise ValueError("Pearson correlation undefined for constant profiles")
    return sum(a * b for a, b in zip(centered_x, centered_y)) / math.sqrt(scale_x * scale_y)


def hungarian_maximum(weights: list[list[float]]) -> list[int]:
    """Solve a square maximum-weight assignment in O(n^3)."""
    n = len(weights)
    if n == 0 or any(len(row) != n for row in weights):
        raise ValueError("Assignment matrix must be non-empty and square")
    row_potential = [0.0] * (n + 1)
    col_potential = [0.0] * (n + 1)
    matched_row = [0] * (n + 1)
    predecessor = [0] * (n + 1)
    for new_row in range(1, n + 1):
        matched_row[0] = new_row
        unused_cost = [math.inf] * (n + 1)
        used = [False] * (n + 1)
        current_col = 0
        while True:
            used[current_col] = True
            current_row = matched_row[current_col]
            delta = math.inf
            next_col = 0
            for col in range(1, n + 1):
                if used[col]:
                    continue
                reduced = -weights[current_row - 1][col - 1] - row_potential[current_row] - col_potential[col]
                if reduced < unused_cost[col]:
                    unused_cost[col] = reduced
                    predecessor[col] = current_col
                if unused_cost[col] < delta:
                    delta = unused_cost[col]
                    next_col = col
            for col in range(n + 1):
                if used[col]:
                    row_potential[matched_row[col]] += delta
                    col_potential[col] -= delta
                else:
                    unused_cost[col] -= delta
            current_col = next_col
            if matched_row[current_col] == 0:
                break
        while current_col:
            previous = predecessor[current_col]
            matched_row[current_col] = matched_row[previous]
            current_col = previous
    assignment = [-1] * n
    for col in range(1, n + 1):
        assignment[matched_row[col] - 1] = col - 1
    if sorted(assignment) != list(range(n)):
        raise AssertionError("Hungarian result is not a bijection")
    return assignment


def exhaustive_audit(weights: list[list[float]], assignment: list[int]) -> tuple[float, float]:
    totals = sorted(
        (sum(weights[row][col] for row, col in enumerate(permutation)) for permutation in itertools.permutations(range(8))),
        reverse=True,
    )
    assigned_total = sum(weights[row][col] for row, col in enumerate(assignment))
    if not math.isclose(assigned_total, totals[0], rel_tol=0.0, abs_tol=1e-12):
        raise AssertionError("Hungarian result disagrees with exhaustive 8! audit")
    return totals[0], totals[1]


def main() -> None:
    rna_names, rna, rna_rows = load_rna()
    gene_to_bin, audit = load_annotation(set(rna))
    atac_names, atac, atac_rows = load_atac(set(gene_to_bin.values()))
    mapped = [(gene, gene_to_bin[gene]) for gene in gene_to_bin if gene_to_bin[gene] in atac]
    audit.update({"rna_rows": rna_rows, "atac_rows": atac_rows, "mapped_genes": len(mapped)})
    if len(mapped) < 2_000:
        raise ValueError(f"Only {len(mapped)} genes map to observed ATAC bins")

    ranked = sorted(
        mapped,
        key=lambda item: (-variance([math.log1p(value) for value in rna[item[0]]]), item[0]),
    )
    selected = ranked[:2_000]
    rna_profiles = [
        [math.log1p(rna[gene][population]) for gene, _ in selected]
        for population in range(8)
    ]
    atac_profiles = [
        [math.log1p(atac[bin_name][column]) for _, bin_name in selected]
        for column in range(8)
    ]
    matrix = [[correlation(rna_profile, atac_profile) for atac_profile in atac_profiles] for rna_profile in rna_profiles]
    assignment = hungarian_maximum(matrix)
    best_total, next_total = exhaustive_audit(matrix, assignment)

    with (OUTPUT / "score_matrix.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["rna_population", *atac_names])
        for name, row in zip(rna_names, matrix):
            writer.writerow([name, *[f"{score:.12f}" for score in row]])

    row_runner_ups: list[float] = []
    with (OUTPUT / "column_mapping.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["rna_population", "atac_column", "match_score", "runner_up_score"])
        for row, name in enumerate(rna_names):
            runner_up = sorted(matrix[row], reverse=True)[1]
            row_runner_ups.append(runner_up)
            assigned_col = assignment[row]
            writer.writerow([
                name,
                atac_names[assigned_col],
                f"{matrix[row][assigned_col]:.12f}",
                f"{runner_up:.12f}",
            ])

    mapping = ", ".join(f"RNA {rna_names[i]} → ATAC {atac_names[assignment[i]]}" for i in range(8))
    row_margins = [matrix[i][assignment[i]] - row_runner_ups[i] for i in range(8)]
    report = f"""# Multiome column matching

## Result

The recovered bijection is {mapping}. The maximum total Pearson correlation is {best_total:.6f}. The next-best complete assignment totals {next_total:.6f}, a global gap of {best_total - next_total:.6f}. The smallest assigned-minus-row-runner-up margin is {min(row_margins):.6f}; a negative row margin can occur because the optimum is global rather than eight independent row choices.

## Shared biological signal and method

The shared signal is coordinated gene activity across molecular layers: accessible chromatin at a gene's transcription start site tends to accompany expression of that gene. Following the supplied rule, Ensembl release 112 symbols were retained only when they occurred exactly once in the complete annotation, appeared in RNA, and lay on chromosomes 1–22, X, or Y. Strand determined whether gene start or end was the TSS. One-based coordinates were converted to zero-based half-open 10 kb ATAC bins. RNA TPM and mapped ATAC values were log1p-transformed; the 2,000 mapped genes with greatest across-RNA-population variance were used for all 64 Pearson correlations. An O(n³) Hungarian maximum-weight solution enforced the bijection, and exhaustive enumeration of all 8! assignments independently confirmed its total.

## Data-quality audit

The inputs contained {audit['rna_rows']:,} unique RNA genes, {audit['atac_rows']:,} unique ATAC bins, and {audit['annotation_rows']:,} annotation rows. After exact-symbol/chromosome filtering there were {audit['eligible_genes']:,} eligible genes; {audit['mapped_genes']:,} had an observed TSS bin. Every matrix row had nine fields, all eight measurement columns were unique, and no blank, negative, or non-finite measurement was found. No zero imputation was performed. `runner_up_score` is the second-highest row-wise correlation, not the second-best global assignment. Correlations describe association, not gene-level causation.

The five prescribed skills guided cross-omics framing, regulatory-coordinate handling, integrity checks, correlation interpretation, and executable-code verification; external reference datasets were neither needed nor used.
"""
    (OUTPUT / "report.md").write_text(report, encoding="utf-8", newline="\n")
    print(f"Audit: {audit}")
    print(f"Best total={best_total:.6f}; global gap={best_total - next_total:.6f}")
    print(mapping)


if __name__ == "__main__":
    main()
