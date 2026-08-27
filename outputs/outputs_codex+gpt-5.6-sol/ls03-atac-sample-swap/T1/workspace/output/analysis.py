#!/usr/bin/env python3
"""Detect a unique ATAC-seq organ-label swap from promoter accessibility.

The analysis is intentionally self-contained and deterministic. It maps gene TSSs
from the whole chromosome-arm GTF coordinates to the split ATAC table contigs,
normalizes every sample by its full-library count, and asks which reciprocal label
exchange most improves tissue-marker promoter coherence.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path


ORGANS = [
    "Bladder", "Brain", "Cloaca", "GallBladder", "Gill", "Heart",
    "Intestine", "Kidney", "Limb", "Liver", "Lung", "Pancreas",
    "Prostate", "Spleen", "Stomach",
]

MARKERS = {
    "Bladder": ["UPK1A", "UPK1B", "UPK2", "UPK3A", "KRT20", "KRT13"],
    "Brain": ["SOX2", "GFAP", "MAP2", "RBFOX3", "NEUROD1", "SLC17A7", "SNAP25"],
    "Cloaca": ["HOXA13", "HOXD13", "KRT13", "KRT14"],
    "GallBladder": ["KRT19", "MUC1", "EPCAM", "SOX17", "KRT8"],
    "Gill": ["FOXI1", "KRT4", "ATP6V1B1", "CA2", "GATA3"],
    "Heart": ["MYH6", "MYH7", "TNNT2", "NKX2-5", "ACTC1", "PLN"],
    "Intestine": ["CDX2", "VIL1", "LGR5", "SI", "MUC2", "FABP2"],
    "Kidney": ["PAX2", "PAX8", "SLC12A1", "UMOD", "NPHS1", "AQP2"],
    "Limb": ["HOXA13", "HOXD13", "TBX5", "TBX4", "FGF8", "PRRX1"],
    "Liver": ["ALB", "AFP", "APOA1", "TTR", "HNF4A", "CYP3A4", "FGA"],
    "Lung": ["SFTPA1", "SFTPB", "SFTPC", "NKX2-1", "SCGB1A1"],
    "Pancreas": ["INS", "GCG", "PDX1", "PRSS1", "AMY2A", "CPA1", "CEL"],
    "Prostate": ["AR", "NKX3-1", "KLK3", "HOXB13"],
    "Spleen": ["SPI1", "PTPRC", "CD3D", "CD79A", "MS4A1"],
    "Stomach": ["GAST", "GIF", "ATP4A", "ATP4B", "MUC5AC", "PGA3"],
}

BIN_SIZE = 10_000
PROMOTER_RADIUS_BINS = 1


def gene_symbols(attribute_text: str) -> set[str]:
    """Return normalized symbols from a GTF gene_name attribute."""
    match = re.search(r'gene_name\s+"([^"]+)"', attribute_text)
    if not match:
        return set()
    symbols = set()
    for token in match.group(1).split("|"):
        symbol = token.strip().split()[0].split("[")[0].upper()
        if symbol:
            symbols.add(symbol)
    return symbols


def read_contigs(path: Path):
    sizes = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            chrom, size = line.rstrip("\n").split("\t")
            sizes.append((chrom, int(size)))

    by_arm = defaultdict(list)
    for chrom, size in sizes:
        arm = re.sub(r"_\d+$", "", chrom)
        by_arm[arm].append((chrom, size))
    return sizes, by_arm


def split_coordinate(arm: str, pos0: int, by_arm):
    """Map a zero-based whole-arm coordinate to a split-contig coordinate."""
    offset = 0
    for chrom, size in by_arm.get(arm, []):
        if pos0 < offset + size:
            return chrom, pos0 - offset
        offset += size
    return None


def collect_marker_tss(gtf_path: Path, by_arm):
    wanted = {symbol for values in MARKERS.values() for symbol in values}
    loci = defaultdict(list)
    with gzip.open(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            symbols = gene_symbols(fields[8]) & wanted
            if not symbols:
                continue
            start, end = int(fields[3]), int(fields[4])
            tss0 = start - 1 if fields[6] == "+" else end - 1
            mapped = split_coordinate(fields[0], tss0, by_arm)
            if mapped is None:
                continue
            chrom, local0 = mapped
            for symbol in symbols:
                loci[symbol].append((chrom, local0))
    return loci


def desired_promoter_bins(loci):
    desired = set()
    marker_bin_sets = defaultdict(list)
    for tissue, symbols in MARKERS.items():
        for symbol in symbols:
            for chrom, pos0 in loci.get(symbol, []):
                center = (pos0 // BIN_SIZE) * BIN_SIZE
                bins = []
                for delta in range(-PROMOTER_RADIUS_BINS, PROMOTER_RADIUS_BINS + 1):
                    start = center + delta * BIN_SIZE
                    if start >= 0:
                        key = (chrom, start)
                        desired.add(key)
                        bins.append(key)
                if bins:
                    marker_bin_sets[tissue].append((symbol, bins))
    return desired, marker_bin_sets


def read_atac(path: Path, desired):
    totals = [0.0] * len(ORGANS)
    captured = {}
    rows = 0
    widths = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[:3] != ["chrom", "start", "end"] or header[3:] != ORGANS:
            raise ValueError(f"Unexpected ATAC header: {header}")
        for fields in reader:
            rows += 1
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            widths.add(end - start)
            values = [float(value) for value in fields[3:]]
            if len(values) != len(ORGANS) or not all(math.isfinite(value) for value in values):
                raise ValueError(f"Invalid numeric row {rows + 1}")
            for i, value in enumerate(values):
                totals[i] += value
            if (chrom, start) in desired:
                captured[(chrom, start)] = values
    if not all(total > 0 for total in totals):
        raise ValueError("All sample library totals must be positive")
    return totals, captured, rows, sorted(widths)


def standardized_promoter_vectors(marker_bin_sets, captured, totals):
    """Create dimensionless promoter profiles with across-organ mean equal to 1."""
    vectors = defaultdict(list)
    matched_genes = defaultdict(set)
    for tissue, loci in marker_bin_sets.items():
        for symbol, bins in loci:
            raw = [0.0] * len(ORGANS)
            observed = 0
            for key in bins:
                values = captured.get(key)
                if values is None:
                    continue
                observed += 1
                for i, value in enumerate(values):
                    raw[i] += value
            if observed == 0:
                continue
            cpm = [raw[i] * 1_000_000.0 / totals[i] for i in range(len(ORGANS))]
            mean = statistics.fmean(cpm)
            if mean <= 0:
                continue
            vectors[tissue].append([value / mean for value in cpm])
            matched_genes[tissue].add(symbol)
    missing_tissues = [tissue for tissue in ORGANS if not vectors[tissue]]
    if missing_tissues:
        raise ValueError(f"No informative marker promoters for: {missing_tissues}")
    return vectors, matched_genes


def mean_profiles(vectors):
    return {
        tissue: [statistics.fmean(row[i] for row in rows) for i in range(len(ORGANS))]
        for tissue, rows in vectors.items()
    }


def rank_swaps(profiles):
    records = []
    for i, organ_a in enumerate(ORGANS):
        for organ_b in ORGANS[i + 1:]:
            ia, ib = ORGANS.index(organ_a), ORGANS.index(organ_b)
            gain = (
                profiles[organ_a][ib] + profiles[organ_b][ia]
                - profiles[organ_a][ia] - profiles[organ_b][ib]
            )
            a, b = sorted((organ_a, organ_b))
            records.append({"organ_a": a, "organ_b": b, "swap_score": float(gain)})
    records.sort(key=lambda row: (-row["swap_score"], row["organ_a"], row["organ_b"]))
    for rank, record in enumerate(records, 1):
        record["rank"] = rank
    return records


def jackknife_stability(vectors, expected_pair):
    """Remove each marker-locus vector once and re-rank all candidate swaps."""
    margins = []
    wins = 0
    trials = 0
    for tissue in ORGANS:
        if len(vectors[tissue]) <= 1:
            continue
        for drop in range(len(vectors[tissue])):
            reduced = {name: list(rows) for name, rows in vectors.items()}
            reduced[tissue] = [row for i, row in enumerate(vectors[tissue]) if i != drop]
            ranked = rank_swaps(mean_profiles(reduced))
            pair = (ranked[0]["organ_a"], ranked[0]["organ_b"])
            if pair == expected_pair:
                wins += 1
            margins.append(ranked[0]["swap_score"] - ranked[1]["swap_score"])
            trials += 1
    return {
        "trials": trials,
        "top_pair_fraction": wins / trials if trials else 0.0,
        "minimum_top_margin": min(margins) if margins else 0.0,
        "median_top_margin": statistics.median(margins) if margins else 0.0,
    }


def main():
    script = Path(__file__).resolve()
    repo = next(parent for parent in script.parents if parent.name == "lifescience-benchmark-827")
    input_dir = repo / "inputs" / "ls03-atac-sample-swap"
    out_dir = script.parent

    _, by_arm = read_contigs(input_dir / "sample.swap.atac.q1.chrom.sizes")
    loci = collect_marker_tss(input_dir / "AmexT_v47-AmexG_v6.0-DD.gtf.gz", by_arm)
    desired, marker_bin_sets = desired_promoter_bins(loci)
    totals, captured, row_count, widths = read_atac(
        input_dir / "sample.swap.atac.q1.tsv.gz", desired
    )
    vectors, matched_genes = standardized_promoter_vectors(marker_bin_sets, captured, totals)
    profiles = mean_profiles(vectors)
    ranked = rank_swaps(profiles)
    top, runner_up = ranked[0], ranked[1]
    top_pair = (top["organ_a"], top["organ_b"])
    margin = top["swap_score"] - runner_up["swap_score"]
    stability = jackknife_stability(vectors, top_pair)

    unique = margin > 0.5 and stability["top_pair_fraction"] >= 0.95
    confidence = "high" if unique and stability["minimum_top_margin"] > 0 else "moderate" if unique else "low"

    for row in ranked:
        row["evidence_type"] = "promoter-marker coherence gain; leave-one-marker-locus-out stability"
    with (out_dir / "sample_similarity.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["organ_a", "organ_b", "swap_score", "rank", "evidence_type"],
            lineterminator="\n",
        )
        writer.writeheader()
        for row in ranked:
            writer.writerow({**row, "swap_score": f'{row["swap_score"]:.12f}'})

    ia, ib = ORGANS.index(top["organ_a"]), ORGANS.index(top["organ_b"])
    evidence = [
        (
            f"The strongest reciprocal promoter-marker coherence gain is {top['swap_score']:.6f} "
            f"for {top['organ_a']} <-> {top['organ_b']}; the runner-up is "
            f"{runner_up['organ_a']} <-> {runner_up['organ_b']} at "
            f"{runner_up['swap_score']:.6f} (margin {margin:.6f})."
        ),
        (
            f"{top['organ_a']} markers score {profiles[top['organ_a']][ib]:.4f} in the "
            f"{top['organ_b']}-labelled sample versus {profiles[top['organ_a']][ia]:.4f} in "
            f"the {top['organ_a']}-labelled sample; {top['organ_b']} markers score "
            f"{profiles[top['organ_b']][ia]:.4f} versus {profiles[top['organ_b']][ib]:.4f}."
        ),
        (
            f"Leave-one-marker-locus-out re-ranking retained the same top pair in "
            f"{stability['top_pair_fraction']:.1%} of {stability['trials']} trials; "
            f"minimum top-versus-runner-up margin was {stability['minimum_top_margin']:.6f}."
        ),
        "Signals were normalized by complete-library CPM and expressed relative to each promoter locus's across-organ mean; total library size alone was not used as swap evidence.",
    ]
    call = {
        "swap_detected": bool(unique),
        "organ_a": top["organ_a"] if unique else None,
        "organ_b": top["organ_b"] if unique else None,
        "confidence": confidence,
        "evidence": evidence,
    }
    with (out_dir / "swap_call.json").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(call, indent=2) + "\n")

    marker_summary = ", ".join(
        f"{tissue}={len(matched_genes[tissue])} genes/{len(vectors[tissue])} loci"
        for tissue in ORGANS
    )
    report = f"""# ATAC-seq organ-label swap analysis

## Result

**Swap detected: {'yes' if unique else 'no'}**

**Pair: {top['organ_a']} ↔ {top['organ_b']}**

**Confidence: {confidence}**

The top swap score is `{top['swap_score']:.6f}`. The runner-up is
{runner_up['organ_a']} ↔ {runner_up['organ_b']} at `{runner_up['swap_score']:.6f}`,
leaving a margin of `{margin:.6f}`. The leading pair is therefore distinct rather
than one of several nearly tied alternatives.

## Evidence

- {evidence[1]}
- {evidence[2]}
- {evidence[3]}

## Method

The GTF gene TSS coordinates were mapped from whole chromosome arms onto the
split-contig coordinate system in the supplied chromosome-size file. For a
predeclared panel of organ identity markers, accessibility was summed over the TSS
bin and one adjacent 10-kb bin on each side. Each sample was normalized by its full
ATAC library count (CPM), then divided by that promoter locus's mean across the 15
organs so that each locus contributed a comparable relative-accessibility profile.
For every one of the 105 unordered organ pairs,
the `swap_score` is the gain in reciprocal marker-to-label coherence after exchanging
that pair. A deterministic leave-one-marker-locus-out analysis re-ranked all pairs
to assess robustness; it is an effect-size stability check, not a p-value.

The table contained {row_count:,} data rows with observed bin widths {widths};
{len(captured):,} requested promoter-region bins were present. Marker coverage:
{marker_summary}.

No figure was generated because the requested ranked CSV directly represents the
complete 105-pair comparison without loss of information.

## Reproduction

Run `python analysis.py` from any directory. The script locates the repository from
its own path and rewrites `swap_call.json`, `sample_similarity.csv`, and this report
deterministically.
"""
    with (out_dir / "report.md").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(report)

    print(json.dumps({
        "top_pair": top_pair,
        "top_score": top["swap_score"],
        "runner_up": (runner_up["organ_a"], runner_up["organ_b"]),
        "margin": margin,
        "stability": stability,
        "swap_detected": unique,
        "confidence": confidence,
    }, indent=2))


if __name__ == "__main__":
    main()
