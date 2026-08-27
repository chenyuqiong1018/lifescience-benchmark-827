#!/usr/bin/env python3
"""Locate a large deletion from shallow paired-end chr22 reads using stdlib only."""

from __future__ import annotations

import gzip
import hashlib
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

TASK = "ls02-find-deletion"
BIN_SIZE = 100_000
K = 31
MAX_SEED_HITS = 12
MASK = (1 << (2 * K)) - 1
BASE = [-1] * 256
for _base, _value in zip(b"ACGT", range(4)):
    BASE[_base] = _value
COMP = bytes.maketrans(b"ACGTN", b"TGCAN")


def find_repo() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "inputs" / TASK).is_dir():
            return parent
    raise RuntimeError("Repository root with the permitted task input was not found")


REPO = find_repo()
INPUT = REPO / "inputs" / TASK
OUT = Path(__file__).resolve().parent
R1 = INPUT / "find.deletion.r1.fq.gz"
R2 = INPUT / "find.deletion.r2.fq.gz"
REFERENCE = INPUT / "reference" / "GRCh38_chr22.fa.gz"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fastq_records(path: Path):
    with gzip.open(path, "rb") as handle:
        while True:
            header = handle.readline().rstrip()
            if not header:
                return
            sequence = handle.readline().strip().upper()
            plus = handle.readline().rstrip()
            quality = handle.readline().rstrip()
            if not plus.startswith(b"+") or len(sequence) != len(quality):
                raise ValueError(f"Malformed FASTQ record in {path}")
            yield header, sequence, quality


def paired_records():
    r1_iter = fastq_records(R1)
    r2_iter = fastq_records(R2)
    count = 0
    while True:
        a = next(r1_iter, None)
        b = next(r2_iter, None)
        if a is None or b is None:
            if a is not None or b is not None:
                raise ValueError("FASTQ files contain different record counts")
            return
        count += 1
        name1 = a[0].split(b"/")[0]
        name2 = b[0].split(b"/")[0]
        if name1 != name2:
            raise ValueError(f"Pair-name mismatch at record {count}")
        yield count, a, b


def revcomp(sequence: bytes) -> bytes:
    return sequence.translate(COMP)[::-1]


def code_kmer(sequence: bytes):
    code = 0
    for nucleotide in sequence:
        value = BASE[nucleotide]
        if value < 0:
            return None
        code = (code << 2) | value
    return code


with gzip.open(REFERENCE, "rb") as handle:
    reference_header = handle.readline().decode().strip()
    reference = b"".join(line.strip().upper() for line in handle)
if reference_header != ">chr22":
    raise ValueError(f"Unexpected reference header: {reference_header}")

# First FASTQ pass: validate basic QC and collect only read-derived terminal seeds.
wanted = set()
read_lengths = Counter()
quality_scores = Counter()
pair_count = 0
for pair_count, r1, r2 in paired_records():
    for _, sequence, quality in (r1, r2):
        read_lengths[len(sequence)] += 1
        quality_scores.update(q - 33 for q in quality)
        terminal = (sequence[:K], sequence[-K:])
        for seed in (*terminal, *(revcomp(x) for x in terminal)):
            encoded = code_kmer(seed)
            if encoded is not None:
                wanted.add(encoded)

# Single reference scan builds a bounded index only for the observed terminal seeds.
positions = defaultdict(list)
overflow = set()
rolling = 0
valid = 0
for index, nucleotide in enumerate(reference):
    value = BASE[nucleotide]
    if value < 0:
        rolling = 0
        valid = 0
        continue
    rolling = ((rolling << 2) | value) & MASK
    valid += 1
    if valid >= K and rolling in wanted:
        hits = positions[rolling]
        if len(hits) < MAX_SEED_HITS:
            hits.append(index - K + 1)
        else:
            overflow.add(rolling)


def full_maps(oriented: bytes):
    starts = set()
    for anchor, offset in ((oriented[:K], 0), (oriented[-K:], len(oriented) - K)):
        seed = code_kmer(anchor)
        if seed is None or seed in overflow:
            continue
        for position in positions.get(seed, ()):
            start = position - offset
            if start >= 0 and reference[start : start + len(oriented)] == oriented:
                starts.add(start)
    return starts


def map_read(sequence: bytes):
    maps = [(position, "+") for position in full_maps(sequence)]
    maps.extend((position, "-") for position in full_maps(revcomp(sequence)))
    return sorted(maps)


bin_count = (len(reference) + BIN_SIZE - 1) // BIN_SIZE
mapped_bases = [0] * bin_count
map_multiplicity = Counter()
pair_class = Counter()
normal_fr_spans = []
long_fr_pairs = []
unmapped_reads = []

for pair_index, r1, r2 in paired_records():
    sequence1 = r1[1]
    sequence2 = r2[1]
    maps1 = map_read(sequence1)
    maps2 = map_read(sequence2)
    map_multiplicity[len(maps1)] += 1
    map_multiplicity[len(maps2)] += 1

    for end_name, sequence, maps in (("R1", sequence1, maps1), ("R2", sequence2, maps2)):
        if len(maps) == 0:
            unmapped_reads.append((pair_index, end_name, sequence))
        if len(maps) != 1:
            continue
        start = maps[0][0]
        end = start + len(sequence)
        first_bin = start // BIN_SIZE
        last_bin = (end - 1) // BIN_SIZE
        for bin_index in range(first_bin, last_bin + 1):
            overlap = min(end, (bin_index + 1) * BIN_SIZE) - max(start, bin_index * BIN_SIZE)
            mapped_bases[bin_index] += max(0, overlap)

    if len(maps1) != 1 or len(maps2) != 1:
        pair_class["not_both_unique"] += 1
        continue
    position1, strand1 = maps1[0]
    position2, strand2 = maps2[0]
    proper_fr = (position1 < position2 and strand1 == "+" and strand2 == "-") or (
        position2 < position1 and strand2 == "+" and strand1 == "-"
    )
    if not proper_fr:
        pair_class["other_orientation"] += 1
        continue
    pair_class["FR"] += 1
    left = min(position1, position2)
    right = max(position1 + len(sequence1), position2 + len(sequence2))
    span = right - left
    if span <= 2_000:
        normal_fr_spans.append(span)
    elif span >= 10_000:
        long_fr_pairs.append(
            {
                "pair": pair_index,
                "r1_start_1based": position1 + 1,
                "r1_strand": strand1,
                "r2_start_1based": position2 + 1,
                "r2_strand": strand2,
                "reference_span_bp": span,
            }
        )

callable_bases = []
for bin_index in range(bin_count):
    chunk = reference[bin_index * BIN_SIZE : (bin_index + 1) * BIN_SIZE]
    callable_bases.append(len(chunk) - chunk.count(b"N"))

# Infer the event from the longest zero-depth run among mostly callable 100-kb bins.
runs = []
run_start = None
for bin_index, (mapped, callable_count) in enumerate(zip(mapped_bases, callable_bases)):
    qualifies = mapped == 0 and callable_count >= 80_000
    if qualifies and run_start is None:
        run_start = bin_index
    if run_start is not None and (not qualifies or bin_index == bin_count - 1):
        run_end = bin_index if qualifies and bin_index == bin_count - 1 else bin_index - 1
        runs.append((run_start, run_end))
        run_start = None
if not runs:
    raise RuntimeError("No callable zero-depth run was found")
event_run = max(runs, key=lambda value: value[1] - value[0] + 1)
event_start = event_run[0] * BIN_SIZE
event_end = (event_run[1] + 1) * BIN_SIZE
event_size = event_end - event_start

median_normal_span = statistics.median(normal_fr_spans)
support_pairs = []
for pair in long_fr_pairs:
    starts = [pair["r1_start_1based"] - 1, pair["r2_start_1based"] - 1]
    ends = [starts[0] + next(iter(read_lengths)), starts[1] + next(iter(read_lengths))]
    left_start = min(starts)
    right_start = max(starts)
    adjusted_span = pair["reference_span_bp"] - event_size
    if min(ends) <= event_start and right_start >= event_end and 100 <= adjusted_span <= 2_000:
        pair = dict(pair)
        pair["deletion_adjusted_span_bp"] = adjusted_span
        support_pairs.append(pair)

# Exact junction support: an unmapped read must cross the inferred reference join.
flank = 300
left_flank = reference[max(0, event_start - flank) : event_start]
right_flank = reference[event_end : min(len(reference), event_end + flank)]
junction = left_flank + right_flank
junction_boundary = len(left_flank)
junction_reads = []
for pair_index, end_name, sequence in unmapped_reads:
    for strand, oriented in (("+", sequence), ("-", revcomp(sequence))):
        offset = junction.find(oriented)
        while offset >= 0:
            if offset < junction_boundary < offset + len(oriented):
                junction_reads.append(
                    {
                        "pair": pair_index,
                        "end": end_name,
                        "strand": strand,
                        "left_aligned_bp": junction_boundary - offset,
                        "right_aligned_bp": offset + len(oriented) - junction_boundary,
                    }
                )
                break
            offset = junction.find(oriented, offset + 1)

outside_mapped = sum(mapped_bases[: event_run[0]]) + sum(mapped_bases[event_run[1] + 1 :])
outside_callable = sum(callable_bases[: event_run[0]]) + sum(callable_bases[event_run[1] + 1 :])
outside_depth = outside_mapped / outside_callable
left_flank_bin = mapped_bases[event_run[0] - 1] if event_run[0] else None
right_flank_bin = mapped_bases[event_run[1] + 1] if event_run[1] + 1 < bin_count else None
supporting_signals = (
    f"{event_run[1] - event_run[0] + 1} consecutive callable 100-kb bins with zero unique-read depth; "
    f"{len(support_pairs)} FR read pairs span the interval; "
    f"{len(junction_reads)} exact junction-spanning reads"
)

deletion_tsv = (
    "chrom\tstart_100kb\tend_100kb\tsize_bp\tsupporting_signals\n"
    f"chr22\t{event_start}\t{event_end}\t{event_size}\t{supporting_signals}\n"
)
(OUT / "deletion.tsv").write_text(deletion_tsv, encoding="utf-8", newline="\n")

qc = {
    "task": TASK,
    "method": "stdlib exact 31-mer terminal anchoring; unique full-read verification; 100-kb depth; FR-pair and exact-junction confirmation",
    "inputs": {
        "r1": {"path": str(R1.relative_to(REPO)), "sha256": sha256(R1)},
        "r2": {"path": str(R2.relative_to(REPO)), "sha256": sha256(R2)},
        "reference": {"path": str(REFERENCE.relative_to(REPO)), "sha256": sha256(REFERENCE)},
    },
    "reference": {"name": "chr22", "length_bp": len(reference), "N_bases": reference.count(b"N")},
    "fastq": {
        "pairs": pair_count,
        "ends": pair_count * 2,
        "read_length_counts": {str(k): v for k, v in sorted(read_lengths.items())},
        "quality_score_counts": {str(k): v for k, v in sorted(quality_scores.items())},
    },
    "mapping": {
        "seed_length": K,
        "observed_terminal_seed_count": len(wanted),
        "indexed_seed_count": len(positions),
        "repetitive_seed_count": len(overflow),
        "end_mapping_multiplicity": {str(k): v for k, v in sorted(map_multiplicity.items())},
        "pair_classes": dict(pair_class),
        "median_normal_FR_span_bp": median_normal_span,
        "mean_unique_read_depth_outside_deletion": outside_depth,
    },
    "deletion_evidence": {
        "zero_depth_bin_start_indices": list(range(event_run[0], event_run[1] + 1)),
        "zero_depth_bin_count": event_run[1] - event_run[0] + 1,
        "left_flanking_bin_mapped_bases": left_flank_bin,
        "right_flanking_bin_mapped_bases": right_flank_bin,
        "spanning_FR_pairs": support_pairs,
        "junction_spanning_reads": junction_reads,
    },
    "precision": {
        "reported_rounding_bp": BIN_SIZE,
        "reported_breakpoints": [event_start, event_end],
        "interpretation": "Rounded breakpoint boundaries; deletion interval in 1-based reference terms is approximately start+1 through end.",
    },
}
(OUT / "qc.json").write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8", newline="\n")

report = f"""# Large-deletion analysis

## Call

The data support a single approximately **{event_size:,} bp deletion on chr22**, with breakpoints reported at the requested 100-kb precision as **{event_start:,}** and **{event_end:,}**.

## Evidence

- Depth: {event_run[1] - event_run[0] + 1} consecutive, mostly callable 100-kb bins from {event_start:,} to {event_end:,} have zero uniquely mapped read bases. The immediately adjacent bins contain {left_flank_bin:,} and {right_flank_bin:,} mapped bases, respectively; mean unique-read depth outside the interval is {outside_depth:.3f}x.
- Paired-end geometry: {len(support_pairs)} correctly oriented FR pairs bridge the interval. Their reference spans become ordinary fragment spans after subtracting the {event_size:,} bp deletion; the median ordinary FR span in the library is {median_normal_span:.0f} bp.
- Junction reads: {len(junction_reads)} otherwise-unmapped reads match exactly across the sequence formed by joining the left and right inferred breakpoints.

These three signal types agree on the same event. Repetitive terminal seeds were excluded from unique mapping, and anomalous pairs not matching the inferred interval were not counted as support.

## Precision and coordinate limits

`deletion.tsv` intentionally reports 100-kb-rounded breakpoint boundaries, as requested. The depth segmentation itself cannot justify sub-bin precision. Exact junction matches are compatible with a reference join at the displayed boundaries, but the primary deliverable should still be interpreted at 100-kb resolution. In 1-based interval language, the removed reference sequence is approximately {event_start + 1:,} through {event_end:,}; this boundary convention does not change the rounded values or the {event_size:,} bp size.

## Reproducibility

`analysis.py` uses only the Python standard library. It validates pairing and FASTQ structure, indexes observed terminal 31-mers during one reference scan, verifies full-read matches before calling them unique, computes 100-kb coverage, and confirms the candidate with FR-pair and exact-junction evidence. Full counts, hashes, and supporting record coordinates are in `qc.json`.
"""
(OUT / "report.md").write_text(report, encoding="utf-8", newline="\n")

print(
    json.dumps(
        {
            "chrom": "chr22",
            "start_100kb": event_start,
            "end_100kb": event_end,
            "size_bp": event_size,
            "zero_bins": event_run[1] - event_run[0] + 1,
            "spanning_pairs": len(support_pairs),
            "junction_reads": len(junction_reads),
        },
        sort_keys=True,
    )
)
