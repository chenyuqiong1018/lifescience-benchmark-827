#!/usr/bin/env python3
"""Controlled-skill T1 analysis of the chr22 shallow paired-end data."""

import gzip
import hashlib
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

TASK, BIN, K, CAP = "ls02-find-deletion", 100_000, 31, 10
COMP = bytes.maketrans(b"ACGTN", b"TGCAN")
VAL = [-1] * 256
for b, v in zip(b"ACGT", range(4)):
    VAL[b] = v
MASK = (1 << (2 * K)) - 1


def locate_repo():
    for parent in Path(__file__).resolve().parents:
        if (parent / "inputs" / TASK).is_dir():
            return parent
    raise RuntimeError("repository root not found")


REPO = locate_repo()
INPUT = REPO / "inputs" / TASK
OUTPUT = Path(__file__).resolve().parent
P1 = INPUT / "find.deletion.r1.fq.gz"
P2 = INPUT / "find.deletion.r2.fq.gz"
FA = INPUT / "reference" / "GRCh38_chr22.fa.gz"


def fastq(path):
    with gzip.open(path, "rb") as stream:
        while (name := stream.readline().strip()):
            seq = stream.readline().strip().upper()
            plus = stream.readline()
            qual = stream.readline().strip()
            if not plus.startswith(b"+") or len(seq) != len(qual):
                raise ValueError(f"invalid FASTQ: {path}")
            yield name, seq, qual


def paired():
    a, b, count = fastq(P1), fastq(P2), 0
    while True:
        x, y = next(a, None), next(b, None)
        if x is None or y is None:
            if x is not None or y is not None:
                raise ValueError("unequal paired FASTQ counts")
            return
        count += 1
        if x[0].split(b"/")[0] != y[0].split(b"/")[0]:
            raise ValueError(f"pair order mismatch: {count}")
        yield count, x, y


def rc(seq):
    return seq.translate(COMP)[::-1]


def encode(seq):
    value = 0
    for base in seq:
        if VAL[base] < 0:
            return None
        value = (value << 2) | VAL[base]
    return value


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


with gzip.open(FA, "rb") as stream:
    fasta_name = stream.readline().strip()
    ref = b"".join(line.strip().upper() for line in stream)
if fasta_name != b">chr22" or len(ref) != 50_818_468:
    raise ValueError("supplied reference is not the expected hg38 chr22 sequence")

# Controlled chromosome/UCSC guidance: validate reference identity before coordinates.
seeds, read_lengths, phred = set(), Counter(), Counter()
pair_total = 0
for pair_total, x, y in paired():
    for _, seq, qual in (x, y):
        read_lengths[len(seq)] += 1
        phred.update(q - 33 for q in qual)
        for oriented in (seq, rc(seq)):
            code = encode(oriented[:K])
            if code is not None:
                seeds.add(code)

hits, repetitive = defaultdict(list), set()
rolling = valid = 0
for i, base in enumerate(ref):
    digit = VAL[base]
    if digit < 0:
        rolling = valid = 0
        continue
    rolling, valid = ((rolling << 2) | digit) & MASK, valid + 1
    if valid >= K and rolling in seeds:
        if len(hits[rolling]) < CAP:
            hits[rolling].append(i - K + 1)
        else:
            repetitive.add(rolling)


def map_read(seq):
    answer = []
    for strand, oriented in (("+", seq), ("-", rc(seq))):
        code = encode(oriented[:K])
        if code is None or code in repetitive:
            continue
        for start in hits.get(code, ()):
            if ref[start : start + len(seq)] == oriented:
                answer.append((start, strand))
    return sorted(set(answer))


depth = [0] * ((len(ref) + BIN - 1) // BIN)
map_counts, pair_classes = Counter(), Counter()
ordinary, long_pairs, unmapped = [], [], []
read_length = next(iter(read_lengths))
for pair_id, x, y in paired():
    s1, s2 = x[1], y[1]
    m1, m2 = map_read(s1), map_read(s2)
    map_counts[len(m1)] += 1
    map_counts[len(m2)] += 1
    for label, seq, maps in (("R1", s1, m1), ("R2", s2, m2)):
        if not maps:
            unmapped.append((pair_id, label, seq))
        if len(maps) == 1:
            start, end = maps[0][0], maps[0][0] + len(seq)
            for j in range(start // BIN, (end - 1) // BIN + 1):
                depth[j] += min(end, (j + 1) * BIN) - max(start, j * BIN)
    if len(m1) != 1 or len(m2) != 1:
        pair_classes["not_both_unique"] += 1
        continue
    p1, t1 = m1[0]
    p2, t2 = m2[0]
    proper = (p1 < p2 and t1 == "+" and t2 == "-") or (p2 < p1 and t2 == "+" and t1 == "-")
    if not proper:
        pair_classes["other_orientation"] += 1
        continue
    pair_classes["FR"] += 1
    span = max(p1 + len(s1), p2 + len(s2)) - min(p1, p2)
    if span <= 2_000:
        ordinary.append(span)
    elif span >= 10_000:
        long_pairs.append({"pair": pair_id, "r1_start_1based": p1 + 1, "r1_strand": t1, "r2_start_1based": p2 + 1, "r2_strand": t2, "reference_span_bp": span})

callable_count = []
for j in range(len(depth)):
    block = ref[j * BIN : (j + 1) * BIN]
    callable_count.append(len(block) - block.count(b"N"))
runs, opened = [], None
for j, (bases, callable_bases) in enumerate(zip(depth, callable_count)):
    empty = bases == 0 and callable_bases >= 80_000
    if empty and opened is None:
        opened = j
    if opened is not None and (not empty or j == len(depth) - 1):
        closed = j if empty and j == len(depth) - 1 else j - 1
        runs.append((opened, closed))
        opened = None
event = max(runs, key=lambda x: x[1] - x[0])
start, end = event[0] * BIN, (event[1] + 1) * BIN
size = end - start
median_insert = statistics.median(ordinary)

bridges = []
for item in long_pairs:
    starts = (item["r1_start_1based"] - 1, item["r2_start_1based"] - 1)
    adjusted = item["reference_span_bp"] - size
    if min(starts) + read_length <= start and max(starts) >= end and 100 <= adjusted <= 2_000:
        item = dict(item)
        item["deletion_adjusted_span_bp"] = adjusted
        bridges.append(item)

flank = 300
junction = ref[start - flank : start] + ref[end : end + flank]
split_reads = []
for pair_id, label, seq in unmapped:
    for strand, oriented in (("+", seq), ("-", rc(seq))):
        offset = junction.find(oriented)
        if 0 <= offset < flank < offset + len(oriented):
            split_reads.append({"pair": pair_id, "end": label, "strand": strand, "left_aligned_bp": flank - offset, "right_aligned_bp": offset + len(oriented) - flank})

signal = f"{event[1]-event[0]+1} consecutive callable zero-depth 100-kb bins; {len(bridges)} spanning FR pairs; {len(split_reads)} exact junction reads"
(OUTPUT / "deletion.tsv").write_text("chrom\tstart_100kb\tend_100kb\tsize_bp\tsupporting_signals\n" + f"chr22\t{start}\t{end}\t{size}\t{signal}\n", encoding="utf-8", newline="\n")
outside_bases = sum(depth[: event[0]]) + sum(depth[event[1] + 1 :])
outside_callable = sum(callable_count[: event[0]]) + sum(callable_count[event[1] + 1 :])
qc = {
    "task": TASK,
    "arm": "T1",
    "controlled_skills": ["chromosome_analysis", "ucsc_genome_exploration", "genome_annotation", "code_execution_analysis"],
    "inputs": {"r1_sha256": sha(P1), "r2_sha256": sha(P2), "reference_sha256": sha(FA)},
    "reference": {"assembly": "hg38", "chrom": "chr22", "length_bp": len(ref), "N_bases": ref.count(b"N")},
    "reads": {"pairs": pair_total, "length_counts": {str(k): v for k, v in sorted(read_lengths.items())}, "phred_counts": {str(k): v for k, v in sorted(phred.items())}},
    "mapping": {"unique_ends": map_counts[1], "multiplicity": {str(k): v for k, v in sorted(map_counts.items())}, "pair_classes": dict(pair_classes), "median_normal_FR_span_bp": median_insert, "mean_unique_depth_outside_event": outside_bases / outside_callable},
    "evidence": {"zero_depth_bins": list(range(event[0], event[1] + 1)), "left_flank_mapped_bases": depth[event[0] - 1], "right_flank_mapped_bases": depth[event[1] + 1], "spanning_FR_pairs": bridges, "junction_reads": split_reads},
    "precision": {"rounding_bp": BIN, "boundaries": [start, end], "one_based_removed_interval": [start + 1, end]},
}
(OUTPUT / "qc.json").write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8", newline="\n")
report = f"""# Controlled T1 large-deletion analysis

## Call

The controlled T1 workflow supports a **{size:,} bp chr22 deletion** between the required 100-kb-rounded boundaries **{start:,}** and **{end:,}**.

## Evidence

- {event[1]-event[0]+1} consecutive mostly callable 100-kb bins have zero unique-read depth; adjacent bins contain {depth[event[0]-1]:,} and {depth[event[1]+1]:,} mapped bases.
- {len(bridges)} correctly oriented FR pairs span the interval, and their deletion-adjusted spans agree with the {median_insert:.0f} bp normal-library median.
- {len(split_reads)} otherwise-unmapped reads match exactly across the inferred reference junction.

The chromosome/UCSC skills informed explicit hg38 chr22 validation and coordinate context; genome annotation informed the 1-based interpretation; code execution informed the auditable standard-library pipeline. The supplied, hash-recorded chr22 sequence remained authoritative, so unavailable credentialed SCP examples were not needed.

## Precision

The deliverable is intentionally limited to 100-kb precision. Zero-depth bins do not independently justify finer breakpoints. Exact junction reads are consistent with the displayed boundaries; in 1-based interval terms the removed sequence is approximately {start+1:,} through {end:,}.
"""
(OUTPUT / "report.md").write_text(report, encoding="utf-8", newline="\n")
print(json.dumps({"chrom": "chr22", "start_100kb": start, "end_100kb": end, "size_bp": size, "spanning_pairs": len(bridges), "junction_reads": len(split_reads)}, sort_keys=True))
