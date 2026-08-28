#!/usr/bin/env python3
"""Independent T0 large-deletion caller using exact read anchors and 100-kb depth."""

import gzip
import hashlib
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

TASK = "ls02-find-deletion"
BIN = 100_000
K = 31
MAX_HITS = 10
COMP = bytes.maketrans(b"ACGTN", b"TGCAN")
ENC = [-1] * 256
for base, value in zip(b"ACGT", range(4)):
    ENC[base] = value
MASK = (1 << (2 * K)) - 1


def repo_root():
    for parent in Path(__file__).resolve().parents:
        if (parent / "inputs" / TASK).is_dir():
            return parent
    raise RuntimeError("task input not found")


REPO = repo_root()
INP = REPO / "inputs" / TASK
OUT = Path(__file__).resolve().parent
R1 = INP / "find.deletion.r1.fq.gz"
R2 = INP / "find.deletion.r2.fq.gz"
REF = INP / "reference" / "GRCh38_chr22.fa.gz"


def records(path):
    with gzip.open(path, "rb") as handle:
        while (header := handle.readline().rstrip()):
            sequence = handle.readline().strip().upper()
            plus = handle.readline()
            quality = handle.readline().strip()
            if not plus.startswith(b"+") or len(sequence) != len(quality):
                raise ValueError(f"malformed FASTQ: {path}")
            yield header, sequence, quality


def pairs():
    left, right = records(R1), records(R2)
    index = 0
    while True:
        a, b = next(left, None), next(right, None)
        if a is None or b is None:
            if a is not None or b is not None:
                raise ValueError("unequal FASTQ record counts")
            return
        index += 1
        if a[0].split(b"/")[0] != b[0].split(b"/")[0]:
            raise ValueError(f"pair mismatch at {index}")
        yield index, a, b


def reverse_complement(sequence):
    return sequence.translate(COMP)[::-1]


def encode(sequence):
    code = 0
    for nucleotide in sequence:
        value = ENC[nucleotide]
        if value < 0:
            return None
        code = (code << 2) | value
    return code


def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


with gzip.open(REF, "rb") as handle:
    header = handle.readline().strip()
    reference = b"".join(line.strip().upper() for line in handle)
if header != b">chr22":
    raise ValueError(f"expected chr22, found {header!r}")

# Index only read-derived prefix anchors, including reverse-complement orientation.
wanted = set()
lengths, qualities = Counter(), Counter()
pair_count = 0
for pair_count, first, second in pairs():
    for _, sequence, quality in (first, second):
        lengths[len(sequence)] += 1
        qualities.update(q - 33 for q in quality)
        for oriented in (sequence, reverse_complement(sequence)):
            seed = encode(oriented[:K])
            if seed is not None:
                wanted.add(seed)

index = defaultdict(list)
repetitive = set()
rolling = valid = 0
for position, nucleotide in enumerate(reference):
    value = ENC[nucleotide]
    if value < 0:
        rolling = valid = 0
        continue
    rolling = ((rolling << 2) | value) & MASK
    valid += 1
    if valid >= K and rolling in wanted:
        hits = index[rolling]
        if len(hits) < MAX_HITS:
            hits.append(position - K + 1)
        else:
            repetitive.add(rolling)


def map_end(sequence):
    result = []
    for strand, oriented in (("+", sequence), ("-", reverse_complement(sequence))):
        seed = encode(oriented[:K])
        if seed is None or seed in repetitive:
            continue
        for start in index.get(seed, ()):
            if reference[start : start + len(oriented)] == oriented:
                result.append((start, strand))
    return sorted(set(result))


depth_bases = [0] * ((len(reference) + BIN - 1) // BIN)
mapping_counts = Counter()
pair_classes = Counter()
normal_spans, large_pairs, unmapped = [], [], []
read_length = next(iter(lengths))

for pair_index, first, second in pairs():
    seq1, seq2 = first[1], second[1]
    maps1, maps2 = map_end(seq1), map_end(seq2)
    mapping_counts[len(maps1)] += 1
    mapping_counts[len(maps2)] += 1
    for end_name, sequence, maps in (("R1", seq1, maps1), ("R2", seq2, maps2)):
        if not maps:
            unmapped.append((pair_index, end_name, sequence))
        if len(maps) == 1:
            start, end = maps[0][0], maps[0][0] + len(sequence)
            for bin_id in range(start // BIN, (end - 1) // BIN + 1):
                depth_bases[bin_id] += min(end, (bin_id + 1) * BIN) - max(start, bin_id * BIN)
    if len(maps1) != 1 or len(maps2) != 1:
        pair_classes["not_both_unique"] += 1
        continue
    p1, s1 = maps1[0]
    p2, s2 = maps2[0]
    proper = (p1 < p2 and s1 == "+" and s2 == "-") or (p2 < p1 and s2 == "+" and s1 == "-")
    if not proper:
        pair_classes["other_orientation"] += 1
        continue
    pair_classes["FR"] += 1
    span = max(p1 + len(seq1), p2 + len(seq2)) - min(p1, p2)
    if span <= 2_000:
        normal_spans.append(span)
    elif span >= 10_000:
        large_pairs.append({"pair": pair_index, "r1_start_1based": p1 + 1, "r1_strand": s1, "r2_start_1based": p2 + 1, "r2_strand": s2, "reference_span_bp": span})

callable_bases = []
for bin_id in range(len(depth_bases)):
    sequence = reference[bin_id * BIN : (bin_id + 1) * BIN]
    callable_bases.append(len(sequence) - sequence.count(b"N"))

runs, start = [], None
for bin_id, (depth, callable_count) in enumerate(zip(depth_bases, callable_bases)):
    empty_callable = depth == 0 and callable_count >= 80_000
    if empty_callable and start is None:
        start = bin_id
    if start is not None and (not empty_callable or bin_id == len(depth_bases) - 1):
        end = bin_id if empty_callable and bin_id == len(depth_bases) - 1 else bin_id - 1
        runs.append((start, end))
        start = None
event_bins = max(runs, key=lambda item: item[1] - item[0])
event_start, event_end = event_bins[0] * BIN, (event_bins[1] + 1) * BIN
event_size = event_end - event_start
median_span = statistics.median(normal_spans)

support_pairs = []
for item in large_pairs:
    starts = (item["r1_start_1based"] - 1, item["r2_start_1based"] - 1)
    adjusted = item["reference_span_bp"] - event_size
    if min(starts) + read_length <= event_start and max(starts) >= event_end and 100 <= adjusted <= 2_000:
        row = dict(item)
        row["deletion_adjusted_span_bp"] = adjusted
        support_pairs.append(row)

flank = 300
junction = reference[event_start - flank : event_start] + reference[event_end : event_end + flank]
junction_reads = []
for pair_index, end_name, sequence in unmapped:
    for strand, oriented in (("+", sequence), ("-", reverse_complement(sequence))):
        offset = junction.find(oriented)
        if 0 <= offset < flank < offset + len(oriented):
            junction_reads.append({"pair": pair_index, "end": end_name, "strand": strand, "left_aligned_bp": flank - offset, "right_aligned_bp": offset + len(oriented) - flank})

signals = f"{event_bins[1]-event_bins[0]+1} consecutive callable 100-kb zero-depth bins; {len(support_pairs)} spanning FR pairs; {len(junction_reads)} exact junction reads"
(OUT / "deletion.tsv").write_text(
    "chrom\tstart_100kb\tend_100kb\tsize_bp\tsupporting_signals\n"
    f"chr22\t{event_start}\t{event_end}\t{event_size}\t{signals}\n",
    encoding="utf-8",
    newline="\n",
)

outside_mapped = sum(depth_bases[: event_bins[0]]) + sum(depth_bases[event_bins[1] + 1 :])
outside_callable = sum(callable_bases[: event_bins[0]]) + sum(callable_bases[event_bins[1] + 1 :])
qc = {
    "task": TASK,
    "arm": "T0",
    "skills_used": ["chromosome_analysis", "genome_annotation", "code_execution_analysis"],
    "input_sha256": {"r1": file_hash(R1), "r2": file_hash(R2), "reference": file_hash(REF)},
    "reference": {"name": "chr22", "length_bp": len(reference), "N_bases": reference.count(b"N")},
    "fastq": {"pairs": pair_count, "read_length_counts": {str(k): v for k, v in sorted(lengths.items())}, "quality_score_counts": {str(k): v for k, v in sorted(qualities.items())}},
    "mapping": {"anchor_bp": K, "unique_end_count": mapping_counts[1], "mapping_multiplicity": {str(k): v for k, v in sorted(mapping_counts.items())}, "pair_classes": dict(pair_classes), "median_normal_FR_span_bp": median_span, "mean_unique_depth_outside_event": outside_mapped / outside_callable},
    "evidence": {"zero_depth_bins": list(range(event_bins[0], event_bins[1] + 1)), "left_flank_mapped_bases": depth_bases[event_bins[0] - 1], "right_flank_mapped_bases": depth_bases[event_bins[1] + 1], "spanning_FR_pairs": support_pairs, "junction_reads": junction_reads},
    "precision": {"rounding_bp": BIN, "start_boundary": event_start, "end_boundary": event_end, "one_based_interval_interpretation": [event_start + 1, event_end]},
}
(OUT / "qc.json").write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8", newline="\n")

report = f"""# chr22 large-deletion call (T0)

## Result

The independent T0 analysis calls an approximately **{event_size:,} bp deletion** between 100-kb-rounded chr22 boundaries **{event_start:,}** and **{event_end:,}**.

## Evidence

- Depth segmentation found {event_bins[1]-event_bins[0]+1} consecutive callable 100-kb bins with no uniquely mapped read bases. The immediately flanking bins contain {depth_bases[event_bins[0]-1]:,} and {depth_bases[event_bins[1]+1]:,} mapped bases.
- {len(support_pairs)} correctly oriented FR pairs bridge the entire interval. Subtracting the deleted reference span restores ordinary library spans; the normal FR median is {median_span:.0f} bp.
- {len(junction_reads)} otherwise-unmapped reads match exactly across the inferred left/right reference join.

The selected chromosome-analysis skill motivated explicit assembly/chromosome validation; genome-annotation guidance motivated unambiguous coordinate reporting; code-execution guidance motivated the self-contained, deterministic standard-library implementation.

## Precision

The required values are rounded to 100 kb. The zero-depth segmentation cannot independently support finer precision, even though exact junction reads are compatible with the same displayed boundaries. Under a 1-based interval convention, the removed bases are approximately {event_start+1:,} through {event_end:,}.
"""
(OUT / "report.md").write_text(report, encoding="utf-8", newline="\n")
print(json.dumps({"chrom": "chr22", "start_100kb": event_start, "end_100kb": event_end, "size_bp": event_size, "spanning_pairs": len(support_pairs), "junction_reads": len(junction_reads)}, sort_keys=True))
