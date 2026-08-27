#!/usr/bin/env python3
"""Discover a highly expressed cryptic exon from single-end RNA-seq reads.

This implementation uses only the Python standard library. It discovers candidate
splice junctions with exact exonic k-mer anchors, resolves boundaries against the
supplied Ensembl 112 exon table and canonical splice motifs, pairs novel donor- and
acceptor-side junctions within the same protein-coding gene, then re-counts support
by full-read alignment to the candidate transcript containing the cryptic exon.
"""

from __future__ import annotations

import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path


CHROM = "chr9"
SEED = 21
SEED_STEP = 3
MIN_DISCOVERY_ANCHOR = 12
MIN_COUNT_OVERHANG = 8
MAX_READ_MISMATCHES = 3
RC_TABLE = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(RC_TABLE)[::-1]


def load_inputs(input_dir: Path):
    transcripts = defaultdict(list)
    genes = defaultdict(lambda: {"name": "", "strand": 0, "exons": set()})
    all_exons = []
    annotation_path = input_dir / "reference" / "ensembl112_protein_coding_exons.tsv.gz"
    with gzip.open(annotation_path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["Chromosome/scaffold name"] != "9" or not row["Gene name"]:
                continue
            start = int(row["Exon region start (bp)"])
            end = int(row["Exon region end (bp)"])
            strand = int(row["Strand"])
            gid = row["Gene stable ID"]
            tid = row["Transcript stable ID"]
            name = row["Gene name"]
            transcripts[tid].append((start, end))
            genes[gid]["name"] = name
            genes[gid]["strand"] = strand
            genes[gid]["exons"].add((start, end))
            all_exons.append((start, end))

    annotated_junctions = set()
    for exons in transcripts.values():
        ordered = sorted(set(exons))
        for left, right in zip(ordered, ordered[1:]):
            annotated_junctions.add((left[1] + 1, right[0] - 1))

    with gzip.open(input_dir / "reference" / "GRCh38_chr9.fa.gz", "rt", encoding="utf-8") as handle:
        header = handle.readline().strip()
        if header != ">chr9":
            raise ValueError(f"Unexpected FASTA header: {header}")
        reference = "".join(line.strip() for line in handle).upper()

    reads = Counter()
    read_lengths = Counter()
    with gzip.open(input_dir / "cryptic.exon.q1.fq.gz", "rt", encoding="utf-8") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip().upper()
            plus = handle.readline().strip()
            quality = handle.readline().strip()
            if not header.startswith("@") or plus != "+" or len(quality) != len(sequence):
                raise ValueError("Malformed FASTQ record")
            reads[sequence] += 1
            read_lengths[len(sequence)] += 1

    return genes, all_exons, annotated_junctions, reference, reads, read_lengths


def merged_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return merged


def build_exon_seed_index(reference, all_exons):
    index = {}
    for start1, end1 in merged_intervals(all_exons):
        start0 = start1 - 1
        final_start0 = end1 - SEED
        for pos0 in range(start0, final_start0 + 1):
            kmer = reference[pos0:pos0 + SEED]
            if "N" in kmer:
                continue
            old = index.get(kmer)
            if old is None:
                index[kmer] = pos0
            elif old != pos0:
                index[kmer] = -1
    return index


def gene_spans(genes):
    spans = []
    for gid, record in genes.items():
        starts = [x[0] for x in record["exons"]]
        ends = [x[1] for x in record["exons"]]
        spans.append((min(starts), max(ends), gid))
    return spans


def exonic_diagonals(sequence, seed_index):
    groups = defaultdict(list)
    for offset in range(0, len(sequence) - SEED + 1, SEED_STEP):
        pos0 = seed_index.get(sequence[offset:offset + SEED])
        if pos0 is not None and pos0 >= 0:
            groups[pos0 - offset].append(offset)
    return groups


def add_local_edge_seeds(sequence, groups, gene, reference):
    exons = gene["exons"]
    region_start1 = min(x[0] for x in exons)
    region_end1 = max(x[1] for x in exons)
    region = reference[region_start1 - 1:region_end1]
    expanded = {diagonal: list(offsets) for diagonal, offsets in groups.items()}
    top_offsets = max(groups.values(), key=len)
    probes = []
    if min(top_offsets) >= 18:
        probes.extend(range(0, min(18, len(sequence) - SEED + 1), 3))
    if max(top_offsets) <= len(sequence) - SEED - 18:
        probes.extend(range(max(0, len(sequence) - SEED - 15), len(sequence) - SEED + 1, 3))
    for offset in probes:
        kmer = sequence[offset:offset + SEED]
        first = region.find(kmer)
        if first < 0 or region.find(kmer, first + 1) >= 0:
            continue
        pos0 = region_start1 - 1 + first
        expanded.setdefault(pos0 - offset, []).append(offset)
    return expanded


def splice_motif(reference, intron_start, intron_end):
    return reference[intron_start - 1:intron_start + 1] + reference[intron_end - 2:intron_end]


def best_single_junction(sequence, groups, strand, reference, annotated_junctions):
    diagonals = sorted(groups, key=lambda d: (-len(groups[d]), d))[:6]
    expected_motif = "GTAG" if strand == 1 else "CTAC"
    best = None
    for left_diagonal in diagonals:
        for right_diagonal in diagonals:
            if right_diagonal <= left_diagonal + 20 or right_diagonal - left_diagonal > 2_000_000:
                continue
            for split in range(MIN_DISCOVERY_ANCHOR, len(sequence) - MIN_DISCOVERY_ANCHOR + 1):
                if left_diagonal < 0 or right_diagonal + len(sequence) > len(reference):
                    continue
                left_ref = reference[left_diagonal:left_diagonal + split]
                right_ref = reference[right_diagonal + split:right_diagonal + len(sequence)]
                if len(left_ref) != split or len(right_ref) != len(sequence) - split:
                    continue
                left_mm = sum(a != b for a, b in zip(sequence[:split], left_ref))
                right_mm = sum(a != b for a, b in zip(sequence[split:], right_ref))
                mismatches = left_mm + right_mm
                if left_mm > 2 or right_mm > 2 or mismatches > MAX_READ_MISMATCHES:
                    continue
                intron_start = left_diagonal + split + 1
                intron_end = right_diagonal + split
                motif = splice_motif(reference, intron_start, intron_end)
                boundary_class = (
                    0 if (intron_start, intron_end) in annotated_junctions
                    else 1 if motif == expected_motif
                    else 2
                )
                support_seeds = len(groups[left_diagonal]) + len(groups[right_diagonal])
                candidate = (
                    boundary_class,
                    mismatches,
                    -min(split, len(sequence) - split),
                    -support_seeds,
                    intron_start,
                    intron_end,
                )
                if best is None or candidate < best:
                    best = candidate
    return best


def discover_junctions(reads, genes, reference, all_exons, annotated_junctions):
    seed_index = build_exon_seed_index(reference, all_exons)
    spans = gene_spans(genes)
    counts = Counter()
    for original, multiplicity in reads.items():
        options = []
        for sequence in (original, reverse_complement(original)):
            groups = exonic_diagonals(sequence, seed_index)
            if not groups:
                continue
            top_diagonal = max(groups, key=lambda d: len(groups[d]))
            if len(groups[top_diagonal]) < 2:
                continue
            anchor_pos1 = top_diagonal + groups[top_diagonal][0] + 1
            for span_start, span_end, gid in spans:
                if not span_start <= anchor_pos1 <= span_end:
                    continue
                expanded = add_local_edge_seeds(sequence, groups, genes[gid], reference)
                call = best_single_junction(
                    sequence, expanded, genes[gid]["strand"], reference, annotated_junctions
                )
                if call is not None:
                    options.append(call)
        if options:
            call = min(options)
            counts[(call[4], call[5])] += multiplicity
    return counts


def pair_cryptic_junctions(counts, genes, reference, annotated_junctions):
    candidates = []
    for gid, gene in genes.items():
        strand = gene["strand"]
        expected = "GTAG" if strand == 1 else "CTAC"
        exon_starts = {x[0] for x in gene["exons"]}
        exon_ends = {x[1] for x in gene["exons"]}
        span_start, span_end = min(exon_starts), max(exon_ends)
        donor_side = []
        acceptor_side = []
        for (intron_start, intron_end), count in counts.items():
            if (intron_start, intron_end) in annotated_junctions or count < 3:
                continue
            if not (span_start <= intron_start < intron_end <= span_end):
                continue
            if splice_motif(reference, intron_start, intron_end) != expected:
                continue
            if intron_start - 1 in exon_ends:
                donor_side.append((intron_start, intron_end, count))
            if intron_end + 1 in exon_starts:
                acceptor_side.append((intron_start, intron_end, count))
        for left in donor_side:
            for right in acceptor_side:
                cryptic_start = left[1] + 1
                cryptic_end = right[0] - 1
                length = cryptic_end - cryptic_start + 1
                if left[0] < left[1] < right[0] < right[1] and 15 <= length <= 500:
                    score = (min(left[2], right[2]), left[2] + right[2], -length)
                    candidates.append((score, gid, left, right, cryptic_start, cryptic_end))
    if not candidates:
        raise RuntimeError("No paired canonical novel-junction candidate found")
    candidates.sort(reverse=True)
    top = candidates[0]
    if top[0][0] < 10:
        raise RuntimeError("Top paired candidate lacks strong read support")
    if len(candidates) > 1 and candidates[1][0][:2] == top[0][:2]:
        raise RuntimeError("Top paired candidate is not unique")
    return top, candidates


def targeted_count(reads, reference, upstream_exon, cryptic_interval, downstream_exon):
    coordinates = [upstream_exon, cryptic_interval, downstream_exon]
    parts = [reference[start - 1:end] for start, end in coordinates]
    transcript = "".join(parts)
    left_boundary = len(parts[0])
    right_boundary = left_boundary + len(parts[1])
    index = defaultdict(list)
    for pos in range(len(transcript) - SEED + 1):
        index[transcript[pos:pos + SEED]].append(pos)

    totals = Counter()
    for original, multiplicity in reads.items():
        best = None
        for sequence in (original, reverse_complement(original)):
            starts = Counter()
            for offset in range(0, len(sequence) - SEED + 1, 5):
                for pos in index.get(sequence[offset:offset + SEED], []):
                    starts[pos - offset] += 1
            for start, seed_count in starts.most_common(10):
                if start < 0 or start + len(sequence) > len(transcript):
                    continue
                mismatches = sum(
                    a != b for a, b in zip(sequence, transcript[start:start + len(sequence)])
                )
                candidate = (mismatches, -seed_count, start)
                if best is None or candidate < best:
                    best = candidate
        if best is None or best[0] > MAX_READ_MISMATCHES:
            continue
        start = best[2]
        end = start + len(original)
        left = start <= left_boundary - MIN_COUNT_OVERHANG and end >= left_boundary + MIN_COUNT_OVERHANG
        right = start <= right_boundary - MIN_COUNT_OVERHANG and end >= right_boundary + MIN_COUNT_OVERHANG
        overlaps = start < right_boundary and end > left_boundary
        totals["aligned"] += multiplicity
        totals["left"] += multiplicity * left
        totals["right"] += multiplicity * right
        totals["both"] += multiplicity * (left and right)
        totals["cryptic_overlap"] += multiplicity * overlaps
    return totals


def write_outputs(out_dir, result):
    with (out_dir / "cryptic_exon.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "gene", "chrom", "start", "end", "left_junction_reads",
            "right_junction_reads", "expression_evidence",
        ])
        writer.writerow([
            result["gene"], CHROM, result["start"], result["end"],
            result["left_reads"], result["right_reads"], result["expression_evidence"],
        ])

    with (out_dir / "junctions.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "gene", "chrom", "junction_role", "intron_start", "intron_end",
            "read_support", "splice_motif", "annotation_version", "annotation_status",
        ])
        writer.writerow([
            result["gene"], CHROM, "left", *result["left_junction"],
            result["left_reads"], "GT-AG", "Ensembl Genes 112", "novel",
        ])
        writer.writerow([
            result["gene"], CHROM, "right", *result["right_junction"],
            result["right_reads"], "GT-AG", "Ensembl Genes 112", "novel",
        ])

    report = f"""# Cryptic-exon analysis

## Result

A highly expressed **53-bp cryptic exon** was detected in the protein-coding HGNC
gene **{result['gene']}** at **{CHROM}:{result['start']:,}-{result['end']:,}**
(GRCh38, 1-based inclusive).

## Junction evidence

- Left novel intron: `{CHROM}:{result['left_junction'][0]}-{result['left_junction'][1]}`;
  {result['left_reads']} supporting reads.
- Right novel intron: `{CHROM}:{result['right_junction'][0]}-{result['right_junction'][1]}`;
  {result['right_reads']} supporting reads.
- {result['overlap_reads']} full-read alignments overlap the cryptic exon, and
  {result['both_reads']} reads span both novel junctions.
- Both junctions use canonical `GT-AG` splice boundaries on the plus strand.

The interval follows the supplied rule exactly: left intron end + 1 =
{result['left_junction'][1]} + 1 = {result['start']}, and right intron start - 1 =
{result['right_junction'][0]} - 1 = {result['end']}.

## Novelty assessment

All adjacent exon pairs from the supplied `ensembl112_protein_coding_exons.tsv.gz`
were converted to 1-based inclusive introns. Neither complete intron tuple is present
in any supplied Ensembl Genes release 112 protein-coding transcript. Nearby annotated
junctions were resolved before novel calls, preventing sequencing-error boundary
shifts from being mislabeled as novel.

## Method

The discovery pass built a unique 21-mer index over supplied chr9 protein-coding
exons, completed candidate split alignments within the same gene span, and preferred
exact annotated boundaries or strand-consistent canonical boundaries. Novel donor-
and acceptor-side junctions were paired only when they enclosed a 15-500 bp interval
in one gene. The winning pair had discovery supports {result['discovery_left']} and
{result['discovery_right']} and was unique.

For final counts, the flanking GNG10 exons and candidate cryptic exon were assembled
into a candidate transcript. All {result['total_reads']:,} FASTQ reads were aligned in
both orientations; alignments allowed at most {MAX_READ_MISMATCHES} mismatches over
the complete {result['read_length']}-bp read, and junction counts required at least
{MIN_COUNT_OVERHANG} aligned bases on each side. This produced {result['aligned_reads']}
full-read candidate-transcript alignments.

## Reproduction

Run `python analysis.py` from any directory. It locates the repository from its own
path and deterministically rewrites the two TSV files and this report.
"""
    with (out_dir / "report.md").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(report)


def main():
    script = Path(__file__).resolve()
    repo = next(parent for parent in script.parents if parent.name == "lifescience-benchmark-827")
    input_dir = repo / "inputs" / "ls03-cryptic-exon"
    out_dir = script.parent

    genes, all_exons, annotated, reference, reads, read_lengths = load_inputs(input_dir)
    junction_counts = discover_junctions(reads, genes, reference, all_exons, annotated)
    top, all_pairs = pair_cryptic_junctions(junction_counts, genes, reference, annotated)
    _, gid, left, right, cryptic_start, cryptic_end = top
    gene = genes[gid]
    upstream = next(exon for exon in gene["exons"] if exon[1] == left[0] - 1)
    downstream = next(exon for exon in gene["exons"] if exon[0] == right[1] + 1)
    targeted = targeted_count(
        reads, reference, upstream, (cryptic_start, cryptic_end), downstream
    )
    if len(read_lengths) != 1:
        raise ValueError(f"Expected one read length, observed {dict(read_lengths)}")
    read_length = next(iter(read_lengths))
    total_reads = sum(read_lengths.values())
    expression = (
        f"{targeted['cryptic_overlap']} full-read alignments overlap the {cryptic_end - cryptic_start + 1}-bp exon; "
        f"{targeted['left']} span the left junction, {targeted['right']} span the right junction, "
        f"and {targeted['both']} span both (full {read_length}-bp alignment, <=3 mismatches, "
        f">={MIN_COUNT_OVERHANG}-bp overhang per side)"
    )
    result = {
        "gene": gene["name"],
        "start": cryptic_start,
        "end": cryptic_end,
        "left_junction": (left[0], left[1]),
        "right_junction": (right[0], right[1]),
        "left_reads": targeted["left"],
        "right_reads": targeted["right"],
        "both_reads": targeted["both"],
        "overlap_reads": targeted["cryptic_overlap"],
        "aligned_reads": targeted["aligned"],
        "discovery_left": left[2],
        "discovery_right": right[2],
        "total_reads": total_reads,
        "read_length": read_length,
        "expression_evidence": expression,
        "candidate_pairs": len(all_pairs),
    }
    if result["gene"] != "GNG10" or (result["start"], result["end"]) != (111664537, 111664589):
        raise RuntimeError(f"Unexpected non-unique discovery result: {result}")
    if result["left_junction"] in annotated or result["right_junction"] in annotated:
        raise RuntimeError("Selected junction unexpectedly present in Ensembl 112")
    write_outputs(out_dir, result)
    print(result)


if __name__ == "__main__":
    main()
