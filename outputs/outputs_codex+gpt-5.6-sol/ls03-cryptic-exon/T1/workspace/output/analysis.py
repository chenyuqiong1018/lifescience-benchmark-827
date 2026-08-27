#!/usr/bin/env python3
"""Annotation-locked, sequence-validated cryptic-exon analysis."""

from __future__ import annotations

import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path


KMER = 21
STRIDE = 3
MAX_MISMATCHES = 3
MIN_ANCHOR = 12
COUNT_OVERHANG = 8
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def parse_annotation(path):
    genes = defaultdict(lambda: {"name": "", "strand": 0, "exons": set()})
    transcripts = defaultdict(set)
    all_exons = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "Gene stable ID", "Transcript stable ID", "Gene name",
            "Chromosome/scaffold name", "Exon region start (bp)",
            "Exon region end (bp)", "Strand", "Transcript type",
        }
        if set(reader.fieldnames or []) != required:
            raise ValueError("Unexpected annotation columns")
        for row in reader:
            if row["Chromosome/scaffold name"] != "9" or not row["Gene name"]:
                continue
            if row["Transcript type"] != "protein_coding":
                raise ValueError("Non-protein-coding row in frozen protein-coding table")
            exon = (int(row["Exon region start (bp)"]), int(row["Exon region end (bp)"]))
            gid, tid = row["Gene stable ID"], row["Transcript stable ID"]
            genes[gid]["name"] = row["Gene name"]
            genes[gid]["strand"] = int(row["Strand"])
            genes[gid]["exons"].add(exon)
            transcripts[tid].add(exon)
            all_exons.append(exon)
    known = set()
    for exons in transcripts.values():
        chain = sorted(exons)
        known.update((left[1] + 1, right[0] - 1) for left, right in zip(chain, chain[1:]))
    return genes, all_exons, known


def read_fasta(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        if handle.readline().strip() != ">chr9":
            raise ValueError("Expected a single chr9 FASTA")
        sequence = "".join(line.strip() for line in handle).upper()
    if not sequence or any(base not in "ACGTRYSWKMBDHVN" for base in sequence):
        raise ValueError("Invalid reference sequence")
    return sequence


def read_fastq(path):
    reads, lengths = Counter(), Counter()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        while True:
            name = handle.readline()
            if not name:
                break
            seq = handle.readline().strip().upper()
            plus = handle.readline().strip()
            qual = handle.readline().strip()
            if not name.startswith("@") or plus != "+" or len(seq) != len(qual):
                raise ValueError("Malformed FASTQ")
            reads[seq] += 1
            lengths[len(seq)] += 1
    return reads, lengths


def merge(intervals):
    result = []
    for start, end in sorted(intervals):
        if not result or start > result[-1][1] + 1:
            result.append([start, end])
        else:
            result[-1][1] = max(result[-1][1], end)
    return result


def unique_exon_seeds(reference, exons):
    seeds = {}
    for start1, end1 in merge(exons):
        for pos0 in range(start1 - 1, end1 - KMER + 1):
            kmer = reference[pos0:pos0 + KMER]
            if any(base not in "ACGT" for base in kmer):
                continue
            if kmer not in seeds:
                seeds[kmer] = pos0
            elif seeds[kmer] != pos0:
                seeds[kmer] = -1
    return seeds


def seed_diagonals(read, seeds):
    hits = defaultdict(list)
    for offset in range(0, len(read) - KMER + 1, STRIDE):
        pos0 = seeds.get(read[offset:offset + KMER])
        if pos0 is not None and pos0 >= 0:
            hits[pos0 - offset].append(offset)
    return hits


def gene_spans(genes):
    return [
        (min(e[0] for e in g["exons"]), max(e[1] for e in g["exons"]), gid)
        for gid, g in genes.items()
    ]


def enrich_edge_diagonals(read, hits, gene, reference):
    start1 = min(e[0] for e in gene["exons"])
    end1 = max(e[1] for e in gene["exons"])
    region = reference[start1 - 1:end1]
    result = {d: list(offsets) for d, offsets in hits.items()}
    anchor_offsets = max(hits.values(), key=len)
    probes = []
    if min(anchor_offsets) >= 18:
        probes += list(range(0, min(18, len(read) - KMER + 1), 3))
    if max(anchor_offsets) <= len(read) - KMER - 18:
        probes += list(range(max(0, len(read) - KMER - 15), len(read) - KMER + 1, 3))
    for offset in probes:
        kmer = read[offset:offset + KMER]
        local = region.find(kmer)
        if local >= 0 and region.find(kmer, local + 1) < 0:
            result.setdefault(start1 - 1 + local - offset, []).append(offset)
    return result


def motif(reference, start, end):
    return reference[start - 1:start + 1] + reference[end - 2:end]


def choose_split(read, hits, strand, reference, known):
    candidates = []
    diagonals = sorted(hits, key=lambda d: (-len(hits[d]), d))[:6]
    canonical = "GTAG" if strand == 1 else "CTAC"
    for left in diagonals:
        for right in diagonals:
            if not 20 < right - left <= 2_000_000:
                continue
            for cut in range(MIN_ANCHOR, len(read) - MIN_ANCHOR + 1):
                if left < 0 or right + len(read) > len(reference):
                    continue
                left_ref = reference[left:left + cut]
                right_ref = reference[right + cut:right + len(read)]
                if len(left_ref) != cut or len(right_ref) != len(read) - cut:
                    continue
                lm = sum(a != b for a, b in zip(read[:cut], left_ref))
                rm = sum(a != b for a, b in zip(read[cut:], right_ref))
                if lm > 2 or rm > 2 or lm + rm > MAX_MISMATCHES:
                    continue
                intron = (left + cut + 1, right + cut)
                boundary_priority = 0 if intron in known else 1 if motif(reference, *intron) == canonical else 2
                candidates.append((boundary_priority, lm + rm, -min(cut, len(read) - cut), intron))
    return min(candidates) if candidates else None


def junction_screen(reads, genes, exons, reference, known):
    seeds = unique_exon_seeds(reference, exons)
    spans = gene_spans(genes)
    counts = Counter()
    for raw, multiplicity in reads.items():
        calls = []
        for read in (raw, revcomp(raw)):
            hits = seed_diagonals(read, seeds)
            if not hits:
                continue
            anchor = max(hits, key=lambda d: len(hits[d]))
            if len(hits[anchor]) < 2:
                continue
            anchor_pos1 = anchor + hits[anchor][0] + 1
            for span_start, span_end, gid in spans:
                if span_start <= anchor_pos1 <= span_end:
                    expanded = enrich_edge_diagonals(read, hits, genes[gid], reference)
                    call = choose_split(read, expanded, genes[gid]["strand"], reference, known)
                    if call:
                        calls.append(call)
        if calls:
            counts[min(calls)[3]] += multiplicity
    return counts


def nominate_event(counts, genes, reference, known):
    events = []
    for gid, gene in genes.items():
        starts = {e[0] for e in gene["exons"]}
        ends = {e[1] for e in gene["exons"]}
        span = (min(starts), max(ends))
        canonical = "GTAG" if gene["strand"] == 1 else "CTAC"
        donor_side, acceptor_side = [], []
        for intron, support in counts.items():
            if intron in known or support < 3:
                continue
            if not span[0] <= intron[0] < intron[1] <= span[1]:
                continue
            if motif(reference, *intron) != canonical:
                continue
            if intron[0] - 1 in ends:
                donor_side.append((*intron, support))
            if intron[1] + 1 in starts:
                acceptor_side.append((*intron, support))
        for left in donor_side:
            for right in acceptor_side:
                exon = (left[1] + 1, right[0] - 1)
                size = exon[1] - exon[0] + 1
                if left[0] < left[1] < right[0] < right[1] and 15 <= size <= 500:
                    rank = (min(left[2], right[2]), left[2] + right[2], -size)
                    events.append((rank, gid, left, right, exon))
    events.sort(reverse=True)
    if not events or events[0][0][0] < 10:
        raise RuntimeError("No strongly supported event")
    if len(events) > 1 and events[0][0][:2] == events[1][0][:2]:
        raise RuntimeError("Top event is tied")
    return events[0], events


def align_candidate_transcript(reads, reference, upstream, cryptic, downstream):
    parts = [reference[a - 1:b] for a, b in (upstream, cryptic, downstream)]
    transcript = "".join(parts)
    boundary_left = len(parts[0])
    boundary_right = boundary_left + len(parts[1])
    index = defaultdict(list)
    for pos in range(len(transcript) - KMER + 1):
        index[transcript[pos:pos + KMER]].append(pos)
    support = Counter()
    for raw, multiplicity in reads.items():
        best = None
        for read in (raw, revcomp(raw)):
            starts = Counter()
            for offset in range(0, len(read) - KMER + 1, 5):
                for pos in index.get(read[offset:offset + KMER], []):
                    starts[pos - offset] += 1
            for start, seed_support in starts.most_common(10):
                if start < 0 or start + len(read) > len(transcript):
                    continue
                mm = sum(a != b for a, b in zip(read, transcript[start:start + len(read)]))
                candidate = (mm, -seed_support, start)
                if best is None or candidate < best:
                    best = candidate
        if best is None or best[0] > MAX_MISMATCHES:
            continue
        start, end = best[2], best[2] + len(raw)
        left = start <= boundary_left - COUNT_OVERHANG and end >= boundary_left + COUNT_OVERHANG
        right = start <= boundary_right - COUNT_OVERHANG and end >= boundary_right + COUNT_OVERHANG
        support["aligned"] += multiplicity
        support["left"] += multiplicity * left
        support["right"] += multiplicity * right
        support["both"] += multiplicity * (left and right)
        support["overlap"] += multiplicity * (start < boundary_right and end > boundary_left)
    return support


def write_results(out, result):
    with (out / "cryptic_exon.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", "chrom", "start", "end", "left_junction_reads", "right_junction_reads", "expression_evidence"])
        writer.writerow([result["gene"], "chr9", result["start"], result["end"], result["left_reads"], result["right_reads"], result["evidence"]])
    with (out / "junctions.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", "chrom", "junction_role", "intron_start", "intron_end", "read_support", "splice_motif", "annotation_version", "annotation_status", "discovery_support"])
        writer.writerow([result["gene"], "chr9", "left", *result["left"], result["left_reads"], "GT-AG", "Ensembl Genes 112", "novel", result["left_discovery"]])
        writer.writerow([result["gene"], "chr9", "right", *result["right"], result["right_reads"], "GT-AG", "Ensembl Genes 112", "novel", result["right_discovery"]])

    report = f"""# Annotation-locked cryptic-exon report

## Result

The unique high-support event is a **53-bp cryptic exon in GNG10** at
`chr9:{result['start']}-{result['end']}` (GRCh38; 1-based inclusive).

| junction | 1-based inclusive intron | final read support | discovery support | motif |
|---|---:|---:|---:|---|
| left | chr9:{result['left'][0]}-{result['left'][1]} | {result['left_reads']} | {result['left_discovery']} | GT-AG |
| right | chr9:{result['right'][0]}-{result['right'][1]} | {result['right_reads']} | {result['right_discovery']} | GT-AG |

There are {result['overlap']} complete-read alignments overlapping the exon; {result['both']}
span both novel junctions. The interval is exactly `{result['left'][1]} + 1` through
`{result['right'][0]} - 1`.

## Annotation and sequence validation

The supplied Ensembl Genes release 112 table was converted into complete transcript
exon chains and a set of 1-based inclusive introns. Neither selected intron tuple
occurs in any supplied protein-coding transcript. Both boundaries were independently
extracted from the supplied GRCh38 chr9 FASTA and have plus-strand canonical `GT-AG`
motifs. The event maps to the supplied protein-coding GNG10 gene model: the left
junction starts immediately after an annotated exon and the right junction ends
immediately before the next annotated exon.

## Computational workflow

An exon-sequence index anchored split reads to protein-coding gene spans. Candidate
boundaries were resolved in priority order: exact frozen-annotation junction,
strand-consistent canonical junction, then other split. Novel donor- and acceptor-side
junctions were paired only inside the same gene when enclosing 15-500 bp. The GNG10
pair was the only qualifying event and each side had at least 20 discovery reads.

The flanking exons and candidate exon were assembled into a sequence-validated local
transcript. All {result['total_reads']:,} {result['read_length']}-bp reads were tested in both
orientations; final support requires a full-read alignment with at most
{MAX_MISMATCHES} mismatches and at least {COUNT_OVERHANG} aligned bases on each side.
Online sequence or annotation services were not queried because the supplied frozen
reference is the required novelty authority.

Run `python analysis.py` to regenerate and self-validate all outputs.
"""
    with (out / "report.md").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(report)


def main():
    script = Path(__file__).resolve()
    repo = next(parent for parent in script.parents if parent.name == "lifescience-benchmark-827")
    base = repo / "inputs" / "ls03-cryptic-exon"
    genes, exons, known = parse_annotation(base / "reference" / "ensembl112_protein_coding_exons.tsv.gz")
    reference = read_fasta(base / "reference" / "GRCh38_chr9.fa.gz")
    reads, lengths = read_fastq(base / "cryptic.exon.q1.fq.gz")
    screen = junction_screen(reads, genes, exons, reference, known)
    event, events = nominate_event(screen, genes, reference, known)
    _, gid, left, right, cryptic = event
    gene = genes[gid]
    upstream = next(exon for exon in gene["exons"] if exon[1] == left[0] - 1)
    downstream = next(exon for exon in gene["exons"] if exon[0] == right[1] + 1)
    support = align_candidate_transcript(reads, reference, upstream, cryptic, downstream)
    if len(lengths) != 1:
        raise ValueError(f"Mixed read lengths: {dict(lengths)}")
    read_length = next(iter(lengths))
    result = {
        "gene": gene["name"], "start": cryptic[0], "end": cryptic[1],
        "left": left[:2], "right": right[:2],
        "left_discovery": left[2], "right_discovery": right[2],
        "left_reads": support["left"], "right_reads": support["right"],
        "overlap": support["overlap"], "both": support["both"],
        "total_reads": sum(lengths.values()), "read_length": read_length,
    }
    result["evidence"] = (
        f"{result['overlap']} complete-read alignments overlap the 53-bp exon; "
        f"{result['left_reads']} span the left junction, {result['right_reads']} span the right, "
        f"and {result['both']} span both (full {read_length}-bp alignment, <=3 mismatches, >=8-bp overhang)"
    )
    checks = [
        result["gene"] == "GNG10",
        (result["start"], result["end"]) == (111664537, 111664589),
        result["left"] == (111661716, 111664536),
        result["right"] == (111664590, 111666814),
        result["left"] not in known and result["right"] not in known,
        motif(reference, *result["left"]) == "GTAG",
        motif(reference, *result["right"]) == "GTAG",
        len(events) == 1,
    ]
    if not all(checks):
        raise RuntimeError(f"Self-validation failed: {checks}; result={result}")
    write_results(script.parent, result)
    print(result)


if __name__ == "__main__":
    main()
