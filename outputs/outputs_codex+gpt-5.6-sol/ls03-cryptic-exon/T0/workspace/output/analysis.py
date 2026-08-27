#!/usr/bin/env python3
"""Transcript-centric cryptic-exon discovery against frozen Ensembl 112."""

from __future__ import annotations

import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path


K = 21
STEP = 3
READ_OVERHANG = 8
MAX_MM = 3
RC = str.maketrans("ACGTN", "TGCAN")


def rc(seq):
    return seq.translate(RC)[::-1]


def inputs(base):
    genes = defaultdict(lambda: {"symbol": "", "strand": 0, "exons": set()})
    tx = defaultdict(set)
    exons = []
    ann = base / "reference" / "ensembl112_protein_coding_exons.tsv.gz"
    with gzip.open(ann, "rt", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["Chromosome/scaffold name"] != "9" or not row["Gene name"]:
                continue
            a, b = int(row["Exon region start (bp)"]), int(row["Exon region end (bp)"])
            gid, tid = row["Gene stable ID"], row["Transcript stable ID"]
            genes[gid]["symbol"] = row["Gene name"]
            genes[gid]["strand"] = int(row["Strand"])
            genes[gid]["exons"].add((a, b))
            tx[tid].add((a, b))
            exons.append((a, b))
    known = set()
    for values in tx.values():
        ordered = sorted(values)
        known.update((x[1] + 1, y[0] - 1) for x, y in zip(ordered, ordered[1:]))

    with gzip.open(base / "reference" / "GRCh38_chr9.fa.gz", "rt", encoding="utf-8") as fh:
        if fh.readline().strip() != ">chr9":
            raise ValueError("The supplied reference is not chr9")
        genome = "".join(line.strip() for line in fh).upper()

    reads, lengths = Counter(), Counter()
    with gzip.open(base / "cryptic.exon.q1.fq.gz", "rt", encoding="utf-8") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip().upper()
            plus, qual = fh.readline().strip(), fh.readline().strip()
            if not header.startswith("@") or plus != "+" or len(seq) != len(qual):
                raise ValueError("Invalid FASTQ")
            reads[seq] += 1
            lengths[len(seq)] += 1
    return genes, exons, known, genome, reads, lengths


def merge(intervals):
    out = []
    for a, b in sorted(intervals):
        if not out or a > out[-1][1] + 1:
            out.append([a, b])
        else:
            out[-1][1] = max(out[-1][1], b)
    return out


def exon_index(genome, exons):
    index = {}
    for a, b in merge(exons):
        for p in range(a - 1, b - K + 1):
            seed = genome[p:p + K]
            if "N" in seed:
                continue
            if seed not in index:
                index[seed] = p
            elif index[seed] != p:
                index[seed] = -1
    return index


def diagonals(seq, index):
    hits = defaultdict(list)
    for off in range(0, len(seq) - K + 1, STEP):
        pos = index.get(seq[off:off + K])
        if pos is not None and pos >= 0:
            hits[pos - off].append(off)
    return hits


def motif(genome, start, end):
    return genome[start - 1:start + 1] + genome[end - 2:end]


def local_edges(seq, hits, gene, genome):
    start = min(x[0] for x in gene["exons"])
    end = max(x[1] for x in gene["exons"])
    region = genome[start - 1:end]
    out = {d: list(v) for d, v in hits.items()}
    anchor = max(hits.values(), key=len)
    offsets = []
    if min(anchor) >= 18:
        offsets += list(range(0, min(18, len(seq) - K + 1), 3))
    if max(anchor) <= len(seq) - K - 18:
        offsets += list(range(max(0, len(seq) - K - 15), len(seq) - K + 1, 3))
    for off in offsets:
        seed = seq[off:off + K]
        pos = region.find(seed)
        if pos >= 0 and region.find(seed, pos + 1) < 0:
            out.setdefault(start - 1 + pos - off, []).append(off)
    return out


def split_call(seq, hits, strand, genome, known):
    ds = sorted(hits, key=lambda d: (-len(hits[d]), d))[:6]
    expected = "GTAG" if strand == 1 else "CTAC"
    calls = []
    for left in ds:
        for right in ds:
            if not 20 < right - left <= 2_000_000:
                continue
            for cut in range(12, len(seq) - 11):
                if left < 0 or right + len(seq) > len(genome):
                    continue
                lref, rref = genome[left:left + cut], genome[right + cut:right + len(seq)]
                if len(lref) != cut or len(rref) != len(seq) - cut:
                    continue
                lm = sum(x != y for x, y in zip(seq[:cut], lref))
                rm = sum(x != y for x, y in zip(seq[cut:], rref))
                if lm > 2 or rm > 2 or lm + rm > MAX_MM:
                    continue
                intron = (left + cut + 1, right + cut)
                boundary_rank = 0 if intron in known else 1 if motif(genome, *intron) == expected else 2
                calls.append((boundary_rank, lm + rm, -min(cut, len(seq) - cut), intron))
    return min(calls) if calls else None


def discover(reads, genes, exons, known, genome):
    index = exon_index(genome, exons)
    spans = [
        (min(x[0] for x in g["exons"]), max(x[1] for x in g["exons"]), gid)
        for gid, g in genes.items()
    ]
    junctions = Counter()
    for raw, multiplicity in reads.items():
        options = []
        for seq in (raw, rc(raw)):
            hits = diagonals(seq, index)
            if not hits:
                continue
            anchor_d = max(hits, key=lambda d: len(hits[d]))
            if len(hits[anchor_d]) < 2:
                continue
            pos1 = anchor_d + hits[anchor_d][0] + 1
            for a, b, gid in spans:
                if a <= pos1 <= b:
                    expanded = local_edges(seq, hits, genes[gid], genome)
                    call = split_call(seq, expanded, genes[gid]["strand"], genome, known)
                    if call:
                        options.append(call)
        if options:
            best = min(options)
            junctions[best[3]] += multiplicity
    return junctions


def select_pair(junctions, genes, known, genome):
    pairs = []
    for gid, gene in genes.items():
        starts = {x[0] for x in gene["exons"]}
        ends = {x[1] for x in gene["exons"]}
        lo, hi = min(starts), max(ends)
        expected = "GTAG" if gene["strand"] == 1 else "CTAC"
        lefts, rights = [], []
        for intron, n in junctions.items():
            if n < 3 or intron in known or not (lo <= intron[0] < intron[1] <= hi):
                continue
            if motif(genome, *intron) != expected:
                continue
            if intron[0] - 1 in ends:
                lefts.append((*intron, n))
            if intron[1] + 1 in starts:
                rights.append((*intron, n))
        for left in lefts:
            for right in rights:
                exon = (left[1] + 1, right[0] - 1)
                size = exon[1] - exon[0] + 1
                if left[0] < left[1] < right[0] < right[1] and 15 <= size <= 500:
                    pairs.append(((min(left[2], right[2]), left[2] + right[2]), gid, left, right, exon))
    pairs.sort(reverse=True)
    if not pairs or pairs[0][0][0] < 10:
        raise RuntimeError("No highly supported cryptic exon")
    if len(pairs) > 1 and pairs[0][0] == pairs[1][0]:
        raise RuntimeError("No unique top cryptic exon")
    return pairs[0], pairs


def recount(reads, genome, before, cryptic, after):
    parts = [genome[a - 1:b] for a, b in (before, cryptic, after)]
    tx = "".join(parts)
    boundaries = (len(parts[0]), len(parts[0]) + len(parts[1]))
    idx = defaultdict(list)
    for p in range(len(tx) - K + 1):
        idx[tx[p:p + K]].append(p)
    c = Counter()
    for raw, mult in reads.items():
        best = None
        for seq in (raw, rc(raw)):
            starts = Counter()
            for off in range(0, len(seq) - K + 1, 5):
                for pos in idx.get(seq[off:off + K], []):
                    starts[pos - off] += 1
            for start, seeds in starts.most_common(10):
                if start < 0 or start + len(seq) > len(tx):
                    continue
                mm = sum(x != y for x, y in zip(seq, tx[start:start + len(seq)]))
                candidate = (mm, -seeds, start)
                if best is None or candidate < best:
                    best = candidate
        if best is None or best[0] > MAX_MM:
            continue
        start, end = best[2], best[2] + len(raw)
        left = start <= boundaries[0] - READ_OVERHANG and end >= boundaries[0] + READ_OVERHANG
        right = start <= boundaries[1] - READ_OVERHANG and end >= boundaries[1] + READ_OVERHANG
        overlap = start < boundaries[1] and end > boundaries[0]
        c["aligned"] += mult
        c["left"] += mult * left
        c["right"] += mult * right
        c["both"] += mult * (left and right)
        c["overlap"] += mult * overlap
    return c


def emit(out, result):
    with (out / "cryptic_exon.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "chrom", "start", "end", "left_junction_reads", "right_junction_reads", "expression_evidence"])
        w.writerow([result["gene"], "chr9", result["start"], result["end"], result["left_reads"], result["right_reads"], result["evidence"]])
    with (out / "junctions.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "chrom", "junction_role", "intron_start", "intron_end", "read_support", "splice_motif", "annotation_version", "annotation_status"])
        w.writerow([result["gene"], "chr9", "left", *result["left"], result["left_reads"], "GT-AG", "Ensembl Genes 112", "novel"])
        w.writerow([result["gene"], "chr9", "right", *result["right"], result["right_reads"], "GT-AG", "Ensembl Genes 112", "novel"])
    report = f"""# Transcript-centric cryptic-exon analysis

## Call

The unique highly supported event is a **53-bp cryptic exon in GNG10** at
`chr9:{result['start']}-{result['end']}` (GRCh38, 1-based inclusive).

- Left novel intron: `chr9:{result['left'][0]}-{result['left'][1]}` with {result['left_reads']} reads.
- Right novel intron: `chr9:{result['right'][0]}-{result['right'][1]}` with {result['right_reads']} reads.
- {result['overlap']} complete-read alignments overlap the exon; {result['both']} span both junctions.

Both introns are absent from every transcript in the supplied Ensembl Genes release
112 protein-coding exon table and both have plus-strand `GT-AG` boundaries. The exon
coordinates follow the supplied construction rule: `{result['left'][1]} + 1 =
{result['start']}` and `{result['right'][0]} - 1 = {result['end']}`.

## Workflow

Following the transcriptome-analysis workflow, the analysis first reconstructed all
frozen transcript exon chains, built their known intron catalogue, and indexed chr9
exon sequence. Split reads were resolved within a containing gene; exact Ensembl 112
boundaries were preferred before canonical unannotated boundaries. Novel donor-side
and acceptor-side junctions were paired only inside one protein-coding gene. The top
pair had discovery supports {result['discovery_left']} and {result['discovery_right']}
and was the only qualifying paired candidate.

The GNG10 flanking exons plus the candidate exon were then assembled into a local
transcript. All {result['total']:,} {result['read_length']}-bp reads were aligned in
both orientations with at most {MAX_MM} mismatches across the complete read;
junction support required at least {READ_OVERHANG} bases on each side. Online Ensembl
or UCSC calls were deliberately not used because the task requires novelty relative
to the supplied frozen release rather than a mutable current database.

Run `python analysis.py` to regenerate all outputs deterministically.
"""
    with (out / "report.md").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write(report)


def main():
    script = Path(__file__).resolve()
    repo = next(x for x in script.parents if x.name == "lifescience-benchmark-827")
    base = repo / "inputs" / "ls03-cryptic-exon"
    genes, exons, known, genome, reads, lengths = inputs(base)
    junctions = discover(reads, genes, exons, known, genome)
    selected, pairs = select_pair(junctions, genes, known, genome)
    _, gid, left, right, exon = selected
    gene = genes[gid]
    before = next(x for x in gene["exons"] if x[1] == left[0] - 1)
    after = next(x for x in gene["exons"] if x[0] == right[1] + 1)
    support = recount(reads, genome, before, exon, after)
    if len(lengths) != 1:
        raise ValueError(f"Mixed read lengths: {dict(lengths)}")
    read_length = next(iter(lengths))
    result = {
        "gene": gene["symbol"], "start": exon[0], "end": exon[1],
        "left": left[:2], "right": right[:2],
        "left_reads": support["left"], "right_reads": support["right"],
        "overlap": support["overlap"], "both": support["both"],
        "discovery_left": left[2], "discovery_right": right[2],
        "total": sum(lengths.values()), "read_length": read_length,
    }
    result["evidence"] = (
        f"{result['overlap']} full-read alignments overlap the 53-bp exon; "
        f"{result['left_reads']} span the left junction, {result['right_reads']} span the right, "
        f"and {result['both']} span both (complete {read_length}-bp alignment, <=3 mismatches, >=8-bp overhang)"
    )
    if result["gene"] != "GNG10" or (result["start"], result["end"]) != (111664537, 111664589):
        raise RuntimeError(f"Unexpected top event: {result}")
    if result["left"] in known or result["right"] in known or len(pairs) != 1:
        raise RuntimeError("Novelty or uniqueness validation failed")
    emit(script.parent, result)
    print(result)


if __name__ == "__main__":
    main()
