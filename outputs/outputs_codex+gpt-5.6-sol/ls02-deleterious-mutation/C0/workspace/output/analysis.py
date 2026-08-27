from __future__ import annotations

import collections
import csv
import gzip
import json
import re
import statistics
from pathlib import Path


OUT = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[6]
DATA = ROOT / "inputs" / "ls02-deleterious-mutation"
FASTQ = DATA / "deleterious.mutation.q2.R1.fq.gz"
FASTA = DATA / "reference" / "GRCh38_chr9.fa.gz"
GTF = DATA / "reference" / "gencode.v47.chr9.annotation.gtf.gz"
CHROM = "chr9"
POS1 = 127_661_125
POS0 = POS1 - 1
REF = "G"
ALT = "T"
MIN_BASE_Q = 20
K = 31
COMP = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMP)[::-1]


def load_reference() -> str:
    with gzip.open(FASTA, "rt") as handle:
        handle.readline()
        return "".join(line.strip().upper() for line in handle)


def read_fastq():
    with gzip.open(FASTQ, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline().strip().upper()
            handle.readline()
            quality = handle.readline().strip()
            yield sequence, quality


def local_alignment(sequence: str, window: str, window_start: int):
    for is_reverse, oriented in ((False, sequence), (True, reverse_complement(sequence))):
        offsets = sorted(set((0, max(0, len(oriented)//3-K//2), max(0, 2*len(oriented)//3-K//2), max(0, len(oriented)-K))))
        votes = collections.Counter()
        for offset in offsets:
            if offset + K > len(oriented):
                continue
            seed = oriented[offset:offset+K]
            first = window.find(seed)
            if first >= 0 and window.find(seed, first+1) < 0:
                votes[window_start + first - offset] += 1
        if not votes:
            continue
        start, support = votes.most_common(1)[0]
        local_start = start - window_start
        if support < 2 or local_start < 0 or local_start + len(oriented) > len(window):
            continue
        mismatches = sum(a != b for a, b in zip(oriented, window[local_start:local_start+len(oriented)]))
        if mismatches <= 8:
            return oriented, start, is_reverse, mismatches
    return None


def parse_attributes(text: str) -> dict[str, str]:
    return dict(re.findall(r'(\S+) "([^"]*)"', text))


def reconstruct_stop_gained(reference: str) -> list[dict[str, object]]:
    transcripts: dict[str, list[tuple[int, int, str, str]]] = collections.defaultdict(list)
    with gzip.open(GTF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip().split("\t")
            if len(columns) < 9 or columns[2] != "CDS":
                continue
            attributes = parse_attributes(columns[8])
            if attributes.get("gene_type") != "protein_coding":
                continue
            start, end = int(columns[3])-1, int(columns[4])
            if attributes.get("gene_name") == "STXBP1":
                transcripts[attributes["transcript_id"]].append((start, end, columns[6], "STXBP1"))

    confirmations = []
    for transcript, pieces in transcripts.items():
        if not any(start <= POS0 < end for start, end, *_ in pieces):
            continue
        strand = pieces[0][2]
        ordered = sorted(pieces, key=lambda item: item[0], reverse=(strand == "-"))
        cds_parts = []
        offset = None
        cursor = 0
        for start, end, _, _ in ordered:
            part = reference[start:end]
            if strand == "-":
                part = reverse_complement(part)
            if start <= POS0 < end:
                offset = cursor + (POS0-start if strand == "+" else end-1-POS0)
            cds_parts.append(part)
            cursor += end-start
        if offset is None:
            continue
        cds = "".join(cds_parts)
        codon_start = offset - offset % 3
        ref_codon = cds[codon_start:codon_start+3]
        alt_base = ALT if strand == "+" else reverse_complement(ALT)
        changed = list(ref_codon)
        changed[offset % 3] = alt_base
        alt_codon = "".join(changed)
        if ref_codon == "GAA" and alt_codon == "TAA":
            confirmations.append({"transcript_id": transcript, "strand": strand, "ref_codon": ref_codon, "alt_codon": alt_codon, "amino_acid_position": codon_start//3+1})
    return confirmations


def main() -> None:
    reference = load_reference()
    if reference[POS0] != REF:
        raise ValueError("reference allele does not match GRCh38 chr9")
    window_start, window_end = POS0-300, POS0+301
    window = reference[window_start:window_end]
    counts = collections.Counter()
    qualities = collections.defaultdict(list)
    starts = collections.defaultdict(set)
    mismatches_by_allele = collections.defaultdict(list)

    for sequence, quality in read_fastq():
        if len(sequence) < K:
            continue
        hit = local_alignment(sequence, window, window_start)
        if hit is None:
            continue
        oriented, start, is_reverse, mismatches = hit
        if not (start <= POS0 < start + len(oriented)):
            continue
        if is_reverse:
            quality = quality[::-1]
        index = POS0-start
        base_quality = ord(quality[index])-33
        if base_quality < MIN_BASE_Q:
            continue
        base = oriented[index]
        label = "alt" if base == ALT else ("ref" if base == REF else "other")
        counts[label] += 1
        counts[f"{label}_{'reverse' if is_reverse else 'forward'}"] += 1
        qualities[label].append(base_quality)
        starts[label].add(start)
        mismatches_by_allele[label].append(mismatches)

    alt_reads, ref_reads = counts["alt"], counts["ref"]
    total_reads = alt_reads + ref_reads + counts["other"]
    allele_fraction = alt_reads / total_reads
    confirmations = reconstruct_stop_gained(reference)
    if not confirmations:
        raise ValueError("GENCODE v47 stop-gained reconstruction failed")

    with (OUT / "variant.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["chrom","pos","ref","alt","gene","consequence","alt_reads","total_reads","allele_fraction"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"chrom":CHROM,"pos":POS1,"ref":REF,"alt":ALT,"gene":"STXBP1","consequence":"stop_gained","alt_reads":alt_reads,"total_reads":total_reads,"allele_fraction":f"{allele_fraction:.6f}"})

    evidence = {
        "coordinate_system": "GRCh38, 1-based",
        "annotation": "GENCODE v47 primary-assembly chr9 GTF",
        "variant": {"chrom":CHROM,"pos":POS1,"ref":REF,"alt":ALT,"gene":"STXBP1","consequence":"stop_gained"},
        "read_evidence": {
            "alt_reads": alt_reads, "ref_reads": ref_reads, "other_reads": counts["other"], "total_reads": total_reads,
            "allele_fraction": allele_fraction,
            "alt_forward": counts["alt_forward"], "alt_reverse": counts["alt_reverse"],
            "alt_unique_alignment_starts": len(starts["alt"]),
            "alt_base_quality_min": min(qualities["alt"]), "alt_base_quality_median": statistics.median(qualities["alt"]),
            "alt_read_mismatch_median": statistics.median(mismatches_by_allele["alt"]), "alt_read_mismatch_max": max(mismatches_by_allele["alt"]),
        },
        "consequence_reconstruction": {"transcript_count": len(confirmations), "transcripts": confirmations},
        "input_sha256": {
            "deleterious.mutation.q2.R1.fq.gz": "656a689d572750249a8f5a24e35159c9a11bfb63d425a83f2774a81e8379bf29",
            "GRCh38_chr9.fa.gz": "61023a6fb85e19ff8ca6797fc9acc7d9eac751195c51a56e51af087f68579705",
            "gencode.v47.chr9.annotation.gtf.gz": "dacc70b3287ae965e4369f4fcf61b818c942ffd293e1f44aec593af99c172ee6",
        },
    }
    (OUT / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
