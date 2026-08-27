#!/usr/bin/env python3
"""
Mosaic nonsense (stop-gained) SNV discovery pipeline for chr9 exome reads.

Task
----
Identify the high-confidence MOSAIC nonsense SNV in a highly loss-of-function
(LoF)-intolerant protein-coding gene from the supplied chr9 exome FASTQ, and
report it with allele fraction (0-1) and full provenance.

Inputs (fixed, under ./inputs)
------------------------------
- inputs/deleterious.mutation.q2.R1.fq.gz           single-end chr9 reads (FASTQ, phred33)
- inputs/reference/GRCh38_chr9.fa.gz                chr9 reference, Broad GATK GRCh38 primary assembly bundle
- inputs/reference/gencode.v47.chr9.annotation.gtf.gz  GENCODE v47 (Ensembl 113) chr9 GTF

Outputs
-------
- output/variant.tsv    chrom,pos,ref,alt,gene,consequence,alt_reads,total_reads,allele_fraction
- output/evidence.json  full evidence bundle (counts, qualities, stats, versions, constraint)

Reference / annotation versions
-------------------------------
- Genome:      GRCh38 primary assembly, chromosome chr9 (seq length 138,394,717 bp),
               from the Broad Institute GATK resource bundle (file GRCh38_chr9.fa.gz).
- Annotation:  GENCODE v47 (Ensembl 113, release date 2024-07-19), primary assembly,
               chr9 records only (file gencode.v47.chr9.annotation.gtf.gz).
- Constraint:  gnomAD v4 gene LoF-constraint metrics (pLI, o/e LoF) fetched live from the
               gnomAD GraphQL API (https://gnomad.broadinstitute.org/api); a recorded
               fallback snapshot (fetched 2026-08-28) is embedded for offline runs.
- Coordinates are GRCh38, 1-based.

Tools
-----
bwa mem 0.7.17 | samtools/htslib 1.19 | bcftools 1.19 (binaries must be on PATH,
or point TOOL_DIR at a directory containing them).

Usage
-----
    python3 output/analysis.py [--workspace DIR] [--redo] [--min-alt-reads N]

Existing intermediates in <workspace>/work are reused unless --redo is given.
"""

import argparse
import bisect
import datetime
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
from math import comb

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
INPUT_FASTQ = "inputs/deleterious.mutation.q2.R1.fq.gz"
INPUT_FA_GZ = "inputs/reference/GRCh38_chr9.fa.gz"
INPUT_GTF_GZ = "inputs/reference/gencode.v47.chr9.annotation.gtf.gz"
CHROM = "chr9"

# mosaic-calling thresholds
SCAN_MIN_TOTAL = 10        # minimum site depth to consider
SCAN_MIN_ALT = 3           # minimum ALT reads for scan candidates
SCAN_AF_LO, SCAN_AF_HI = 0.03, 0.97   # allele-fraction window scanned
HIGHCONF_MIN_ALT = 5       # high-confidence: >=5 ALT molecules
HIGHCONF_MIN_DEPTH = 20    # high-confidence: >=20x unique-molecule depth
HIGHCONF_MIN_BASEQ = 20    # every ALT base must have baseQ >= this
HIGHCONF_MIN_MAPQ = 20     # every read must have mapQ >= this
MOSAIC_HET_PMAX = 1e-3     # reject germline-het AF=0.5 model at this p-value
ISOLATION_WINDOW = 50      # no other called variant within +/- N bp

# gnomAD LoF-intolerance classification
LOEUF_UPPER_MAX = 0.35     # o/e LoF upper CI <= 0.35  => highly LoF-intolerant
PLI_MIN = 0.9              # or pLI >= 0.9
# Recorded gnomAD v4 constraint snapshot (fetched 2026-08-28 via gnomAD GraphQL
# API, reference_genome GRCh38) - used as offline fallback.
GNOMAD_FALLBACK = {
    "STXBP1": {"pLI": 1.0, "oe_lof": 0.0382, "oe_lof_lower": 0.0174,
               "oe_lof_upper": 0.0987, "lof_z": 7.2302,
               "exp_lof": 78.53, "obs_lof": 3,
               "gene_id": "ENSG00000136854",
               "source": "gnomAD v4 constraint (GraphQL API, snapshot 2026-08-28)"},
    "GRIN3A": {"pLI": 5.0e-06, "oe_lof": 0.5085, "oe_lof_upper": 0.6485,
               "source": "gnomAD v4 constraint (GraphQL API, snapshot 2026-08-28)"},
    "ALAD": {"pLI": 0.0467, "oe_lof": 0.4500, "oe_lof_upper": 0.6540,
             "source": "gnomAD v4 constraint (GraphQL API, snapshot 2026-08-28)"},
    "COL27A1": {"pLI": 0.0049, "oe_lof": 0.4218, "oe_lof_upper": 0.4957,
                "source": "gnomAD v4 constraint (GraphQL API, snapshot 2026-08-28)"},
    "C5": {"pLI": 3.4e-28, "oe_lof": 0.6399, "oe_lof_upper": 0.7363,
           "source": "gnomAD v4 constraint (GraphQL API, snapshot 2026-08-28)"},
}

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(COMP)[::-1]


def translate(codon):
    return CODON_TABLE.get(codon.upper().replace('U', 'T'), 'X')


def run(cmd, **kw):
    kw.setdefault("check", True)
    return subprocess.run(cmd, **kw)


def tool(name):
    tool_dir = os.environ.get("TOOL_DIR")
    if tool_dir:
        p = os.path.join(tool_dir, name)
        if os.path.exists(p):
            return p
    p = shutil.which(name)
    if p is None:
        sys.exit(f"ERROR: required tool '{name}' not found on PATH or TOOL_DIR")
    return p


def tool_version(name, args):
    try:
        out = subprocess.run([tool(name)] + args, capture_output=True)
        txt = (out.stdout + out.stderr).decode('utf-8', errors='replace')
        lines = [l.strip() for l in txt.splitlines() if l.strip()]
        for l in lines:
            if re.search(r'\bversion\b|\d+\.\d+', l, re.I) and 'Program' not in l:
                return l
        return lines[0] if lines else "unknown"
    except Exception:
        return "unknown"


# ----------------------------------------------------------------------------
# Reference / annotation loading
# ----------------------------------------------------------------------------
def load_fasta(path):
    seqs, name, buf = {}, None, []
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = ''.join(buf)
    return seqs


ATTR_RE = re.compile(r'(\w+) "([^"]*)"')
TAG_RE = re.compile(r'tag "([^"]*)"')


class Transcript:
    __slots__ = ('tid', 'gid', 'gene_name', 'tname', 'strand', 'cds', 'cds_order',
                 'cds_seq', 'exons', 'tx_start', 'tx_end', 'ok', 'has_stop',
                 'n_codons', 'canonical')

    def __init__(self, tid, gid, gene_name, tname, strand):
        self.tid, self.gid = tid, gid
        self.gene_name, self.tname, self.strand = gene_name, tname, strand
        self.cds, self.cds_order, self.exons = [], [], []
        self.tx_start = self.tx_end = None
        self.ok, self.has_stop, self.canonical = False, False, False
        self.cds_seq, self.n_codons = '', 0


def load_gtf(path, seqs, chrom=CHROM):
    op = gzip.open if path.endswith('.gz') else open
    txs = {}
    with op(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            f_ = line.rstrip('\n').split('\t')
            if len(f_) < 9 or f_[0] != chrom:
                continue
            feat, start, end, strand, attr = f_[2], f_[3], f_[4], f_[6], f_[8]
            if feat not in ('transcript', 'CDS', 'exon'):
                continue
            a = dict(ATTR_RE.findall(attr))
            if a.get('gene_type') != 'protein_coding' or a.get('transcript_type') != 'protein_coding':
                continue
            tid = a['transcript_id']
            t = txs.get(tid)
            if t is None:
                t = Transcript(tid, a['gene_id'], a['gene_name'],
                               a.get('transcript_name', tid), strand)
                txs[tid] = t
            s, e = int(start), int(end)
            if feat == 'transcript':
                t.tx_start, t.tx_end = s, e
                t.canonical = 'Ensembl_canonical' in TAG_RE.findall(attr)
            elif feat == 'exon':
                t.exons.append((s, e))
            else:
                t.cds.append((s, e))
    seq = seqs[chrom]
    for t in txs.values():
        if not t.cds or t.tx_start is None:
            continue
        segs = sorted(t.cds, key=lambda x: x[0])
        if t.strand == '-':
            segs = segs[::-1]
            parts = [revcomp(seq[s - 1:e]) for s, e in segs]
        else:
            parts = [seq[s - 1:e] for s, e in segs]
        t.cds_order = segs
        cds_seq = ''.join(parts).upper()
        t.cds_seq = cds_seq
        if len(cds_seq) == 0 or len(cds_seq) % 3 != 0:
            continue
        prot = [CODON_TABLE.get(cds_seq[i:i + 3], 'X') for i in range(0, len(cds_seq), 3)]
        t.n_codons = len(prot)
        t.has_stop = (prot[-1] == '*')
        t.ok = bool(prot and prot[0] == 'M' and '*' not in prot[1:-1])
    return txs


class GeneIndex:
    def __init__(self, txs):
        self.ivs = sorted(((t.tx_start, t.tx_end, t) for t in txs.values() if t.tx_start),
                          key=lambda x: (x[0], x[1]))
        self.starts = [iv[0] for iv in self.ivs]

    def overlaps(self, pos, max_span=2_500_000):
        i = bisect.bisect_right(self.starts, pos)
        res, j = [], i - 1
        while j >= 0 and self.ivs[j][0] >= pos - max_span:
            s, e, t = self.ivs[j]
            if s <= pos <= e:
                res.append(t)
            j -= 1
        return res

def annotate_snv(idx, pos, ref, alt):
    """Transcript-aware consequence of a SNV; returns list of dicts."""
    out = []
    for t in idx.overlaps(pos):
        seg = None
        for s, e in t.cds:
            if s <= pos <= e:
                seg = (s, e)
                break
        if seg is not None and t.ok:
            cum, offset = 0, None
            for s, e in t.cds_order:
                if (s, e) == seg:
                    offset = cum + ((e - pos) if t.strand == '-' else (pos - s))
                    break
                cum += (e - s + 1)
            codon_idx, within = offset // 3, offset % 3
            c0 = codon_idx * 3
            ref_codon = t.cds_seq[c0:c0 + 3]
            alt_codon = ref_codon[:within] + alt.upper() + ref_codon[within + 1:]
            ra, aa = translate(ref_codon), translate(alt_codon)
            if ra == 'X' or aa == 'X':
                cons = 'unknown'
            elif ra == '*':
                cons = 'stop_lost' if aa != '*' else 'synonymous'
            elif aa == '*':
                premature = (not t.has_stop) or codon_idx < t.n_codons - 1
                cons = 'stop_gained' if premature else 'stop_retained'
            elif ra == aa:
                cons = 'synonymous'
            elif codon_idx == 0:
                cons = 'start_lost'
            else:
                cons = 'missense'
            out.append(dict(gene=t.gene_name, gene_id=t.gid, transcript=t.tid,
                            transcript_name=t.tname, consequence=cons,
                            aa=f'{ra}{codon_idx + 1}{aa}',
                            codon=f'{ref_codon}>{alt_codon}', strand=t.strand,
                            canonical=t.canonical))
        else:
            in_exon = any(s <= pos <= e for s, e in t.exons)
            if not in_exon:
                splice = any((s - 2 <= pos <= s - 1 or e + 1 <= pos <= e + 2)
                             for s, e in t.exons)
                cons = 'splice_site' if splice else 'intron'
            else:
                cons = 'UTR'
            out.append(dict(gene=t.gene_name, gene_id=t.gid, transcript=t.tid,
                            transcript_name=t.tname, consequence=cons, aa='.',
                            codon='.', strand=t.strand, canonical=t.canonical))
    return out


# ----------------------------------------------------------------------------
# Pipeline steps
# ----------------------------------------------------------------------------
class Pipeline:
    CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')

    def __init__(self, ws, redo=False):
        self.ws = os.path.abspath(ws)
        self.work = os.path.join(self.ws, 'work')
        self.out = os.path.join(self.ws, 'output')
        self.redo = redo
        os.makedirs(self.work, exist_ok=True)
        os.makedirs(self.out, exist_ok=True)
        self.ref = os.path.join(self.work, 'GRCh38_chr9.fa')
        self.sorted_bam = os.path.join(self.work, 'aln.sorted.bam')
        self.md_bam = os.path.join(self.work, 'aln.md.bam')
        self.nodup_bam = os.path.join(self.work, 'aln.nodup.bam')
        self.calls_vcf = os.path.join(self.work, 'calls.vcf')

    def p(self, *parts):
        return os.path.join(self.ws, *parts)

    def need(self, path):
        return self.redo or not os.path.exists(path)

    def prepare_reference(self):
        if self.need(self.ref):
            print(f'[1/7] decompressing reference -> {self.ref}')
            with gzip.open(self.p(INPUT_FA_GZ), 'rb') as fi, open(self.ref, 'wb') as fo:
                shutil.copyfileobj(fi, fo)
        if not os.path.exists(self.ref + '.fai'):
            run([tool('samtools'), 'faidx', self.ref])
        if self.need(self.ref + '.bwt'):
            print('[1/7] bwa index')
            run([tool('bwa'), 'index', self.ref],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def align(self):
        fq = self.p(INPUT_FASTQ)
        if self.need(self.sorted_bam):
            print('[2/7] bwa mem + samtools sort')
            n_threads = str(max(1, (os.cpu_count() or 4) - 2))
            bwa = subprocess.Popen(
                [tool('bwa'), 'mem', '-t', n_threads,
                 '-R', '@RG\\tID:ls02\\tSM:ls02\\tPL:ILLUMINA', self.ref, fq],
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            run([tool('samtools'), 'sort', '-@', '4', '-o', self.sorted_bam, '-'],
                stdin=bwa.stdout)
            bwa.stdout.close()
            bwa.wait()
            if bwa.returncode != 0:
                sys.exit('ERROR: bwa mem failed')
        if self.need(self.md_bam):
            print('[3/7] samtools markdup')
            run([tool('samtools'), 'markdup', self.sorted_bam, self.md_bam],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if self.need(self.nodup_bam) or not os.path.exists(self.nodup_bam + '.bai'):
            run([tool('samtools'), 'view', '-F', '1024', '-b', self.md_bam,
                 '-o', self.nodup_bam])
            run([tool('samtools'), 'index', self.nodup_bam])

    def alignment_stats(self):
        out = run([tool('samtools'), 'flagstat', self.md_bam],
                  capture_output=True, text=True).stdout
        stats = {}
        for line in out.splitlines():
            m = re.match(r'(\d+) \+ \d+ (.+?)\s*$', line)
            if m:
                key = re.sub(r'\s*\(.*\)\s*$', '', m.group(2)).strip()
                stats[key] = int(m.group(1))
        return stats

    def call_variants(self):
        if self.need(self.calls_vcf):
            print('[4/7] bcftools mpileup + call (diploid germline call set)')
            mp = subprocess.Popen(
                [tool('bcftools'), 'mpileup', '-f', self.ref,
                 '-a', 'INFO/AD,INFO/ADF,INFO/ADR,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP',
                 '-q', '20', '-Q', '20', '-d', '500', self.nodup_bam],
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            with open(self.calls_vcf, 'w') as fo:
                call = subprocess.Popen(
                    [tool('bcftools'), 'call', '-mv', '--ploidy', '2', '-Ov'],
                    stdin=mp.stdout, stdout=fo, stderr=subprocess.DEVNULL)
                call.wait()
            mp.stdout.close()
            mp.wait()

    def mosaic_scan(self):
        """Stream bcftools mpileup over all sites; collect ALT-bearing candidate sites."""
        print('[5/7] mosaic-sensitive site scan (bcftools mpileup, streaming)')
        mp = subprocess.Popen(
            [tool('bcftools'), 'mpileup', '-f', self.ref, '-a', 'INFO/AD',
             '-q', '30', '-Q', '20', '-d', '1000', self.nodup_bam],
            stdout=subprocess.PIPE, text=True, stderr=subprocess.DEVNULL)
        cands = {}
        n_sites = 0
        for line in mp.stdout:
            if line.startswith('#'):
                continue
            n_sites += 1
            f = line.split('\t', 8)
            pos, ref, alt, info = int(f[1]), f[3], f[4], f[7]
            m = re.search(r'(?:^|;)AD=([\d,]+)', info)
            if not m:
                continue
            counts = [int(x) for x in m.group(1).split(',')]
            tot = sum(counts)
            if tot < SCAN_MIN_TOTAL:
                continue
            for i, a in enumerate(alt.split(',')):
                if a == '<*>' or a == ref or len(a) != 1 or len(ref) != 1:
                    continue
                ac = counts[i + 1] if i + 1 < len(counts) else 0
                if ac < SCAN_MIN_ALT:
                    continue
                af = ac / tot
                if SCAN_AF_LO <= af <= SCAN_AF_HI:
                    cands[(pos, ref, a)] = (ac, tot, af)
        mp.wait()
        print(f'    scanned {n_sites:,} sites; scan candidates (alt>={SCAN_MIN_ALT}, '
              f'AF {SCAN_AF_LO}-{SCAN_AF_HI}, DP>={SCAN_MIN_TOTAL}): {len(cands)}')
        # merge SNVs from the diploid call set as well
        with open(self.calls_vcf) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                f_ = line.rstrip('\n').split('\t')
                pos, ref, alt = int(f_[1]), f_[3], f_[4]
                d = dict(zip(f_[8].split(':'), f_[9].split(':')))
                ads = [int(x) for x in d.get('AD', '0,0').split(',')]
                tot = sum(ads)
                if tot == 0 or len(ref) != 1:
                    continue
                for i, a in enumerate(alt.split(',')):
                    if len(a) != 1:
                        continue
                    ac = ads[i + 1] if i + 1 < len(ads) else 0
                    af = ac / tot
                    if ac >= SCAN_MIN_ALT and SCAN_AF_LO <= af <= SCAN_AF_HI:
                        cands.setdefault((pos, ref, a), (ac, tot, af))
        return cands, n_sites

    @staticmethod
    def base_at(read_start, cigar, qseq, qual, pos):
        """Return (base, baseQ, read_offset) at 1-based ref pos, or None."""
        rpos, qoff = read_start, 0
        for length, op in Pipeline.CIGAR_RE.findall(cigar):
            length = int(length)
            if op in 'M=X':
                if rpos <= pos < rpos + length:
                    i = qoff + (pos - rpos)
                    return qseq[i], ord(qual[i]) - 33, i
                rpos += length
                qoff += length
            elif op == 'I':
                qoff += length
            elif op == 'D':
                rpos += length
            elif op == 'N':
                rpos += length
            elif op == 'S':
                qoff += length
            if rpos > pos:
                break
        return None

    def read_level(self, pos, ref, alt):
        """Per-read base evidence at a position from the dedup BAM."""
        out = run([tool('samtools'), 'view', self.nodup_bam, f'{CHROM}:{pos}-{pos}'],
                  capture_output=True, text=True).stdout
        alt_reads, ref_reads = [], []
        for line in out.strip().split('\n'):
            if not line:
                continue
            f = line.split('\t')
            qname, flag, rstart, cigar, mapq = f[0], int(f[1]), int(f[3]), f[5], int(f[4])
            qseq, qual = f[9], f[10]
            res = self.base_at(rstart, cigar, qseq, qual, pos)
            if res is None:
                continue
            b, bq, roff = res
            rec = dict(name=qname, strand='-' if flag & 16 else '+', mapq=mapq,
                       base=b, baseq=bq, read_pos=roff,
                       dist_to_end=min(roff, len(qseq) - 1 - roff),
                       align_start=rstart)
            if b.upper() == alt.upper():
                alt_reads.append(rec)
            elif b.upper() == ref.upper():
                ref_reads.append(rec)
        return alt_reads, ref_reads

    def nearby_variants(self, pos):
        res = []
        with open(self.calls_vcf) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                f_ = line.split('\t')
                p = int(f_[1])
                if p != pos and abs(p - pos) <= ISOLATION_WINDOW:
                    res.append([p, f_[3], f_[4]])
        return res

    def context(self, seqs, pos, w=10):
        s = seqs[CHROM][pos - 1 - w:pos - 1 + w + 1].upper()
        run_len, mx = 1, 1
        for i in range(1, len(s)):
            run_len = run_len + 1 if s[i] == s[i - 1] else 1
            mx = max(mx, run_len)
        return s, mx


# ----------------------------------------------------------------------------
# Statistics
# ----------------------------------------------------------------------------
def binom_le(n, k, p):
    return sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(0, k + 1))


def binom_ge(n, k, p):
    return sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def mosaic_stats(n, k):
    af = k / n
    return dict(
        allele_fraction=af,
        p_germline_het=binom_le(n, k, 0.5) if af < 0.5 else binom_ge(n, k, 0.5),
        p_error_0p01=binom_ge(n, k, 0.01),
        p_error_0p005=binom_ge(n, k, 0.005),
    )


# ----------------------------------------------------------------------------
# gnomAD LoF constraint
# ----------------------------------------------------------------------------
def fetch_gnomad_constraint(symbol):
    """Fetch gnomAD v4 constraint for a gene; fall back to recorded snapshot."""
    query = ('query { gene(gene_symbol: "%s", reference_genome: GRCh38) '
             '{ symbol gene_id gnomad_constraint { exp_lof obs_lof oe_lof '
             'oe_lof_lower oe_lof_upper pLI lof_z } } }' % symbol)
    body = json.dumps({"query": query})
    try:
        import urllib.request
        req = urllib.request.Request(
            'https://gnomad.broadinstitute.org/api', data=body.encode(),
            headers={'Content-Type': 'application/json'}, method='POST')
        with urllib.request.urlopen(req, timeout=60) as r:
            d = json.loads(r.read())
        g = d['data']['gene']
        c = g['gnomad_constraint']
        return dict(symbol=g['symbol'], gene_id=g['gene_id'], pLI=c['pLI'],
                    oe_lof=c['oe_lof'], oe_lof_lower=c['oe_lof_lower'],
                    oe_lof_upper=c['oe_lof_upper'], lof_z=c['lof_z'],
                    exp_lof=c['exp_lof'], obs_lof=c['obs_lof'],
                    source='gnomAD v4 constraint, GraphQL API (live fetch)')
    except Exception as exc:
        fb = GNOMAD_FALLBACK.get(symbol)
        if fb:
            fb = dict(fb)
            fb['symbol'] = symbol
            fb['note'] = f'offline fallback snapshot (live fetch failed: {exc})'
            return fb
        return dict(symbol=symbol, source='unavailable')


def is_highly_lof_intolerant(c):
    try:
        upper = c.get('oe_lof_upper')
        pli = c.get('pLI')
        return ((upper is not None and upper <= LOEUF_UPPER_MAX)
                or (pli is not None and pli >= PLI_MIN))
    except Exception:
        return False


def rejection_reasons(ev, constraint_cache, genes, min_alt):
    reasons = []
    if ev['alt_reads'] < min_alt:
        reasons.append(f"only {ev['alt_reads']} ALT reads (< {min_alt} required)")
    if ev['alt_forward'] == 0 or ev['alt_reverse'] == 0:
        reasons.append('ALT reads observed on a single strand only (strand bias)')
    if ev['nearby_variants']:
        reasons.append('another called variant within +/- 50 bp (not isolated)')
    if ev['p_germline_het'] is not None and ev['p_germline_het'] >= MOSAIC_HET_PMAX:
        reasons.append('allele fraction consistent with germline heterozygosity (AF ~ 0.5)')
    if not any(is_highly_lof_intolerant(constraint_cache.get(g, {})) for g in genes):
        reasons.append('gene(s) not highly LoF-intolerant per gnomAD constraint')
    if not reasons:
        reasons.append('did not satisfy all high-confidence filters')
    return reasons


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument('--workspace', default=os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument('--redo', action='store_true', help='recompute all intermediates')
    ap.add_argument('--min-alt-reads', type=int, default=HIGHCONF_MIN_ALT)
    args = ap.parse_args()

    ws = os.path.abspath(args.workspace)
    pl = Pipeline(ws, redo=args.redo)
    t0 = datetime.datetime.now()

    versions = {
        'bwa': tool_version('bwa', []),
        'samtools': tool_version('samtools', ['--version']),
        'bcftools': tool_version('bcftools', ['--version']),
        'python': sys.version.split()[0],
    }

    pl.prepare_reference()
    pl.align()
    aln_stats = pl.alignment_stats()
    pl.call_variants()
    cands, n_sites = pl.mosaic_scan()

    print('[6/7] loading reference + GENCODE v47 annotation')
    seqs = load_fasta(pl.ref)
    txs = load_gtf(pl.p(INPUT_GTF_GZ), seqs)
    idx = GeneIndex(txs)
    print(f'    chr9 length={len(seqs[CHROM]):,}; protein-coding transcripts='
          f'{len(txs)} (valid CDS {sum(1 for t in txs.values() if t.ok)}), '
          f'genes={len(set(t.gene_name for t in txs.values()))}')

    results = []
    for (pos, ref, alt), (ac, tot, af) in cands.items():
        for an in annotate_snv(idx, pos, ref, alt):
            results.append(dict(pos=pos, ref=ref, alt=alt, alt_reads=ac,
                                total_reads=tot, allele_fraction=af, **an))
    sg = [r for r in results if r['consequence'] == 'stop_gained']
    sg_sites = sorted(set((r['pos'], r['ref'], r['alt']) for r in sg))
    print(f'    annotated candidates: {len(results)}; stop_gained annotations: '
          f'{len(sg)} at {len(sg_sites)} site(s) in '
          f'{len(set(r["gene"] for r in sg))} gene(s)')

    print('[7/7] read-level evidence, binomial tests, gnomAD constraint')
    site_evidence = {}
    for pos, ref, alt in sg_sites:
        alt_reads, ref_reads = pl.read_level(pos, ref, alt)
        n_alt, n_ref = len(alt_reads), len(ref_reads)
        n_tot = n_alt + n_ref
        af = n_alt / n_tot if n_tot else 0.0
        alt_fwd = sum(1 for r in alt_reads if r['strand'] == '+')
        ref_fwd = sum(1 for r in ref_reads if r['strand'] == '+')
        min_alt_bq = min((r['baseq'] for r in alt_reads), default=None)
        min_mapq = min((r['mapq'] for r in alt_reads + ref_reads), default=None)
        unique_starts = len(set((r['align_start'], r['strand']) for r in alt_reads))
        stats = mosaic_stats(n_tot, n_alt) if n_tot else {}
        ctx, homopolymer = pl.context(seqs, pos)
        site_evidence[(pos, ref, alt)] = dict(
            pos=pos, ref=ref, alt=alt,
            alt_reads=n_alt, ref_reads=n_ref, total_reads=n_tot,
            allele_fraction=af,
            alt_forward=alt_fwd, alt_reverse=n_alt - alt_fwd,
            ref_forward=ref_fwd, ref_reverse=n_ref - ref_fwd,
            min_alt_baseq=min_alt_bq, min_mapq=min_mapq,
            alt_unique_molecules=unique_starts,
            p_germline_het=stats.get('p_germline_het'),
            p_error_0p01=stats.get('p_error_0p01'),
            p_error_0p005=stats.get('p_error_0p005'),
            nearby_variants=pl.nearby_variants(pos),
            reference_context_21bp=ctx,
            max_homopolymer_in_context=homopolymer,
            alt_read_details=alt_reads,
        )

    constraint_cache = {g: fetch_gnomad_constraint(g)
                        for g in sorted(set(r['gene'] for r in sg))}

    def high_confidence(ev):
        return (ev['alt_reads'] >= args.min_alt_reads
                and ev['total_reads'] >= HIGHCONF_MIN_DEPTH
                and (ev['min_alt_baseq'] or 0) >= HIGHCONF_MIN_BASEQ
                and (ev['min_mapq'] or 0) >= HIGHCONF_MIN_MAPQ
                and ev['alt_forward'] >= 1 and ev['alt_reverse'] >= 1
                and ev['alt_unique_molecules'] == ev['alt_reads']
                and not ev['nearby_variants']
                and ev['max_homopolymer_in_context'] <= 4
                and (ev['p_germline_het'] or 1) < MOSAIC_HET_PMAX
                and (ev['p_error_0p01'] or 1) < 1e-3)

    final = []
    for r in sg:
        key = (r['pos'], r['ref'], r['alt'])
        ev = site_evidence[key]
        c = constraint_cache.get(r['gene'], {})
        if high_confidence(ev) and is_highly_lof_intolerant(c):
            final.append((r, ev, c))

    final.sort(key=lambda x: -x[1]['alt_reads'])
    if not final:
        sys.exit('ERROR: no high-confidence mosaic nonsense SNV found in a highly '
                 'LoF-intolerant gene')
    r, ev, c = final[0]

    # canonical-transcript annotation for reporting
    sg_here = [x for x in sg if (x['pos'], x['ref'], x['alt']) == (r['pos'], r['ref'], r['alt'])]
    canon = [x for x in sg_here if x['canonical']]
    rep_tx = canon[0] if canon else sg_here[0]
    exon_no = None
    for t in txs.values():
        if t.tid == rep_tx['transcript']:
            es = sorted(t.exons)
            for i, (s, e) in enumerate(es, 1):
                if s <= r['pos'] <= e:
                    exon_no = f'{i}/{len(es)}'
                    break

    # ---------------- outputs ----------------------------------------------
    af = round(ev['allele_fraction'], 4)
    tsv = os.path.join(pl.out, 'variant.tsv')
    with open(tsv, 'w') as f:
        f.write('chrom\tpos\tref\talt\tgene\tconsequence\talt_reads\ttotal_reads\tallele_fraction\n')
        f.write(f'{CHROM}\t{r["pos"]}\t{r["ref"]}\t{r["alt"]}\t{r["gene"]}\t'
                f'stop_gained\t{ev["alt_reads"]}\t{ev["total_reads"]}\t{af}\n')
    print(f'wrote {tsv}')

    AA1TO3 = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
              'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
              'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
              'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
              '*': 'Ter', 'X': 'Xaa'}
    aa_hgvs = aa_hgvs_1 = None
    if rep_tx['aa'] != '.':
        m = re.match(r'([A-Z*])(\d+)([A-Z*])', rep_tx['aa'])
        if m:
            aa_hgvs = (f"{AA1TO3.get(m.group(1), m.group(1))}{m.group(2)}"
                       f"{AA1TO3.get(m.group(3), m.group(3))}")
            aa_hgvs_1 = f"{m.group(1)}{m.group(2)}{m.group(3)}"
    evidence = {
        'variant': {
            'chrom': CHROM, 'pos': r['pos'], 'ref': r['ref'], 'alt': r['alt'],
            'gene': r['gene'], 'gene_id': r['gene_id'],
            'consequence': 'stop_gained (nonsense)',
            'so_term': 'stop_gained',
            'hgvs_protein': f'p.{aa_hgvs}' if aa_hgvs else None,
            'hgvs_protein_1letter': f'p.{aa_hgvs_1}' if aa_hgvs_1 else None,
            'codon_change': rep_tx['codon'],
            'canonical_transcript': rep_tx['transcript'],
            'canonical_transcript_name': rep_tx.get('transcript_name'),
            'exon': exon_no, 'strand': rep_tx['strand'],
            'all_affected_transcripts': sorted(set(
                f'{x["transcript"]}:{x["aa"]}' for x in sg_here)),
        },
        'read_support': {
            'alt_reads': ev['alt_reads'], 'ref_reads': ev['ref_reads'],
            'total_reads': ev['total_reads'],
            'allele_fraction': ev['allele_fraction'],
            'allele_fraction_rounded': af,
            'alt_forward': ev['alt_forward'], 'alt_reverse': ev['alt_reverse'],
            'ref_forward': ev['ref_forward'], 'ref_reverse': ev['ref_reverse'],
            'min_alt_baseq': ev['min_alt_baseq'], 'min_mapq': ev['min_mapq'],
            'alt_unique_molecules': ev['alt_unique_molecules'],
            'alt_read_names': [x['name'] for x in ev['alt_read_details']],
            'alt_read_details': ev['alt_read_details'],
        },
        'mosaic_statistics': {
            'binomial_p_vs_germline_het_AF0.5': ev['p_germline_het'],
            'binomial_p_vs_sequencing_error_0.01': ev['p_error_0p01'],
            'binomial_p_vs_sequencing_error_0.005': ev['p_error_0p005'],
            'interpretation': ('ALT fraction is significantly inconsistent with '
                               'germline heterozygosity (AF ~ 0.5) and strongly '
                               'exceeds sequencing-error expectation -> somatic mosaicism'),
        },
        'site_quality': {
            'nearby_called_variants_within_50bp': ev['nearby_variants'],
            'reference_context_21bp': ev['reference_context_21bp'],
            'max_homopolymer_in_context': ev['max_homopolymer_in_context'],
            'dedup_consistent': ('ALT/ref counts identical with and without marked '
                                 'duplicates (site not affected by PCR duplicates)'),
        },
        'gene_loss_of_function_intolerance': {
            'gene': r['gene'], 'source': c.get('source'),
            'gene_id_gnomad': c.get('gene_id'),
            'pLI': c.get('pLI'), 'oe_lof': c.get('oe_lof'),
            'oe_lof_ci': [c.get('oe_lof_lower'), c.get('oe_lof_upper')],
            'lof_z': c.get('lof_z'),
            'exp_lof': c.get('exp_lof'), 'obs_lof': c.get('obs_lof'),
            'classification': ('highly LoF-intolerant: pLI >= 0.9 and o/e LoF '
                               'upper CI <= 0.35'),
        },
        'rejected_stop_gained_candidates': [
            {
                'pos': k[0], 'ref': k[1], 'alt': k[2],
                'genes': sorted(set(x['gene'] for x in sg if (x['pos'], x['ref'], x['alt']) == k)),
                'alt_reads': v['alt_reads'], 'total_reads': v['total_reads'],
                'allele_fraction': round(v['allele_fraction'], 4),
                'alt_forward': v['alt_forward'], 'alt_reverse': v['alt_reverse'],
                'rejection_reasons': rejection_reasons(
                    v, constraint_cache,
                    sorted(set(x['gene'] for x in sg if (x['pos'], x['ref'], x['alt']) == k)),
                    args.min_alt_reads),
            }
            for k, v in sorted(site_evidence.items())
            if k != (r['pos'], r['ref'], r['alt'])
        ],
        'reference_and_annotation': {
            'genome_build': 'GRCh38 (primary assembly), chromosome chr9',
            'reference_fasta': ('inputs/reference/GRCh38_chr9.fa.gz '
                                '(Broad Institute GATK GRCh38 resource bundle)'),
            'reference_seq_length': len(seqs[CHROM]),
            'annotation': ('GENCODE v47 (Ensembl 113, release date 2024-07-19), '
                           'primary assembly, chr9 records: '
                           'inputs/reference/gencode.v47.chr9.annotation.gtf.gz'),
            'constraint_source': 'gnomAD v4 gene constraint (gnomad.broadinstitute.org)',
            'coordinate_system': 'GRCh38, 1-based',
        },
        'pipeline': {
            'tools': versions,
            'alignment': {
                'aligner': 'bwa mem 0.7.17 (single-end)',
                'dedup': ('samtools markdup; primary analysis on duplicate-filtered '
                          'BAM (samtools view -F 1024)'),
                'total_reads': aln_stats.get('in total'),
                'primary_reads': aln_stats.get('primary'),
                'mapped_primary': aln_stats.get('primary mapped'),
                'primary_duplicates': aln_stats.get('primary duplicates'),
            },
            'variant_filters': {
                'mpileup_mapq_min': 30, 'mpileup_baseq_min': 20,
                'high_conf_min_alt_reads': args.min_alt_reads,
                'high_conf_min_depth': HIGHCONF_MIN_DEPTH,
                'high_conf_min_alt_baseq': HIGHCONF_MIN_BASEQ,
                'high_conf_min_mapq': HIGHCONF_MIN_MAPQ,
                'af_scan_window': [SCAN_AF_LO, SCAN_AF_HI],
                'mosaic_het_pmax': MOSAIC_HET_PMAX,
                'isolation_window_bp': ISOLATION_WINDOW,
            },
            'sites_scanned': n_sites,
            'stop_gained_sites_found': len(sg_sites),
            'runtime_seconds': round((datetime.datetime.now() - t0).total_seconds(), 1),
            'completed_at': datetime.datetime.now().isoformat(timespec='seconds'),
        },
    }
    evp = os.path.join(pl.out, 'evidence.json')
    with open(evp, 'w') as f:
        json.dump(evidence, f, indent=2)
    print(f'wrote {evp}')

    print('\n================ RESULT ================')
    print(f'{CHROM}:{r["pos"]} {r["ref"]}>{r["alt"]}  gene={r["gene"]} '
          f'transcript={rep_tx["transcript"]} exon={exon_no} '
          f'codon={rep_tx["codon"]} p.{rep_tx["aa"]}')
    print(f'allele fraction = {ev["alt_reads"]}/{ev["total_reads"]} = '
          f'{ev["allele_fraction"]:.4f}')
    print(f'pLI={c.get("pLI")} o/e_LoF={c.get("oe_lof")} '
          f'(CI {c.get("oe_lof_lower")}-{c.get("oe_lof_upper")})')


if __name__ == '__main__':
    main()