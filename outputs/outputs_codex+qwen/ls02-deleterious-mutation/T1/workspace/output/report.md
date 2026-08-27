# Mosaic nonsense SNV analysis - chr9 exome

**Date:** 2026-08-28  
**Task:** Identify the high-confidence mosaic nonsense (stop-gained) SNV in a highly
loss-of-function (LoF)-intolerant protein-coding gene from the supplied chr9 exome reads.

---

## 1. Result

| Field | Value |
|---|---|
| **Variant (GRCh38, 1-based)** | `chr9:127661125 G>T` |
| **Gene** | **STXBP1** (ENSG00000136854.26), chr9q34.11, + strand |
| **Consequence** | **stop_gained (nonsense)** |
| **Canonical transcript** | ENST00000373299.5 (STXBP1-201, Ensembl_canonical) |
| **HGVS (protein)** | **p.Glu117Ter** (p.E117*), codon 117, `GAA>TAA`, exon 6/19 |
| **Allele fraction** | **0.1277** (6 ALT / 47 total unique-molecule reads) |
| **Support** | 6 independent ALT molecules (5 forward / 1 reverse), all baseQ 33, all MAPQ 60 |
| **gnomAD v4 LoF constraint** | pLI = 1.0; o/e LoF = 0.038 (95% CI 0.017-0.099); lof_z = 7.23 |

`STXBP1` is an archetypal highly LoF-intolerant gene (pLI ~ 1, LOEUF ~ 0.04);
heterozygous loss-of-function variants cause STXBP1-related developmental and
epileptic encephalopathy (DEE4, OMIM #612164) via haploinsufficiency. The variant
creates a premature termination codon in exon 6 of the canonical transcript and is
therefore a canonical nonsense (stop-gained) SNV. Its allele fraction of ~12.8% is
significantly inconsistent with germline heterozygosity and far above sequencing
error, establishing it as **mosaic**.

Deliverables written by `output/analysis.py`:

- `output/variant.tsv` - chrom, pos, ref, alt, gene, consequence, alt_reads, total_reads, allele_fraction
- `output/evidence.json` - complete evidence bundle (read-level support, statistics, constraint, provenance)
- `output/analysis.py` - the reproducible end-to-end pipeline
- `output/report.md` - this report

---

## 2. Inputs and reference/annotation versions

| Resource | Version / provenance |
|---|---|
| Reads | `inputs/deleterious.mutation.q2.R1.fq.gz` - single-end chr9 FASTQ, 431,396 primary reads (phred33) |
| Reference genome | **GRCh38 primary assembly**, chromosome `chr9` (138,394,717 bp), from `inputs/reference/GRCh38_chr9.fa.gz` (Broad Institute GATK resource bundle) |
| Gene annotation | **GENCODE v47** (Ensembl 113, release date 2024-07-19), primary assembly, chr9 records only: `inputs/reference/gencode.v47.chr9.annotation.gtf.gz` |
| LoF constraint | **gnomAD v4** gene constraint metrics, fetched live from the gnomAD GraphQL API (`https://gnomad.broadinstitute.org/api`, reference_genome GRCh38) on 2026-08-28; a recorded snapshot is embedded in `analysis.py` as offline fallback |
| Coordinate system | GRCh38, **1-based** (chromosome name `chr9`) |
---

## 3. Methods

All steps are implemented in `output/analysis.py` (Python 3 + subprocess calls to
bwa/samtools/bcftools). Intermediate files live in `work/` and are reused if present.

1. **Alignment.** Reads were aligned to the chr9 reference with `bwa mem 0.7.17-r1188`
   (single-end, with an @RG read group), sorted with `samtools sort`, and PCR duplicates
   were marked with `samtools markdup` (92,995 duplicates; 21.6%). All downstream
   quantification was performed on the duplicate-filtered BAM
   (`samtools view -F 1024`), i.e. counts reflect unique molecules.
   Alignment rate: 429,399/431,396 primary reads mapped (99.54%), essentially all at MAPQ 60.

2. **Germline call set (context).** `bcftools mpileup (-q 20, -Q 20) | bcftools call -mv --ploidy 2`
   produced 18,787 SNVs/indels. Allele-fraction spectrum: 17,426 homozygous-alt sites (AF ~ 1),
   ~755 heterozygous sites centered at AF ~ 0.5, and a small low-AF tail - the expected
   germline + mosaic mixture.

3. **Mosaic-sensitive scan.** `bcftools mpileup (-q 30, -Q 20)` was streamed over all
   16,318,674 covered sites; every site with total depth >= 10 carrying an alternate base
   with >= 3 supporting reads and allele fraction in [0.03, 0.97] was retained (514 scan
   candidates; merged with SNVs from the diploid call set).

4. **Transcript-aware annotation.** A GENCODE v47 protein-coding annotation engine was built:
   for each of the 2,998 protein-coding transcripts (2,500 with valid, in-frame CDS
   reconstruction starting with Met and containing no internal stop), the CDS was
   reconstructed from the reference, and each candidate SNV was projected onto codon
   coordinates (strand- and phase-aware) to classify consequences
   (stop_gained / missense / synonymous / stop_lost / UTR / intron / splice_site).
   A stop codon was called *premature* (nonsense) when it precedes the annotated
   natural stop of the transcript.

5. **High-confidence mosaic filtering.** A stop-gained site is reported only if:
   >= 5 ALT unique molecules, depth >= 20, every ALT base baseQ >= 20 (observed: 33),
   every read MAPQ >= 20 (observed: 60), ALT seen on both strands, all ALT reads have
   distinct alignment starts (no duplicate reads), no other called variant within +/- 50 bp,
   no homopolymer run > 4 in the 21 bp context, binomial P(observed | germline het AF=0.5)
   < 1e-3 (mosaic, not germline), and binomial P(observed | error rate 1%) < 1e-3
   (not sequencing noise).

6. **LoF-intolerance gate.** For every gene carrying a stop-gained candidate, gnomAD v4
   constraint was retrieved; a gene is "highly LoF-intolerant" if pLI >= 0.9 or the
   o/e LoF 95% CI upper bound <= 0.35.

---

## 4. Evidence for the reported variant

### 4.1 Read-level support at chr9:127,661,125

| Metric | Value |
|---|---|
| Reference base / ALT | G / T |
| Depth (unique molecules) | 47 (identical with and without marked duplicates) |
| ALT reads | 6 (5 forward, 1 reverse) |
| REF reads | 41 (23 forward, 18 reverse) |
| Allele fraction | 6/47 = **0.1277** |
| ALT base qualities | all 33 (>= 30 filter: counts unchanged) |
| Mapping qualities | all 60 (unique alignments); MQ0 fraction 0 |
| ALT molecules | 6 distinct alignment starts -> 6 independent molecules |
| Isolation | no other called variant within +/- 50 bp |
| Sequence context | `CCTGTTTAATGAACTGGTAAA` (max homopolymer 3; not a low-complexity site) |
| CIGAR | all ALT reads are full-length matches (no indels/soft clips at the site) |

ALT reads: READ_204796, READ_378102, READ_301814, READ_14027, READ_428597, READ_191105.

### 4.2 Statistical support for mosaicism

| Test | Result |
|---|---|
| Binomial P(X <= 6 \| n = 47, AF = 0.5) - germline het | **8.86e-08** -> rejects germline heterozygosity |
| Binomial P(X >= 6 \| n = 47, error = 0.01) | **7.55e-06** -> rejects sequencing error |
| Binomial P(X >= 6 \| n = 47, error = 0.005) | 1.41e-07 |
| bcftools diploid caller | site NOT called (diploid prior rejects AF ~ 0.13) - consistent with mosaicism |

The allele fraction sits squarely in the mosaic range and is reproducibly estimated
from unique molecules, so 0.1277 is reported (0-1 scale).
### 4.3 Consequence in STXBP1

- Canonical transcript **ENST00000373299.5** (STXBP1-201): exon 6 of 19, coding codon 117,
  `GAA` (Glu) -> `TAA` (stop): **p.Glu117Ter** - a premature termination codon
  ~477 residues upstream of the annotated C-terminus (CDS: 594 codons; GENCODE does not annotate the natural stop codon for this transcript, so any in-CDS stop is premature by definition).
- The same nucleotide change is stop-gained in **17** STXBP1 protein-coding transcripts
  (p.Glu103Ter or p.Glu117Ter depending on alternative first exons), i.e. the nonsense
  effect is transcript-robust.
- STXBP1 (syntaxin-binding protein 1, Munc18-1) is intolerant to haploinsufficiency;
  heterozygous LoF is the established mechanism of STXBP1 encephalopathy.

### 4.4 LoF-intolerance (gnomAD v4, fetched 2026-08-28)

| Gene | pLI | o/e LoF | o/e LoF 95% CI upper | Classification |
|---|---|---|---|---|
| **STXBP1** | **1.00** | **0.038** | **0.099** | **highly LoF-intolerant** |
| GRIN3A | 5.0e-06 | 0.509 | 0.649 | tolerant |
| ALAD | 0.047 | 0.450 | 0.654 | tolerant |
| COL27A1 | 0.0049 | 0.422 | 0.496 | tolerant |
| C5 | ~0 | 0.640 | 0.736 | tolerant |

STXBP1 is the only stop-gained gene that is highly LoF-intolerant; it also has the
strongest constraint of the set by a wide margin (observed LoF variants 3 vs 78.5 expected).

---

## 5. Rejected stop-gained candidates

| Site | Gene | AF | Reads (ALT/DP) | Rejection reason(s) |
|---|---|---|---|---|
| chr9:101670953 C>T | GRIN3A p.G487Ter | 0.550 | 11/20 | AF consistent with germline het (binomial vs 0.5 not significant); gene LoF-tolerant |
| chr9:113391620 A>G | ALAD p.Y56Ter | 0.619 | 13/21 | AF consistent with germline het; gene LoF-tolerant |
| chr9:114196007 C>T | COL27A1 p.R707Ter | 0.174 | 4/23 | only 4 ALT reads (< 5); ALT on reverse strand only (strand-bias artifact); gene LoF-tolerant |
| chr9:121017727 G>A | C5 p.Y544Ter | 0.463 | 19/41 | AF consistent with germline het; gene LoF-tolerant |

A lower-sensitivity scan (ALT >= 2 reads) additionally surfaced chr9:104569892 (OR13C8,
2/30), chr9:128908217 (LRRC8A, 2/20) and chr9:131085628 (LAMC3, 2/62); all have < 5 ALT
reads (below the high-confidence threshold) and none of these genes is highly
LoF-intolerant. No other mosaic nonsense candidate exists in the data.

---

## 6. Tools and software versions

| Tool | Version | Role |
|---|---|---|
| bwa | 0.7.17-r1188 | alignment (bwa mem, single-end) |
| samtools / htslib | 1.19.2 / 1.19 | sorting, markdup, indexing, pileup, read-level evidence |
| bcftools | 1.19 | mpileup, diploid calling |
| Python | 3.12.3 | pipeline orchestration, annotation engine, statistics (stdlib only) |

Tool binaries: Ubuntu 24.04 (noble) `samtools 1.19.2-1build2`, `bwa 0.7.17-7`,
`bcftools 1.19-1build2` packages, executed via WSL2.

---

## 7. Reproducibility

```bash
# from the workspace root; requires bwa, samtools, bcftools on PATH (or set TOOL_DIR)
python3 output/analysis.py --workspace .
# add --redo to recompute all intermediates; --min-alt-reads N to change the ALT-read threshold
```

The script re-uses intermediate files in `work/` (reference index, BAMs, VCF) and
regenerates `output/variant.tsv` and `output/evidence.json`. gnomAD constraint is
fetched live when network is available; otherwise the embedded 2026-08-28 snapshot is used.

## 8. Limitations

- Single-end 150 bp reads: no paired-end consistency checks were possible; confidence
  was instead established via mapping quality, base quality, strand balance, molecule
  independence and binomial testing.
- The mosaic fraction estimate (0.1277) is based on 47 unique molecules; the 95% binomial
  CI is approximately 0.05-0.26.
- Constraint metrics reflect gnomAD v4 population data and are used as a gene-level
  intolerance gate, not as evidence about this specific sample.
- Annotation is limited to GENCODE v47 protein-coding transcripts on chr9; non-coding or
  regulatory consequences were not assessed (out of scope for a nonsense-SNV task).
