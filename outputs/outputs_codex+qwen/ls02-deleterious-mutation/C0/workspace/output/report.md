# Mosaic nonsense SNV in a highly LoF-intolerant gene — analysis report

**Called variant:** `chr9:127661125 G>T` in **STXBP1** — ENST00000373299.5 c.349 p.E117* (GAA > TAA), nonsense (stop_gained).

| field | value |
|---|---|
| chrom | chr9 |
| pos (GRCh38, 1-based) | 127661125 |
| ref / alt | G / T |
| gene | STXBP1 (ENSG00000136854) |
| consequence | stop_gained (nonsense) |
| alt reads | 18 (fwd 14 / rev 4) |
| total reads (depth) | 93 |
| **allele fraction** | **0.1935** (0-1 scale) |
| mean alt base quality | 33.0 |
| alt reads near read ends (<=3 bp) | 0 |

## Reference and annotation versions

* Genome build: **GRCh38** (Broad GATK resource bundle, `Homo_sapiens_assembly38`), chromosome chr9 only, 138394717 bp. Coordinates are GRCh38, 1-based.
* Annotation: **GENCODE v47** (primary assembly), chr9 records, unmodified extract (2998 protein-coding transcripts with CDS used).
* Aligner: **minimap2 2.31-r1302**, preset `-x sr` (single-end short reads); primary alignments with MAPQ >= 20 used.
* LoF-constraint evidence: **gnomAD v4 (GRCh38)** constraint metrics via the gnomAD GraphQL API, accessed 2026-08-27 (gnomAD GraphQL API (accessed 2026-08-27), cached in work/constraint_raw.json).

## Methods

1. Reads (431396 single-end reads) were aligned to GRCh38 chr9 with minimap2 `-ax sr` (427460 primary alignments; 417753 used after MAPQ >= 20).
2. A strand-aware pileup counted bases per reference position (baseQ >= 20). Mosaic SNV candidates required: depth >= 10, >= 4 alt reads, AF in [0.02, 0.50], alt on both strands, mean alt baseQ >= 25, and no read-end artifact enrichment.
3. Candidates were annotated codon-by-codon against all GENCODE v47 protein-coding transcripts (ref codon -> alt codon). GTF reading-frame fields were cross-checked against the reconstructed CDS phase: all CDS-complete transcripts are fully consistent (GRCh38 coordinates match GENCODE v47); the only mismatches occur in 5'-truncated CDS annotations that lack an annotated start codon.
4. Stop-gained candidates were filtered for a mosaic AF (<= 0.45) and for genes that are highly LoF-intolerant per gnomAD (LOEUF < 0.35 and pLI >= 0.9).
5. The selected site was independently recounted with a separate parser (see evidence.json `independent_recount`) and checked against dbSNP.

## Results

Pileup covered 15991973 positions; 328 mosaic SNV candidates passed filters; 2 were stop-gained.

* `chr9:122477222 G>T` (OR1J1): **REJECTED** — AF 0.462 above mosaic window 0.45 (germline-het-like); no hit gene is highly LoF-intolerant (LOEUF<0.35 & pLI>=0.9)
* `chr9:127661125 G>T` (STXBP1): **SELECTED-ELIGIBLE**

### Evidence for mosaicism

* Allele fraction 0.1935 (alt 18 / depth 93) is far below the 0.5 expected for a heterozygous germline variant and above sequencing-error background; consistent with somatic mosaicism.
* Alt allele supported by 18 independent reads on both strands (fwd 14, rev 4), all with MAPQ 60 and baseQ 33; alt bases distributed across read positions (no terminal artifact).
* Independent recount: 94 overlapping primary reads; filtered depth 93; alt 18; AF 0.1935 (matches pileup).
* No dbSNP record at chr9:127661125, consistent with a rare mosaic event rather than a common polymorphism.

### Evidence for loss-of-function intolerance

* STXBP1 (ENSG00000136854): gnomAD v4 observed LoF variants = 3 vs expected = 78.5 (oe_lof = 0.038); **LOEUF (oe_lof_upper) = 0.0987**; **pLI = 1.0000**; lof_z = 7.23. This places STXBP1 among the most LoF-constrained human genes.
* Heterozygous LoF variants in STXBP1 cause developmental and epileptic encephalopathy 4 (DEE4; OMIM 612164; gene entry OMIM 602926; NCBI Gene ID 6812), an autosomal-dominant disorder due to haploinsufficiency - an established dominant LoF disease mechanism matching the 'highly LoF-intolerant' criterion. NCBI Gene places STXBP1 on NC_000009.12 (GRCh38 chr9: 127,611,911-127,696,028), spanning the called position.

### Consequence details

* Representative transcript: ENST00000373299.5 [MANE Select] (gene strand +); CDS position c.349, protein position 117, E (Glu) -> stop.
* Codon change: GAA -> TAA at CDS position 349 (variant position within codon: 0, 0-based; 0 = first codon base).
* Stop-gained in 17 STXBP1 transcripts total; alternative transcripts report c.307 p.Glu103* or c.349 p.Glu117* depending on 5' UTR/CDS start.

### Rejected alternative

* `chr9:122477222 G>T` in OR1J1: AF 0.462 above mosaic window 0.45 (germline-het-like); no hit gene is highly LoF-intolerant (LOEUF<0.35 & pLI>=0.9)
* OR1J1 is an olfactory receptor gene (no gnomAD LoF constraint data, flagged `no_exp_lof`); OR genes are LoF-tolerant and prone to multi-mapping, and AF 0.462 resembles a heterozygous germline variant rather than mosaicism.

## Deliverables

* `output/variant.tsv` — called variant (chrom, pos, ref, alt, gene, consequence, alt_reads, total_reads, allele_fraction)
* `output/evidence.json` — machine-readable evidence bundle
* `output/analysis.py` — this reproducible pipeline
* `output/report.md` — this report
* Intermediates: `work/aln/aln.sam`, `work/candidates.tsv`, `work/candidates.json`, `work/constraint_raw.json`
