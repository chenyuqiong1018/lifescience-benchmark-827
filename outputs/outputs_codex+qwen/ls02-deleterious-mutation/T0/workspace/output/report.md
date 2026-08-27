# Mosaic nonsense SNV in a LoF-intolerant gene on chr9 — analysis report

## Result

| Field | Value |
|---|---|
| chrom | chr9 |
| pos | 127,661,125 (GRCh38, 1-based) |
| ref / alt | G / T |
| gene | **STXBP1** (ENSG00000136854.26) |
| consequence | **stop_gained** (nonsense) |
| protein effect | p.(Glu117Ter) — codon GAA (Glu) → TAA (stop), CDS position 349 |
| reported transcript | ENST00000373299.5 (MANE Select; identical call in all 17 STXBP1 protein-coding transcripts, p.Glu117Ter / p.Glu103Ter) |
| alt reads / total reads | 18 / 93 (Q>=20 bases) |
| **allele fraction** | **0.1935** (reported on a 0-1 scale) |
| strand support | 14 forward / 4 reverse (two-sided binomial p = 0.031) |
| alt-read cleanliness | 18/18 alt reads carry no other candidate allele (clean_frac = 1.0) |
| mean alt base quality | 33.0 |

The data contain exactly one high-confidence mosaic nonsense SNV in a highly
loss-of-function (LoF)-intolerant protein-coding gene: a post-zygotic
STXBP1 p.Glu117Ter variant present at ~19% allele fraction.

## Reference and annotation versions

| Item | Version |
|---|---|
| Genome assembly | GRCh38 primary assembly, chromosome 9 (file `inputs/reference/GRCh38_chr9.fa.gz`, chromosome 9 from the Broad Institute GATK GRCh38 resource bundle; 138,394,717 bp) |
| Coordinate system | chromosome name `chr9`, 1-based |
| Gene annotation | GENCODE v47 (Ensembl 113, 2024-07-19), chr9 records unmodified (`inputs/reference/gencode.v47.chr9.annotation.gtf.gz`); only `protein_coding` genes/transcripts were used for the target model and consequence calling |
| Constraint prior | gnomAD v2.1.1 LoF-constraint metrics (Karczewski et al., Nature 2020), approximate published values used as external prior knowledge for gene interpretation only (no external files were read) |

## Methods (implemented in `output/analysis.py`, pure Python + numpy)

No aligner binaries were available, so the pipeline implements its own
mapper and caller:

1. **Target model.** All CDS intervals of GENCODE v47 protein-coding
   transcripts on chr9 were merged (1,423,418 bp in 8,303 intervals) and
   flanked by +/-100 bp to form the alignment target (3,016,279 bp).
2. **Seeding.** A 16-mer index (step 1, k-mers with >256 hits masked as
   repetitive) was built over the target. Each read is scanned in both
   orientations with seeds every 3 bases; genomic diagonals are voted on.
3. **Placement.** The best diagonal needs >=5 seed votes and must dominate
   any competing diagonal (>=0.8x votes at a distinct location marks the
   read as a multimapper and drops it). The placement is verified ungapped;
   the largest sub-window with <=5% mismatches (>=50 bp) is kept as the
   aligned segment.
4. **Pileup / calling.** Bases with quality >=20 (Phred+33) overlapping CDS
   positions are counted per base and strand. A position with depth >=8 and
   >=3 alt reads (AF 0.02-0.98) becomes an SNV candidate (968 candidates).
5. **Consequence.** Each candidate is reconstructed in every overlapping
   protein-coding transcript (frame-aware codon reconstruction from the
   genome) and classified (stop_gained, etc.). 26 candidates are
   stop-gained.
6. **Artifact suppression by allele co-occurrence.** The reads are scanned
   a second time: for every candidate we count the reads carrying the alt
   allele and how many of those reads simultaneously carry any other
   candidate allele. True mosaic alleles sit on reads that otherwise match
   the reference (clean_frac -> 1); paralog/pseudogene mismapping produces
   reads carrying several candidate alleles at once (clean_frac -> 0).
7. **High-confidence mosaic selection.** stop_gained AND >=5 alt reads AND
   clean_frac >=0.8 AND AF in [0.05, 0.45] AND alt reads on both strands
   AND mean alt base quality >=25. Exactly one candidate passes.

## Global statistics

- Reads: 431,396 single-end (<=151 bp); 266,881 (61.9%) uniquely placed;
  the remainder were rejected (no seed support, multimapping, or <50 bp
  aligning at <=5% mismatches — e.g. non-captured/UTR-only/repetitive
  fragments).
- Coverage of the protein-coding CDS: mean 15.1x (capture is highly
  non-uniform); 41.6% of CDS positions at >=10x; STXBP1 mean CDS depth
  62.0x, 93x at the variant position.

## Rejected alternative stop-gained candidates

| Candidate | AF | Depth | clean_frac | Rejection reason |
|---|---|---|---|---|
| FAM78A chr9:131263558 A>G | 0.331 | 254 | 0.0 | paralog decoy: ~50 co-occurring candidate alleles across chr9:131,263,486-131,263,597 at 600-800x depth; every alt read carries ~2.4 other alleles |
| C5 chr9:121017727 G>A p.Tyr550Ter | 0.525 | 59 | 1.0 | AF in germline-heterozygous range (not mosaic); C5 is recessive, not LoF-intolerant |
| ALAD chr9:113391620 A>G p.Tyr65Ter | 0.571 | 35 | 1.0 | AF > 0.5 and strong strand bias (2+/18-); ALAD recessive, not LoF-intolerant |
| GRIN3A chr9:101670953 C>T p.Gly487Ter | 0.552 | 29 | 1.0 | AF in germline-heterozygous range (not mosaic) |
| ENG chr9:127811261 G>A | 0.119 | 168 | 0.35 | alt reads co-occur with other candidate alleles (mismapping); ENG itself is LoF-intolerant but the signal is not clean |
| COL27A1 chr9:114196007 C>T p.Arg707Ter | 0.143 | 28 | 1.0 | only 4 alt reads, all on one strand; gene not highly LoF-intolerant |
| RALGPS1, OR13D1, PAEP, LAMC3 (x2), ANKRD18A (x3), ANKRD18B, GLT6D1, BAG1, C9orf85, FAM78A (x2), C5 (x2), PTPA, SEC16A, GLDC | 0.02-0.48 | 11-783 | <=0.5 | low-level artifacts: alt reads co-occur with other candidate alleles, too few reads, or single strand |

The dominant decoy in this dataset is the FAM78A region, where a dense
cluster of candidate SNVs at extreme, spike-like depth (up to ~800x) and
near-zero clean_frac marks reads from a diverged paralogous copy
mismapped onto the locus. The co-occurrence filter removes all of these.

## Evidence supporting the selected variant

- Depth 93x (56 fwd / 37 rev); 18 alt reads (14 fwd / 4 rev), all with
  base quality 33 at the variant position.
- 18/18 alt reads are "clean": none carries any of the other 967 candidate
  alleles, and alt reads average only 1.5 mismatches across ~150 bp (the
  variant itself plus background), i.e. they align uniquely and perfectly
  to STXBP1 outside the variant. Mean seed support 37.8 votes.
- Example alignment (READ_6572, genomic 127,661,070-127,661,217): the
  read differs from GRCh38 at exactly one base, the G->T at 127,661,125
  (context ...TTTAAT**G/T**AACT...; GAA -> TAA).
- A whole-chromosome scan found no second chr9 locus matching the
  surrounding 120 bp window within 10 mismatches (no local paralog).
- Consequence is a canonical nonsense substitution (Glu -> TAA stop) early
  in the open reading frame (codon 117 of 594 in the MANE Select
  transcript), guaranteeing nonsense-mediated decay / truncation —
  a true loss-of-function allele.

## Why STXBP1 is "highly LoF-intolerant"

STXBP1 (syntaxin-binding protein 1, chr9q34.11) is one of the most
LoF-constrained human genes: gnomAD v2.1.1 pLI = 1.0 and LOEUF ~0.05
(essentially depleted of protein-truncating variation in the general
population). Heterozygous loss-of-function of STXBP1 is an established
disease mechanism (STXBP1 developmental and epileptic encephalopathy,
OMIM 612164). Among the genes carrying stop-gained candidates in this
dataset, STXBP1 is the only one that is both highly LoF-intolerant and
supported by a clean, mosaic (AF < 0.5) signal.

## Mosaicism interpretation

The variant allele fraction of 0.1935 is far below the 0.5 expected for a
constitutional heterozygous variant and far above the ~1% sequencing-error
floor; all supporting reads are high quality and mapping-clean. This is
consistent with a post-zygotic (somatic) mosaic nonsense mutation present
in roughly 39% of the sampled cells (AF ~0.19 ~ heterozygous in ~39% of
cells), i.e. a high-confidence mosaic call.

## Reproducibility

```powershell
# from the workspace root (Python >= 3.10 with numpy; ~6 min)
python output/analysis.py
```

regenerates `output/variant.tsv` and `output/evidence.json` deterministically
from `inputs/` only. `evidence.json` contains the full parameter set,
read/coverage statistics, all 26 stop-gained candidates with their
co-occurrence metrics, and the read-level evidence for the final call.

## Limitations

- The mapper is ungapped: reads whose only placement would require an
  indel are dropped (acceptable for SNV calling; biases against indel-rich
  regions).
- Coverage is non-uniform (capture-like); 58% of CDS positions are below
  10x, so low-level mosaicism (<~5% AF) elsewhere cannot be excluded.
- LoF-intolerance is assessed from published gnomAD constraint metrics
  recalled as prior knowledge (no external data files were accessed), as
  the input bundle contains no constraint table.
