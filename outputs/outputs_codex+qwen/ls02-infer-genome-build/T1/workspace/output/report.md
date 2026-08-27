# Genome-build inference for `vcf.infer.build.q1.vcf.gz` (chr20)

**Call: hg19 (GRCh37) — high confidence**

## Question

Determine whether the supplied chromosome-20 VCF uses hg18, hg19, hg38 or
T2T coordinates. The call must rest on reproducible allele/coordinate
checks; chromosome naming alone is not accepted as proof.

## Inputs

| File | Role |
|---|---|
| `inputs/vcf.infer.build.q1.vcf.gz` | Query VCF (gzipped, 84,664 variant records, no `##` meta lines) |
| `inputs/references/hg18_chr20.fa.gz` | UCSC hg18 chr20 (62,435,964 bp) |
| `inputs/references/hg19_chr20.fa.gz` | Broad hg19/GRCh37 v0 chr20 (63,025,520 bp) |
| `inputs/references/hg38_chr20.fa.gz` | Broad hg38/GRCh38 v0 chr20 (64,444,167 bp) |
| `inputs/references/reference_manifest.json` | Sources + SHA-256 for each reference |

All three reference FASTAs were SHA-256 verified against the manifest before
use (all `verified: true`), so the comparison targets are known-good.

## Method

For every VCF record the REF allele was fetched from each supplied reference
at the record's declared **1-based** `POS` (slice `seq[POS-1 : POS-1+len(REF)]`,
case-insensitive comparison, `N` treated as uninformative). A VCF whose
coordinates belong to a given build must reproduce the reference at
essentially every site, so the build whose reference matches the REF
alleles is the correct one. This is implemented in `output/analysis.py`
(standard library only; rerun with `python output/analysis.py`).

Chromosome naming was recorded but deliberately **not** used as evidence:
all 84,664 records use the bare contig name `20` (no `chr` prefix), which
is compatible with either convention and proves nothing by itself.
An independent random spot-check of 10 variants (separate code path,
seed = 42) reproduced the same per-build verdicts.

## Results

| Build | REF matches | REF mismatches | Uninformative | Out of range | Match rate |
|---|---:|---:|---:|---:|---:|
| hg18 | 20,758 | 62,951 | 349 | 606 | 24.80% |
| **hg19** | **84,664** | **0** | **0** | **0** | **100.00%** |
| hg38 | 20,966 | 63,517 | 181 | 0 | 24.82% |
| T2T (hs1) | — | — | — | — | excluded (see below) |

* VCF records analyzed: **84,664** (0 skipped), POS range 61,795–62,965,185.
* Against hg19, every single REF allele matches the reference exactly:
  84,664/84,664 (100%), with zero mismatches, zero uninformative sites and
  zero out-of-range positions. 48,200 variants match hg19 *uniquely*.
* Against hg18 and hg38 only ~24.8% of REF alleles match — the rate
  expected by random chance for single-base alleles — i.e. the coordinates
  are wrong for those builds. The first mismatching sites against both
  appear at the very first record (POS 61,795: VCF REF `G`, hg18/hg38 have
  `C`), showing the divergence is systematic, not edge noise.
* Additional coordinate evidence against hg18: 606 variant REF intervals
  extend beyond the end of hg18 chr20 (length 62,435,964) while the VCF
  reaches POS 62,965,185; all positions fit within hg19 (63,025,520 bp)
  and hg38 (64,444,167 bp).
* Match-rate margin of the winner over the runner-up: **0.7518**.

### T2T candidate

No T2T/hs1 chr20 FASTA is supplied under `inputs/references` (the manifest
notes the chromosome-only endpoint was unavailable). Following the
references README, T2T is **explicitly excluded** from scoring rather than
counted as a mismatch; it cannot win, and nothing here is treated as
evidence against it.

## Decision rule

The build with the highest REF match rate is called, requiring ≥99% match
rate and ≥5-point margin over the runner-up for a "high" confidence call.
hg19 satisfies this with 100% matches and a 75.2-point margin.

## Deliverables

* `output/build_call.json` — machine-readable call:
  `build=hg19`, `confidence=high`, `n_variants_checked=84664`,
  `n_ref_matches=84664`, `n_ref_mismatches=0`, plus full per-build evidence.
* `output/analysis.py` — reproducible analysis (stdlib-only).
* `output/report.md` — this report.
