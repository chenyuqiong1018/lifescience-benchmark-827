# Genome-build inference for chr20 VCF

**Call: `hg19` (confidence: high).**

## Input

- VCF: `inputs/vcf.infer.build.q1.vcf.gz` - 84664 variant records, all on
  chromosome 20 (CHROM label '20'), POS 61795-62965185,
  0 multi-base REF alleles, no `##contig` header lines.
- References: chr20 FASTAs for hg18, hg19, hg38 under `inputs/references/`.
  SHA-256 checksums of all three files match `reference_manifest.json`
  (verified), so the references are authentic.

## Method

For every one of the 84664 VCF records the complete REF allele was extracted
from each supplied chr20 reference at the record's declared 1-based POS and
compared exactly. Each variant is scored per build as match, mismatch, or
out of range (coordinate/allele extends past the end of that reference's
chr20). Chromosome naming was deliberately ignored as evidence; only
allele/coordinate agreement decides the call.

## Results

| Build | chr20 length | REF matches | REF mismatches | Out of range | Match rate (in range) |
|---|---|---|---|---|---|
| hg18 | 62435964 | 20758 | 63300 | 606 | 24.6949% |
| hg19 | 63025520 | 84664 | 0 | 0 | 100.0000% |  <- called
| hg38 | 64444167 | 20966 | 63698 | 0 | 24.7638% |

- `hg19` reproduces **84664/84664** REF alleles exactly, with
  0 mismatches and 0 out-of-range positions.
- `hg38` (next best) matches only
  20966/84664 in-range alleles (24.7638%),
  and every other build is excluded by
  127604 discriminative variants whose REF allele matches `hg19`
  at the declared coordinates but not the alternative build
  (63906 vs hg18, 63698 vs hg38).
- Variants with POS/REF beyond a reference's chr20 length are structurally
  impossible in that build (e.g. hg18 chr20 is 62435964 bp
  while VCF positions reach 62965185).

## T2T

T2T/hs1: no T2T chr20 reference file is supplied with this task; the reference_manifest.json records that the authoritative chromosome-only UCSC T2T endpoint does not exist. Per inputs/references/README.md, T2T is excluded explicitly and its absence is not counted as a mismatch. The call below is made solely on the supplied references.

Because `hg19` matches 100% of REF alleles at the declared coordinates, a
T2T coordinate interpretation is also empirically excluded: T2T chr20
(~66,210,255 bp) differs from the GRCh37/hg19 arrangement by
structural changes, which would break REF agreement for variants near the
rearranged regions.

## Reproducibility

`python output/analysis.py` re-runs the entire analysis (stdlib only) and
regenerates `output/build_call.json` and this report.
