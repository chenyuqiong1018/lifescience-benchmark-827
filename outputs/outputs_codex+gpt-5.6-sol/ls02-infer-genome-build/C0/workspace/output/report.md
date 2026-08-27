# Genome-build inference

## Call

The VCF uses **hg19** coordinates, with **high confidence**. All 84,664 VCF REF alleles match the hg19 chr20 reference at their declared 1-based coordinates; the selected call has 0 mismatch(es).

## Reproducible allele checks

| Candidate | chr20 length | REF matches | REF mismatches | Match fraction |
|---|---:|---:|---:|---:|
| hg18 | 62,435,964 | 20,758 | 63,906 | 0.245181 |
| hg19 | 63,025,520 | 84,664 | 0 | 1.000000 |
| hg38 | 64,444,167 | 20,966 | 63,698 | 0.247638 |

The code checks the complete VCF, including multibase REF alleles, using `sequence[POS-1:POS-1+len(REF)]`. Chromosome naming (`20`) is recorded only as QC and is not used as build proof. The winner exceeds the next-best candidate by 63,698 REF matches.

## T2T limitation

T2T/hs1 cannot be tested because no corresponding reference file is supplied. It is explicitly marked unavailable, not assigned a mismatch count. The three supplied references nevertheless yield a unique dominant call.

## Reproduction

Run `python analysis.py` in this directory. Input/reference SHA-256 values, mismatch examples, coordinate convention, and every candidate's counts are stored in `build_call.json`.
