# Controlled T1 genome-build call

## Result

The chr20 VCF uses **hg19** coordinates with **high confidence**. All 84,664 REF alleles match hg19; mismatches: 0.

| Candidate | Matches | Mismatches | Out of range | Fraction |
|---|---:|---:|---:|---:|
| hg18 | 20,758 | 63,906 | 606 | 0.245181 |
| hg19 | 84,664 | 0 | 0 | 1.000000 |
| hg38 | 20,966 | 63,698 | 0 | 0.247638 |

The analysis exhaustively checked `FASTA[POS-1:POS-1+len(REF)]`; chromosome naming was recorded only as QC. T2T/hs1 is explicitly excluded because no T2T reference was supplied, not scored as a mismatch. The controlled chromosome/UCSC skills informed assembly-aware coordinate validation and the code-execution skill informed the reproducible implementation.
