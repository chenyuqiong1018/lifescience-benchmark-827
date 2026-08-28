# T0 genome-build inference

## Call

The VCF uses **hg19** coordinates with **high confidence**: 84,664/84,664 REF alleles match and 0 mismatch.

| Candidate | REF matches | REF mismatches | Fraction |
|---|---:|---:|---:|
| hg18 | 20,758 | 63,906 | 0.245181 |
| hg19 | 84,664 | 0 | 1.000000 |
| hg38 | 20,966 | 63,698 | 0.247638 |

Every VCF REF was checked at its declared 1-based coordinate; chromosome naming was QC only. T2T is explicitly untested because no T2T reference was supplied. The selected chromosome/UCSC skills informed assembly-aware coordinate checks, while code-execution guidance informed the auditable exhaustive script.
