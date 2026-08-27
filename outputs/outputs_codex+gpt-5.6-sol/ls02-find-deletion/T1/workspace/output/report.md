# Controlled T1 large-deletion analysis

## Call

The controlled T1 workflow supports a **1,000,000 bp chr22 deletion** between the required 100-kb-rounded boundaries **20,000,000** and **21,000,000**.

## Evidence

- 10 consecutive mostly callable 100-kb bins have zero unique-read depth; adjacent bins contain 256,104 and 261,150 mapped bases.
- 2 correctly oriented FR pairs span the interval, and their deletion-adjusted spans agree with the 500 bp normal-library median.
- 4 otherwise-unmapped reads match exactly across the inferred reference junction.

The chromosome/UCSC skills informed explicit hg38 chr22 validation and coordinate context; genome annotation informed the 1-based interpretation; code execution informed the auditable standard-library pipeline. The supplied, hash-recorded chr22 sequence remained authoritative, so unavailable credentialed SCP examples were not needed.

## Precision

The deliverable is intentionally limited to 100-kb precision. Zero-depth bins do not independently justify finer breakpoints. Exact junction reads are consistent with the displayed boundaries; in 1-based interval terms the removed sequence is approximately 20,000,001 through 21,000,000.
