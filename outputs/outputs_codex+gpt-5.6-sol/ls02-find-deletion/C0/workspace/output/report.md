# Large-deletion analysis

## Call

The data support a single approximately **1,000,000 bp deletion on chr22**, with breakpoints reported at the requested 100-kb precision as **20,000,000** and **21,000,000**.

## Evidence

- Depth: 10 consecutive, mostly callable 100-kb bins from 20,000,000 to 21,000,000 have zero uniquely mapped read bases. The immediately adjacent bins contain 264,504 and 271,200 mapped bases, respectively; mean unique-read depth outside the interval is 2.665x.
- Paired-end geometry: 2 correctly oriented FR pairs bridge the interval. Their reference spans become ordinary fragment spans after subtracting the 1,000,000 bp deletion; the median ordinary FR span in the library is 500 bp.
- Junction reads: 4 otherwise-unmapped reads match exactly across the sequence formed by joining the left and right inferred breakpoints.

These three signal types agree on the same event. Repetitive terminal seeds were excluded from unique mapping, and anomalous pairs not matching the inferred interval were not counted as support.

## Precision and coordinate limits

`deletion.tsv` intentionally reports 100-kb-rounded breakpoint boundaries, as requested. The depth segmentation itself cannot justify sub-bin precision. Exact junction matches are compatible with a reference join at the displayed boundaries, but the primary deliverable should still be interpreted at 100-kb resolution. In 1-based interval language, the removed reference sequence is approximately 20,000,001 through 21,000,000; this boundary convention does not change the rounded values or the 1,000,000 bp size.

## Reproducibility

`analysis.py` uses only the Python standard library. It validates pairing and FASTQ structure, indexes observed terminal 31-mers during one reference scan, verifies full-read matches before calling them unique, computes 100-kb coverage, and confirms the candidate with FR-pair and exact-junction evidence. Full counts, hashes, and supporting record coordinates are in `qc.json`.
