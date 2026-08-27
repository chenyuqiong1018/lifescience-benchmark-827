# Primer/transcript audit

The controlled T1 skill set was freshly installed and opened in the prescribed order. It guided local DNA validation, transcript/isoform auditing, and reproducible computation. Ensembl, UCSC, and other remote tools were not called because the prompt restricts analysis to the supplied sequences.

- **p01:** `TX_CANONICAL` is the only matched transcript. The computed product is 102 bp, not the declared 108 bp. Its `CDS=101-700` annotation lies outside the 102-nt sequence. Status: fail; CDS compatible: false.
- **p02:** `TX_ALT` is the only matched transcript. The computed product is 99 bp, not the declared 104 bp. Its `CDS=101-640` annotation lies outside the 102-nt sequence. Status: fail; CDS compatible: false.
- **p03:** no valid amplicon is present in either supplied transcript. Its expected transcript also has out-of-range CDS metadata. Status: fail; CDS compatible: not applicable.

All malformed or inconsistent metadata is retained and reported; none is silently repaired.
