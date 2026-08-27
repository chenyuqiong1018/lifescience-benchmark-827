# Primer/transcript audit

The selected whitelist skills were used for local DNA validation, transcript-aware sequence matching, and reproducible code structure. Their remote tools were not called because the task permits only the supplied sequences.

All six primer oligonucleotides are valid 20-nt DNA strings, and both transcript sequences contain only A/C/G/T. However, both FASTA CDS declarations are internally inconsistent with the 102-nt sequences.

- **p01:** matches only `TX_CANONICAL`; observed amplicon 102 bp versus expected 108 bp. `CDS=101-700` exceeds the 102-nt transcript. Status: fail; CDS compatible: false.
- **p02:** matches only `TX_ALT`; observed amplicon 99 bp versus expected 104 bp. `CDS=101-640` exceeds the 102-nt transcript. Status: fail; CDS compatible: false.
- **p03:** no valid amplicon in either transcript. The expected `TX_CANONICAL` metadata is also out of bounds. Status: fail; CDS compatible: not applicable.

No metadata was clipped, repaired, or inferred from external sources.
