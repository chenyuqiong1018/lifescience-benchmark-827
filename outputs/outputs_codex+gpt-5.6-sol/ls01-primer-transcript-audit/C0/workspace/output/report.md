# Primer/transcript audit

All three primer pairs were audited against both supplied transcript sequences. No sequence or metadata was repaired.

- **p01:** One valid amplicon was found in `TX_CANONICAL`, spanning 102 bp. This disagrees with the declared 108 bp product. The transcript is 102 nt long, but its header declares `CDS=101-700`, which is out of bounds. Status: fail; CDS compatibility cannot be established and is recorded as false for the matched transcript.
- **p02:** One valid amplicon was found in `TX_ALT`, spanning 99 bp. This disagrees with the declared 104 bp product. The transcript is 102 nt long, but its header declares `CDS=101-640`, which is out of bounds. Status: fail; CDS compatibility cannot be established and is recorded as false.
- **p03:** No valid amplicon was found in either transcript. Its expected transcript `TX_CANONICAL` also carries the out-of-bounds CDS declaration. Status: fail; CDS compatibility is not applicable.

The malformed internal CDS coordinates are reported as supplied rather than silently clipped or reinterpreted.
