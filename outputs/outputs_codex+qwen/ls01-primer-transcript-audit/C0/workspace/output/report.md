# Primer audit report

Audit of every primer pair in `inputs/primer_candidates.csv` against the supplied
transcript isoforms in `inputs/transcripts.fa`. Only the supplied sequences were
used; no external references were consulted. Malformed or internally inconsistent
metadata is reported below and was **not** silently repaired.

## 1. Method summary

- Primers were aligned with **exact (zero-mismatch) matching** only.
- Forward primer matched against the transcript sense strand; reverse primer matched as its reverse complement.
- A valid amplicon requires the forward site upstream of, and non-overlapping with, the reverse site; amplicon length spans both primer sites inclusive.
- Every pair was screened against **all** supplied isoforms to detect cross-isoform / off-target binding.
- CDS compatibility was evaluated strictly against the annotated CDS coordinates as given. Because the annotations are internally inconsistent (Section 2), compatibility could not be established for any pair.

## 2. Input validation and sequence-metadata findings

| transcript | header | sequence length | metadata findings |
|---|---|---|---|
| TX_CANONICAL | `TX_CANONICAL exon_joined CDS=101-700` | 102 bp | METADATA INCONSISTENCY: annotated CDS end 700 exceeds sequence length 102 (header 'TX_CANONICAL exon_joined CDS=101-700') |
| TX_ALT | `TX_ALT exon_joined CDS=101-640` | 102 bp | METADATA INCONSISTENCY: annotated CDS end 640 exceeds sequence length 102 (header 'TX_ALT exon_joined CDS=101-640') |

**Key metadata inconsistency (reported, not repaired).** Both transcripts are 102 bp
long, yet their headers annotate `CDS=101-700` (TX_CANONICAL) and `CDS=101-640`
(TX_ALT). In each case the annotated CDS end lies far beyond the end of the supplied
sequence, and only 2 nt (positions 101-102) of the annotated CDS exist in the
sequence. Under the `exon_joined` interpretation the coordinates apply to the spliced
transcript, so the mismatch cannot be explained away by genomic coordinates.

- Observation for TX_CANONICAL: the supplied 102 bp sequence itself is a complete in-frame ORF (starts ATG at position 1, ends with stop codon TGA at positions 100-102), which contradicts the annotated CDS coordinates. This observation is reported only; it was **not** used to overwrite or re-anchor the supplied CDS annotation.
- Observation for TX_ALT: the supplied 102 bp sequence itself is a complete in-frame ORF (starts ATG at position 1, ends with stop codon TAG at positions 100-102), which contradicts the annotated CDS coordinates. This observation is reported only; it was **not** used to overwrite or re-anchor the supplied CDS annotation.

Consequences for the audit: CDS compatibility cannot be affirmed for any amplicon,
and every pair targeting these transcripts inherits the metadata inconsistency in its
verdict.

## 3. Per-pair results

| pair_id | expected transcript | expected product | transcripts with amplicon | observed amplicon length | CDS compatible | status |
|---|---|---|---|---|---|---|
| p01 | TX_CANONICAL | 108 bp | TX_CANONICAL | 102 | no | fail |
| p02 | TX_ALT | 104 bp | TX_ALT | 99 | no | fail |
| p03 | TX_CANONICAL | 120 bp | none | na | na | fail |

### p01

- forward: `ATGGCTGACCTGACTGACCT` (20 nt)
- reverse: `TCATTCAGGTCAGTCAGGTC` (20 nt)
- expected transcript: `TX_CANONICAL`; expected product: 108 bp
- findings:
  - on TX_CANONICAL: forward primer match(es) at 1; reverse primer (reverse complement) match(es) at 83
  - amplicon on TX_CANONICAL: positions 1-102, length 102 bp
  - AMPLICON LENGTH MISMATCH: observed 102 bp vs expected 108 bp
  - expected product 108 bp exceeds the supplied transcript length (102 bp); it cannot be produced from the supplied sequence
  - cross-isoform binding without product on TX_ALT: forward primer binds at 1
  - amplicon spans positions 1-102 and is NOT contained within the annotated CDS 101-700; additionally the CDS annotation itself is internally inconsistent (see metadata findings), so CDS compatibility cannot be established

### p02

- forward: `GCTGACCTGACTGACCTGAA` (20 nt)
- reverse: `CTATTCAGGTCAGTCAGGTC` (20 nt)
- expected transcript: `TX_ALT`; expected product: 104 bp
- findings:
  - on TX_ALT: forward primer match(es) at 4; reverse primer (reverse complement) match(es) at 83
  - amplicon on TX_ALT: positions 4-102, length 99 bp
  - AMPLICON LENGTH MISMATCH: observed 99 bp vs expected 104 bp
  - expected product 104 bp exceeds the supplied transcript length (102 bp); it cannot be produced from the supplied sequence
  - cross-isoform binding without product on TX_CANONICAL: forward primer binds at 4
  - amplicon spans positions 4-102 and is NOT contained within the annotated CDS 101-640; additionally the CDS annotation itself is internally inconsistent (see metadata findings), so CDS compatibility cannot be established

### p03

- forward: `AAAAAAAAAAAAAAAAAAAA` (20 nt)
- reverse: `TTTTTTTTTTTTTTTTTTTT` (20 nt)
- expected transcript: `TX_CANONICAL`; expected product: 120 bp
- findings:
  - on TX_CANONICAL: forward primer match(es) at none; reverse primer (reverse complement) match(es) at none
  - expected product 120 bp exceeds the supplied transcript length (102 bp); it cannot be produced from the supplied sequence
  - no amplicon on the expected transcript; CDS compatibility not assessable
  - NO MATCH: neither primer has an exact binding site on any supplied transcript

## 4. Summary

- **0 of 3 primer pairs pass** the audit.
- **p01** fails: the only product on TX_CANONICAL is 102 bp, not the expected 108 bp (108 bp exceeds the supplied transcript length entirely); the amplicon is not within the annotated CDS; TX_CANONICAL CDS metadata is internally inconsistent.
- **p02** fails: the only product on TX_ALT is 99 bp, not the expected 104 bp; the amplicon is not within the annotated CDS; TX_ALT CDS metadata is internally inconsistent. The forward primer also binds TX_CANONICAL (no product forms there because the reverse primer is isoform-specific).
- **p03** fails: neither primer matches any supplied transcript; no product is possible. Note also the expected 120 bp product exceeds the 102 bp transcript length.

## 5. Recommendations

1. Correct the CDS annotations in `transcripts.fa` (as supplied they cannot be reconciled with the 102 bp sequences) before any CDS-dependent interpretation.
2. Re-derive expected product sizes from the actual supplied sequences (102 bp for p01, 99 bp for p02) or supply the transcript versions on which 108/104/120 bp products are genuine.
3. Redesign p03: poly-A/poly-T primers have no binding sites in these GC-rich transcripts.
