# Primer pair audit against supplied transcript isoforms

Generated: 2026-08-27 19:39 by `output/analysis.py`.

**Inputs**

- `inputs/transcripts.fa` - 2 transcript record(s)
- `inputs/primer_candidates.csv` - 3 primer pair(s)

## 1. Scope and method

- Only the supplied sequences were used; no external data were fetched or inferred.
- Primer binding was assessed by **exact match**: the forward primer must occur on the
  sense strand, and the reverse primer must occur (as its reverse complement) on the
  sense strand downstream of the forward primer, without overlap.
- Amplicon coordinates are 1-based inclusive, spanning the forward primer start to the
  end of the reverse-primer binding site (both primers included in the length).
- A transcript counts as *matched* only when a productive amplicon forms on it.
- `cds_compatible = yes` only when the amplicon lies fully inside the annotated CDS
  range **and** that CDS annotation is internally consistent.
- Malformed or internally inconsistent metadata are reported below and were **not repaired**.

## 2. Transcript inventory and metadata validation

| transcript | length (nt) | header (verbatim) | CDS annotation | metadata verdict |
|---|---|---|---|---|
| TX_CANONICAL | 102 | `>TX_CANONICAL exon_joined CDS=101-700` | 101-700 | INCONSISTENT: CDS end 700 exceeds sequence length 102 (annotation claims 700 nt, only 102 nt supplied) |
| TX_ALT | 102 | `>TX_ALT exon_joined CDS=101-640` | 101-640 | INCONSISTENT: CDS end 640 exceeds sequence length 102 (annotation claims 640 nt, only 102 nt supplied) |

**Metadata findings (reported, not repaired):**

- `TX_CANONICAL`: header states `CDS=101-700` but the supplied sequence is only 102 nt long; the CDS end exceeds the sequence length by 598 nt. No CDS-based conclusion for this transcript can be trusted; the annotation was left as supplied.
- `TX_ALT`: header states `CDS=101-640` but the supplied sequence is only 102 nt long; the CDS end exceeds the sequence length by 538 nt. No CDS-based conclusion for this transcript can be trusted; the annotation was left as supplied.

Both headers carry the tag `exon_joined`; no exon/intron structure was supplied, so no
junction-specific checks are possible.

## 3. Primer audit results

| pair_id | transcripts_matched | amplicon_length | cds_compatible | status |
|---|---|---|---|---|
| p01 | TX_CANONICAL | 102 | no | length_mismatch |
| p02 | TX_ALT | 99 | no | length_mismatch |
| p03 | none | NA | no | no_amplicon |

### Per-pair details

#### p01

- Forward (5'->3'): `ATGGCTGACCTGACTGACCT` (20 nt)
- Reverse (5'->3'): `TCATTCAGGTCAGTCAGGTC` (20 nt)
- Expected transcript: `TX_CANONICAL`; expected product: 108 bp
- TX_CANONICAL: forward exact match at 1; reverse primer (as revcomp GACCTGACTGACCTGAATGA) exact match at 83; amplicon 1-102 (102 bp)
- TX_ALT: no productive amplicon (forward binds exactly at 1; reverse primer (revcomp) has no exact match (best near-match 2/20 mismatches at pos 83))
- expected product 108 bp exceeds TX_CANONICAL sequence length 102 nt (infeasible with supplied sequence)
- TX_CANONICAL: amplicon 1-102 is NOT contained in stated CDS 101-700
- TX_CANONICAL: CDS annotation is internally inconsistent (CDS end 700 exceeds sequence length 102 (annotation claims 700 nt, only 102 nt supplied)); reported, not repaired
- observed amplicon 102 bp differs from expected 108 bp (delta -6 bp)

#### p02

- Forward (5'->3'): `GCTGACCTGACTGACCTGAA` (20 nt)
- Reverse (5'->3'): `CTATTCAGGTCAGTCAGGTC` (20 nt)
- Expected transcript: `TX_ALT`; expected product: 104 bp
- TX_CANONICAL: no productive amplicon (forward binds exactly at 4; reverse primer (revcomp) has no exact match (best near-match 2/20 mismatches at pos 7))
- TX_ALT: forward exact match at 4; reverse primer (as revcomp GACCTGACTGACCTGAATAG) exact match at 83; amplicon 4-102 (99 bp)
- expected product 104 bp exceeds TX_ALT sequence length 102 nt (infeasible with supplied sequence)
- TX_ALT: amplicon 4-102 is NOT contained in stated CDS 101-640
- TX_ALT: CDS annotation is internally inconsistent (CDS end 640 exceeds sequence length 102 (annotation claims 640 nt, only 102 nt supplied)); reported, not repaired
- observed amplicon 99 bp differs from expected 104 bp (delta -5 bp)

#### p03

- Forward (5'->3'): `AAAAAAAAAAAAAAAAAAAA` (20 nt)
- Reverse (5'->3'): `TTTTTTTTTTTTTTTTTTTT` (20 nt)
- Expected transcript: `TX_CANONICAL`; expected product: 120 bp
- TX_CANONICAL: no productive amplicon (forward has no exact match (best near-match 14/20 mismatches at pos 8); reverse primer (revcomp) has no exact match (best near-match 14/20 mismatches at pos 8))
- TX_ALT: no productive amplicon (forward has no exact match (best near-match 14/20 mismatches at pos 8); reverse primer (revcomp) has no exact match (best near-match 14/20 mismatches at pos 8))
- expected product 120 bp exceeds TX_CANONICAL sequence length 102 nt (infeasible with supplied sequence)

## 4. Cross-file consistency findings

- `p01`: expected product 108 bp is longer than the entire `TX_CANONICAL` transcript (102 nt); the stated expectation cannot be satisfied by the supplied sequences.
- `p02`: expected product 104 bp is longer than the entire `TX_ALT` transcript (102 nt); the stated expectation cannot be satisfied by the supplied sequences.
- `p03`: expected product 120 bp is longer than the entire `TX_CANONICAL` transcript (102 nt); the stated expectation cannot be satisfied by the supplied sequences.
- `p03` uses homopolymer primers (20x dA / 20x dT). No 20-mer A/T homopolymer run exists in either
  supplied transcript, so the pair cannot prime anywhere; this is flagged as a design defect, not
  something to be repaired.

## 5. Warnings raised during parsing

- transcripts.fa TX_CANONICAL: internally inconsistent CDS annotation 'CDS=101-700' vs sequence length 102: CDS end 700 exceeds sequence length 102 (annotation claims 700 nt, only 102 nt supplied). Reported as-is; NOT repaired.
- transcripts.fa TX_ALT: internally inconsistent CDS annotation 'CDS=101-640' vs sequence length 102: CDS end 640 exceeds sequence length 102 (annotation claims 640 nt, only 102 nt supplied). Reported as-is; NOT repaired.

## 6. Conclusions

- 0 of 3 primer pair(s) fully pass (amplicon on the expected transcript, expected size, CDS-compatible).
- 2 pair(s) amplify the expected transcript but produce an amplicon whose length differs from the expected product size.
- 1 pair(s) produce no amplicon on any supplied transcript.
- No primer pair is CDS-compatible under the supplied (unrepaired) annotations.
- Both transcript headers carry CDS ranges that exceed their sequence lengths; these internally
  inconsistent annotations were reported and deliberately left unrepaired.

