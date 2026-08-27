# Primer audit report

Generated: 2026-08-27 20:06 by `output/analysis.py`

Scope: 3 primer pair(s) audited against 2 supplied transcript isoform(s). Only the supplied files `inputs/primer_candidates.csv` and `inputs/transcripts.fa` were used; no external sequences or databases were consulted. Matching is exact (no mismatches/gaps); the reverse primer is matched as its reverse complement against the sense strand.

## 1. Transcript metadata integrity check

**Internally inconsistent metadata was detected and is reported below. Per the audit policy it was NOT repaired, clamped, or ignored.**

| transcript | header | seq length (bp) | CDS annotation | metadata verdict |
|---|---|---|---|---|
| TX_CANONICAL | `TX_CANONICAL exon_joined CDS=101-700` | 102 | 101..700 (span 600 bp) | INCONSISTENT |
| TX_ALT | `TX_ALT exon_joined CDS=101-640` | 102 | 101..640 (span 540 bp) | INCONSISTENT |

- **TX_CANONICAL**: annotated CDS=101-700 extends past the end of the 102 bp sequence (annotated CDS span of 600 bp cannot fit inside a 102 bp transcript). The sequence was used exactly as supplied.
- **TX_ALT**: annotated CDS=101-640 extends past the end of the 102 bp sequence (annotated CDS span of 540 bp cannot fit inside a 102 bp transcript). The sequence was used exactly as supplied.

Consequence: for any amplicon mapping to a transcript with an inconsistent CDS annotation, `cds_compatible` is reported as `inconsistent_metadata` instead of yes/no, because the annotation cannot support a reliable compatibility call.

## 2. Method

- Forward primer: exact substring search on the sense strand of each transcript.
- Reverse primer: reverse-complemented, then exact substring search on the sense strand.
- A pair matches a transcript only if a forward site lies upstream of a reverse-complement site; amplicon = forward_start .. rev_site_end, length = end - start + 1.
- `cds_compatible` = `yes` only if the amplicon lies fully inside a self-consistent annotated CDS; `inconsistent_metadata` if the CDS annotation contradicts the sequence; `na` if there is no amplicon / no CDS annotation.
- Status vocabulary: `pass`, `length_mismatch`, `no_binding`, `off_target`, `ambiguous`, `cds_incompatible`, `metadata_conflict`, `expected_transcript_missing`, `malformed_input`.

## 3. Per-pair audit results

### p01 (expected transcript: TX_CANONICAL, expected product: 108 bp)

- forward: `ATGGCTGACCTGACTGACCT` (20 nt)
- reverse: `TCATTCAGGTCAGTCAGGTC` (20 nt), RC = `GACCTGACTGACCTGAATGA`

| transcript | forward hits | RC(reverse) hits | amplicon(s) |
|---|---|---|---|
| TX_CANONICAL | 1 | 83 | 1..102 (102 bp) |
| TX_ALT | 1 | - | - |

- transcripts_matched: **TX_CANONICAL**
- amplicon_length: **102**
- cds_compatible: **inconsistent_metadata**
- status: **length_mismatch**
- reason: forward primer binds TX_CANONICAL uniquely at 1..20. reverse primer (reverse-complemented) binds TX_CANONICAL uniquely at 83..102. predicted amplicon TX_CANONICAL:1..102, length 102 bp. amplicon length 102 bp differs from expected_product_bp 108 (delta -6 bp). cds_compatible=inconsistent_metadata: CDS annotation CDS=101-700 is internally inconsistent with the 102 bp sequence (annotated CDS=101-700 extends past the end of the 102 bp sequence (annotated CDS span of 600 bp cannot fit inside a 102 bp transcript)); compatibility cannot be assessed; metadata reported, not repaired; even taking the annotation at face value, amplicon 1..102 is not contained within CDS 101..700. no amplicon on other supplied transcript(s): TX_ALT. TX_CANONICAL header metadata issue (reported, not repaired): annotated CDS=101-700 extends past the end of the 102 bp sequence (annotated CDS span of 600 bp cannot fit inside a 102 bp transcript).

### p02 (expected transcript: TX_ALT, expected product: 104 bp)

- forward: `GCTGACCTGACTGACCTGAA` (20 nt)
- reverse: `CTATTCAGGTCAGTCAGGTC` (20 nt), RC = `GACCTGACTGACCTGAATAG`

| transcript | forward hits | RC(reverse) hits | amplicon(s) |
|---|---|---|---|
| TX_CANONICAL | 4 | - | - |
| TX_ALT | 4 | 83 | 4..102 (99 bp) |

- transcripts_matched: **TX_ALT**
- amplicon_length: **99**
- cds_compatible: **inconsistent_metadata**
- status: **length_mismatch**
- reason: forward primer binds TX_ALT uniquely at 4..23. reverse primer (reverse-complemented) binds TX_ALT uniquely at 83..102. predicted amplicon TX_ALT:4..102, length 99 bp. amplicon length 99 bp differs from expected_product_bp 104 (delta -5 bp). cds_compatible=inconsistent_metadata: CDS annotation CDS=101-640 is internally inconsistent with the 102 bp sequence (annotated CDS=101-640 extends past the end of the 102 bp sequence (annotated CDS span of 540 bp cannot fit inside a 102 bp transcript)); compatibility cannot be assessed; metadata reported, not repaired; even taking the annotation at face value, amplicon 4..102 is not contained within CDS 101..640. no amplicon on other supplied transcript(s): TX_CANONICAL. TX_ALT header metadata issue (reported, not repaired): annotated CDS=101-640 extends past the end of the 102 bp sequence (annotated CDS span of 540 bp cannot fit inside a 102 bp transcript).

### p03 (expected transcript: TX_CANONICAL, expected product: 120 bp)

- forward: `AAAAAAAAAAAAAAAAAAAA` (20 nt)
- reverse: `TTTTTTTTTTTTTTTTTTTT` (20 nt), RC = `AAAAAAAAAAAAAAAAAAAA`

| transcript | forward hits | RC(reverse) hits | amplicon(s) |
|---|---|---|---|
| TX_CANONICAL | - | - | - |
| TX_ALT | - | - | - |

- transcripts_matched: **none**
- amplicon_length: **na**
- cds_compatible: **na**
- status: **no_binding**
- reason: forward primer has no exact match in TX_CANONICAL. reverse primer (reverse-complemented) has no exact match in TX_CANONICAL. forward primer matches none of the 2 supplied transcripts. reverse-complemented reverse primer matches none of the 2 supplied transcripts. no amplicon can be formed on any supplied isoform. no amplicon on other supplied transcript(s): TX_ALT. TX_CANONICAL header metadata issue (reported, not repaired): annotated CDS=101-700 extends past the end of the 102 bp sequence (annotated CDS span of 600 bp cannot fit inside a 102 bp transcript).

## 4. Summary (contents of `output/primer_audit.csv`)

| pair_id | transcripts_matched | amplicon_length | cds_compatible | status | reason |
|---|---|---|---|---|---|
| p01 | TX_CANONICAL | 102 | inconsistent_metadata | length_mismatch | forward primer binds TX_CANONICAL uniquely at 1..20. reverse primer (reverse-complemented) binds TX_CANONICAL uniquely at 83..102. predicted amplicon TX_CANONICAL:1..102, length 102 bp. amplicon length 102 bp differs from expected_product_bp 108 (delta -6 bp). cds_compatible=inconsistent_metadata: CDS annotation CDS=101-700 is internally inconsistent with the 102 bp sequence (annotated CDS=101-700 extends past the end of the 102 bp sequence (annotated CDS span of 600 bp cannot fit inside a 102 bp transcript)); compatibility cannot be assessed; metadata reported, not repaired; even taking the annotation at face value, amplicon 1..102 is not contained within CDS 101..700. no amplicon on other supplied transcript(s): TX_ALT. TX_CANONICAL header metadata issue (reported, not repaired): annotated CDS=101-700 extends past the end of the 102 bp sequence (annotated CDS span of 600 bp cannot fit inside a 102 bp transcript). |
| p02 | TX_ALT | 99 | inconsistent_metadata | length_mismatch | forward primer binds TX_ALT uniquely at 4..23. reverse primer (reverse-complemented) binds TX_ALT uniquely at 83..102. predicted amplicon TX_ALT:4..102, length 99 bp. amplicon length 99 bp differs from expected_product_bp 104 (delta -5 bp). cds_compatible=inconsistent_metadata: CDS annotation CDS=101-640 is internally inconsistent with the 102 bp sequence (annotated CDS=101-640 extends past the end of the 102 bp sequence (annotated CDS span of 540 bp cannot fit inside a 102 bp transcript)); compatibility cannot be assessed; metadata reported, not repaired; even taking the annotation at face value, amplicon 4..102 is not contained within CDS 101..640. no amplicon on other supplied transcript(s): TX_CANONICAL. TX_ALT header metadata issue (reported, not repaired): annotated CDS=101-640 extends past the end of the 102 bp sequence (annotated CDS span of 540 bp cannot fit inside a 102 bp transcript). |
| p03 | none | na | na | no_binding | forward primer has no exact match in TX_CANONICAL. reverse primer (reverse-complemented) has no exact match in TX_CANONICAL. forward primer matches none of the 2 supplied transcripts. reverse-complemented reverse primer matches none of the 2 supplied transcripts. no amplicon can be formed on any supplied isoform. no amplicon on other supplied transcript(s): TX_ALT. TX_CANONICAL header metadata issue (reported, not repaired): annotated CDS=101-700 extends past the end of the 102 bp sequence (annotated CDS span of 600 bp cannot fit inside a 102 bp transcript). |

## 5. Conclusions

- 0 of 3 primer pairs fully pass the audit.
- **p01**: primers bind the expected transcript TX_CANONICAL uniquely, but the predicted product (102 bp) does not match expected_product_bp (108); the pair should not be used as-is.
- **p02**: primers bind the expected transcript TX_ALT uniquely, but the predicted product (99 bp) does not match expected_product_bp (104); the pair should not be used as-is.
- **p03**: neither primer has an exact binding site in any supplied transcript; the pair cannot amplify the supplied isoforms and should be discarded.
- Metadata: transcript(s) TX_CANONICAL, TX_ALT carry CDS annotations that exceed the actual sequence length. Either the FASTA sequences are truncated excerpts or the CDS coordinates are wrong; this cannot be resolved from the supplied data alone, so the inconsistency is reported rather than repaired.

