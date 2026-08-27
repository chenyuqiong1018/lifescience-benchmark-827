Audit complete. All three deliverables are in `output/`, built using only the supplied `inputs/` files.

## Deliverables
- **`output/analysis.py`** — reproducible audit script (parses FASTA + primer CSV, validates metadata, exact-matches forward primer and RC(reverse primer) on each isoform, computes amplicons, writes CSV + report)
- **`output/primer_audit.csv`** — exact requested columns: `pair_id,transcripts_matched,amplicon_length,cds_compatible,status,reason`
- **`output/report.md`** — full report with metadata-integrity check, per-pair binding tables, method, and conclusions

## Key findings

| pair | transcripts_matched | amplicon | cds_compatible | status |
|---|---|---|---|---|
| p01 | TX_CANONICAL | 102 bp (expected 108) | inconsistent_metadata | length_mismatch |
| p02 | TX_ALT | 99 bp (expected 104) | inconsistent_metadata | length_mismatch |
| p03 | none | na | na | no_binding |

- **p01**: forward binds TX_CANONICAL uniquely at 1..20, reverse (RC) at 83..102 → product 102 bp, **6 bp shorter** than the expected 108. Isoform-specific: no amplicon on TX_ALT.
- **p02**: forward binds TX_ALT at 4..23, reverse (RC) at 83..102 → product 99 bp, **5 bp shorter** than the expected 104. The single-nucleotide tail difference (…GAATAG vs …GAATGA) correctly confines the reverse primer to TX_ALT.
- **p03**: homopolymer A₂₀/T₂₀ primers have no exact binding site in either isoform → cannot amplify anything supplied; discard.

## Metadata inconsistency (reported, not repaired)
Both headers are internally inconsistent with their sequences:
- **TX_CANONICAL**: `CDS=101-700` (600 bp span) on a **102 bp** sequence
- **TX_ALT**: `CDS=101-640` (540 bp span) on a **102 bp** sequence

Per your instruction I did not clamp or repair this — the sequences were used exactly as supplied, `cds_compatible` is flagged `inconsistent_metadata` for both amplicon-bearing pairs (even taken at face value, the amplicons aren't contained in the annotated CDS), and the report notes that either the FASTA entries are truncated excerpts or the CDS coordinates are wrong — undecidable from the supplied data alone.