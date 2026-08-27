Audit complete. All three deliverables are in `output/` (workspace-relative):

## Deliverables
- **`output/analysis.py`** — stdlib-only, deterministic audit script (exact-match primer placement, amplicon calling, CDS containment + metadata validation; regenerates both outputs below)
- **`output/primer_audit.csv`** — `pair_id,transcripts_matched,amplicon_length,cds_compatible,status,reason`
- **`output/report.md`** — full report with per-pair coordinates, near-match diagnostics, and data-quality findings

## Key results
| pair | matched | amplicon | expected | status |
|---|---|---|---|---|
| p01 | TX_CANONICAL | 102 bp (pos 1–102) | 108 bp | length_mismatch |
| p02 | TX_ALT | 99 bp (pos 4–102) | 104 bp | length_mismatch |
| p03 | none | NA | 120 bp | no_amplicon |

- **p01**: forward exact at 1–20, reverse (revcomp `GACCTGACTGACCTGAATGA`) exact at 83–102 of TX_CANONICAL only — TX_ALT's `…GAATAG` tail gives 2 mismatches, so it is isoform-specific. Product is 102 bp, not the expected 108.
- **p02**: forward exact at 4–23 on both isoforms, but the reverse primer only resolves on TX_ALT (canonical tail differs by 2 mismatches). Product is 99 bp, not 104.
- **p03**: homopolymer 20A/20T pair has no exact match anywhere (best near-match 14/20 mismatches); it cannot prime these transcripts at all.
- **cds_compatible = no** for all pairs: no amplicon is contained in the stated CDS ranges, and the CDS annotations themselves are broken.

## Metadata inconsistencies (reported, not repaired)
- `TX_CANONICAL` header claims `CDS=101-700` but the sequence is only **102 nt** (end exceeds length by 598 nt).
- `TX_ALT` claims `CDS=101-640` against the same **102 nt** sequence (exceeds by 538 nt).
- Cross-file: all three expected product sizes (108/104/120 bp) exceed the full 102 nt transcript length, so none of the stated expectations is satisfiable from the supplied sequences.

Per the task rules, these were flagged verbatim in the CSV reasons, report warnings, and parsing log rather than silently corrected.