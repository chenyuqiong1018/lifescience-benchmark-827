Audit complete. All three deliverables are in `output/`, generated entirely from the supplied inputs.

## Deliverables
| File | Contents |
|---|---|
| `output/primer_audit.csv` | One row per pair with `pair_id,transcripts_matched,amplicon_length,cds_compatible,status,reason` |
| `output/analysis.py` | Self-contained audit script (stdlib only); regenerates both the CSV and the report deterministically |
| `output/report.md` | Full audit report: methods, metadata findings, per-pair detail, recommendations |

## Key findings

**Metadata inconsistency (reported, not repaired).** Both supplied transcripts are **102 bp** long, yet their headers annotate `CDS=101-700` (TX_CANONICAL) and `CDS=101-640` (TX_ALT) — the CDS end lies far beyond the sequence, and only positions 101–102 of the annotated CDS exist. Ironically, each 102 bp sequence itself is a complete in-frame ORF (ATG start, TGA/TAG stop, length divisible by 3), contradicting the annotated CDS. Per your instructions this was flagged, not silently re-anchored or clipped.

**Per-pair verdicts (0 of 3 pass):**

- **p01 → fail.** Both primers bind TX_CANONICAL uniquely (positions 1 and 83), giving a **102 bp** amplicon — not the expected 108 bp, which exceeds the transcript length entirely and is impossible from the supplied sequence. Amplicon (1–102) falls outside the annotated CDS. Forward primer also cross-binds TX_ALT (no product there, since the reverse primer is isoform-specific).
- **p02 → fail.** Unique sites on TX_ALT (positions 4 and 83) give **99 bp**, not the expected 104 bp; outside the annotated CDS; TX_ALT metadata is likewise inconsistent. Forward primer cross-binds TX_CANONICAL without forming a product.
- **p03 → fail.** The poly-A/poly-T primers have **no exact binding site on any supplied transcript** — no amplicon is possible (and the claimed 120 bp product exceeds the 102 bp transcript anyway).

## Verification performed
- `analysis.py` compiles clean; amplicon lengths cross-checked by an independent recomputation (102 / 99 / none — matched).
- CSV header and row structure validated; metadata inconsistency is explicitly flagged in the reasons.
- Re-running `analysis.py` from a clean state regenerates byte-identical CSV and report; no external sequences were used.