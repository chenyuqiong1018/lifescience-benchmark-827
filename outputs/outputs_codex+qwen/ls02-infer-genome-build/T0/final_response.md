**Call: the VCF is on hg19 (GRCh37) coordinates, confidence: high.**

## How it was determined
`output/analysis.py` (stdlib-only, rerunnable via `python output/analysis.py`) compared **every full REF allele** of all 84,664 chr20 records against each supplied chr20 FASTA at the declared 1-based POS:

| Build | REF matches | Mismatches | Out of range | Match rate |
|---|---|---|---|---|
| **hg19** | **84,664** | **0** | 0 | **100.0000%** |
| hg38 | 20,966 | 63,698 | 0 | 24.76% |
| hg18 | 20,758 | 63,300 | 606 | 24.69% |

## Key evidence
- **hg19 reproduces 100% of REF alleles exactly** across all 84,664 variants (POS 61,795–62,965,185); the wrong builds fall to ~25% (≈ chance base agreement), with 63,698/63,906 discriminative variants refuting hg38/hg18 respectively.
- **hg18 is structurally impossible**: 606 variants lie beyond hg18 chr20's length (62,435,964 bp vs. VCF positions reaching 62,965,185).
- **T2T excluded explicitly, not as mismatches**: no T2T reference is supplied (manifest records the UCSC T2T chr20 endpoint as nonexistent); a 100% hg19 REF match is also empirically incompatible with T2T coordinates.
- **Naming was not used as proof**: the VCF has bare `20` labels and no `##contig` lines; only allele/coordinate agreement drove the call.
- **Reference integrity verified**: SHA-256 of all three FASTAs matches `reference_manifest.json`, and loaded sequence lengths equal the public chr20 lengths — ruling out mislabeled references.
- An independent manual spot-check of 7 loci (including the two tail positions beyond hg18's end) confirmed the pipeline.

## Deliverables (all in `output/`)
- `build_call.json` — build, confidence, n_variants_checked=84664, n_ref_matches=84664, n_ref_mismatches=0, evidence, plus per-build breakdown and reference-integrity details
- `analysis.py` — the reproducible analysis script
- `report.md` — full methods, results table, T2T rationale, and reproduction instructions