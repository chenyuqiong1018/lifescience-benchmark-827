**Call: the VCF is on hg19 (GRCh37) — high confidence.**

**Method** (allele/coordinate-based, naming ignored as proof): for each of the 84,664 variant records, the REF allele was fetched from every supplied chr20 reference at the declared 1-based `POS` and compared case-insensitively. All three references were SHA-256 verified against `reference_manifest.json` first. An independent seeded spot-check (10 random variants, separate code path) reproduced the same verdicts.

**Results**

| Build | REF matches | Mismatches | Out of range | Match rate |
|---|---:|---:|---:|---:|
| hg18 | 20,758 | 62,951 | 606 | 24.80% |
| **hg19** | **84,664** | **0** | **0** | **100.00%** |
| hg38 | 20,966 | 63,517 | 0 | 24.82% |
| T2T | excluded — no reference supplied (not counted as a mismatch, per the references README) ||||

Key evidence:
- hg19 reproduces every REF allele exactly (100%, zero mismatches/uninformative/out-of-range; 48,200 variants match hg19 uniquely) — margin over runner-up: 75.2 points.
- hg18/hg38 match only ~25%, the chance rate for single-base alleles; mismatches start at the very first record (POS 61,795: REF `G` vs `C` in both), showing systematic coordinate misalignment.
- Extra coordinate evidence against hg18: 606 REF intervals run past the end of hg18 chr20 (62,435,964 bp) while the VCF reaches POS 62,965,185; all positions fit within hg19/hg38.
- Contig naming (`20`, no `chr` prefix) was recorded but used only descriptively, not as proof.

**Deliverables** (all in `output/`, re-running `python output/analysis.py` is deterministic):
- `output/build_call.json` — `build: hg19`, `confidence: high`, `n_variants_checked: 84664`, `n_ref_matches: 84664`, `n_ref_mismatches: 0`, full per-build evidence incl. integrity checks and T2T exclusion.
- `output/analysis.py` — stdlib-only reproducible analysis.
- `output/report.md` — full report.