## Call: **hg19** — confidence **high**

**Method (reproducible, allele-based):** For every one of the 84,664 variant records in `inputs/vcf.infer.build.q1.vcf.gz`, the VCF REF allele was compared against each supplied chr20 reference at the declared 1-based coordinate (full-string comparison, so indels are verified too). References were first SHA-256-verified against `reference_manifest.json`. Chromosome naming (`20`) was recorded but deliberately not used as evidence — although it happens to match the hg19 FASTA header convention.

| Build | Contig length | REF matches | Match rate |
|---|---|---|---|
| hg18 | 62,435,964 bp | 20,758 | 24.52% |
| **hg19** | 63,025,520 bp | **84,664 / 84,664** | **100.0000%** |
| hg38 | 64,444,167 bp | 20,966 | 24.76% |
| T2T (hs1) | — | — | explicitly excluded (no reference file supplied; per the references README, an unavailable candidate is excluded, not counted as a mismatch) |

hg19 reproduces 100% of REF alleles with 0 mismatches and 0 unverifiable positions (POS range 61,795–62,965,185 fits within the hg19 chr20 contig), while both alternative builds mismatch ~75% of sites — an unambiguous coordinate signature.

**Verification:** an independent spot-check (different code path) on 5 random records confirmed hg19 bases; two consecutive runs of the analysis produce byte-identical outputs.

**Deliverables (all in `./output`):**
- `build_call.json` — `build: hg19`, `confidence: high`, `n_variants_checked: 84664`, `n_ref_matches: 84664`, `n_ref_mismatches: 0`, evidence array, per-build stats, integrity results, T2T status
- `analysis.py` — self-contained stdlib-only script; `python output/analysis.py` regenerates both outputs
- `report.md` — full report with match table, mismatch examples, T2T exclusion rationale, and reproducibility notes

(One note: my first ad-hoc spot-check used slow string concatenation and is still finishing in the background; per the isolation rules I left it running — it doesn't affect the deliverables.)