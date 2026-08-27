# Genome build inference report

**Input:** `inputs/vcf.infer.build.q1.vcf.gz` (84664 chr20 variant records, POS 61795–62965185; 84664 SNVs, 0 indels/other).

## Call: **hg19** (confidence: high)

The VCF REF allele of every record was compared against each supplied chr20 reference at the declared 1-based VCF coordinate (full-string comparison, so indels are verified as well). Chromosome naming was recorded but deliberately not used as evidence.

| Build | Contig length (bp) | REF matches | Mismatches | Unverifiable | Match rate |
|---|---|---|---|---|---|
| hg18 | 62,435,964 | 20,758 | 62,951 | 955 | 24.5181% |
| hg19 | 63,025,520 | 84,664 | 0 | 0 | 100.0000% |
| hg38 | 64,444,167 | 20,966 | 63,517 | 181 | 24.7638% |
| T2T (hs1) | — | — | — | — | excluded (no reference supplied) |

## Interpretation

- All 84664 REF alleles match **hg19** exactly; no other build comes close (hg38 matches only 24.76% of sites).
- A VCF aligned to a different build would show widespread REF mismatches at these same coordinates, so the coordinate system is unambiguous.

Example REF mismatches under the wrong builds (VCF REF vs reference base at the same coordinate):

- **hg18**: pos 61795: VCF=G vs hg18=C; pos 63360: VCF=C vs hg18=A; pos 63799: VCF=C vs hg18=A
- **hg38**: pos 61795: VCF=G vs hg38=A; pos 65288: VCF=G vs hg38=T; pos 65900: VCF=G vs hg38=C

## T2T/hs1

No T2T chr20 FASTA is present in `inputs/references` (its README records that the chromosome-only endpoint was unavailable). Following the supplied instructions, T2T is explicitly excluded as a candidate rather than counted as a mismatch. The exact genome-wide REF agreement with hg19 further argues against T2T coordinates.

## Reproducibility

- Reference integrity: SHA-256 of all three FASTAs matches `reference_manifest.json`.
- Re-run: `python output/analysis.py` regenerates `output/build_call.json` and this report.
- Only files under `./inputs` were read; only `./output` was written.
