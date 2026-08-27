Done. The deletion is located and all four deliverables are written under `output/`.

## Result

**chr22:20,000,000–21,000,000 (GRCh38), ~1 Mb homozygous deletion** — reported (per task convention) with breakpoints rounded to the nearest 100 kb; the exact deleted interval is chr22:20,000,001–21,000,000 (1,000,000 bp), junction at 20,000,000|21,000,001.

## How it was found

1. **Alignment** — no aligner was installed, so I downloaded Bowtie 2 2.5.5 (mingw build) into `tools/`, built an index from the supplied chr22 FASTA, and aligned the 332,123 pairs (150 bp, ~2.7× over mappable chr22) with `--very-sensitive`.
2. **Read depth** — 100 kb bins normalized by mappable (non-N) bp show exactly one zero-depth stretch: ten consecutive bins across 20.0–21.0 Mb (baseline 2.74×, inside 0.05× ≈ 1.7% → homozygous, not heterozygous). At 1 kb resolution both edges are sharp.
3. **Junction (split) reads** — 3 reads (`READ_250584/2`, `READ_15241/2`, `READ_152441/2`) align exactly as left-flank + right-flank sequence, all agreeing on the single-base join 20,000,000|21,000,001 (base-resolution evidence).
4. **Discordant spanning pairs** — 2 FR pairs (`READ_220881`, `READ_256753`) with mates on opposite flanks, spans 1,000,563/1,000,492 bp ⇒ implied fragment lengths 563/492 bp, consistent with the library (median 499 bp).

## Evidence vs. precision limits (as requested)

- **Evidence**: depth edges localize breakpoints to ~1 kb; junction reads fix them to base resolution; spanning pairs independently agree.
- **Reporting limit**: the TSV rounds each breakpoint to the nearest 100 kb (intrinsic ±50 kb band, unrelated to data quality). Here rounding lands on the Mb grid, so reported size equals the exact size.
- **Rejected distractors**: long-span FR clusters elsewhere (~6.2 Mb, ~2.4 Mb, etc.) are segmental-duplication/paralogy artifacts (normal depth across each interval); the 10.5–13 Mb raw-depth dip is explained by reference N-gaps after mappability normalization; 12.9–15.2 Mb is the centromeric N-gap.

## Deliverables

| file | content |
|---|---|
| `output/deletion.tsv` | `chr22  20000000  21000000  1000000` + supporting_signals (depth, junction reads, spanning pairs) |
| `output/qc.json` | mapping/orientation/insert/coverage QC, exact vs. rounded call, per-read evidence, artifacts reviewed |
| `output/analysis.py` | self-contained reproducible pipeline (stdlib + Bowtie 2): decompress → index → align → depth/junction/span analysis → writes all outputs (`python output/analysis.py`, verified re-runnable) |
| `output/report.md` | full methods, evidence tables, zygosity, evidence-vs-precision section, limitations |

Intermediate alignment/index files live in `work/` (referenced by `analysis.py`); only files under `./inputs` were used as data sources.