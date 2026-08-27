Done. I located the large deletion and produced all four deliverables in `output/`.

## Result

**chr22:20,000,000–21,000,000 (rounded 100 kb grid), a ~1,000,000 bp homozygous deletion**

| | |
|---|---|
| chrom | chr22 |
| start_100kb / end_100kb | 20,000,000 / 21,000,000 |
| size_bp | 1,000,000 |
| refined breakpoints (1-based) | 20,000,001 – 21,000,000 |

## How it was found (no external aligner — pure NumPy pipeline in `output/analysis.py`)

1. **Seeding/alignment** — built an exact 25‑mer index of the chr22 reference (2‑bit hashes, radix-sorted), seeded each read on both strands (offsets 0/62/125), and classified reads (unique/ambiguous/split/unmapped). 96.9% of the 664,246 reads placed with ≥2 consistent seeds; a strand-aware subsample check showed reads are effectively error-free (99.7% match with ≤2 mismatches).
2. **Read depth (primary evidence)** — 100 kb bins normalized by per-bin mappable sequence: **10 consecutive bins collapse to CN 0.002 vs 1.028 elsewhere** (homozygous-like), with sharp 10 kb-resolution boundaries at 20.0 Mb and 21.0 Mb. N-gap/centromere bins are excluded (null in the QC profile).
3. **Discordant pairs** — insert median 500 bp; 1 FR pair spans the region with mates flanking the depth edges (span 1,000,563 bp).
4. **Split reads** — 2 junction reads verified base-by-base against the implied junction, both giving **exactly D = 1,000,000 bp**. The junction carries a 1 bp microhomology (A), so the left breakpoint is a 2 bp window (first deleted base 20,000,001 or 20,000,002); all resolutions round to identical 100 kb values.

## Evidence vs precision limits
- Depth alone localizes breakpoints to ~10 kb and establishes interval + zygosity; the discordant pair narrows to ~fragment length; split reads give near-single-base resolution.
- The 100 kb grid in the TSV is the **mandated reporting precision**, not the measurement limit — refined coordinates, per-signal counts, CN profiles, and the microhomology window are all preserved in `qc.json`/`report.md`.

**Deliverables:** `output/deletion.tsv` (one call row with supporting signals), `output/qc.json` (alignment/QC metrics, insert stats, CN profiles, split clusters, breakpoint windows), `output/analysis.py` (self-contained, deterministic, reruns in ~35 s), `output/report.md` (narrative). All work stayed inside the workspace using only `./inputs`.