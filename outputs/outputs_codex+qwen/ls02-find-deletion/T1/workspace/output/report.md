# Deletion call report -- shallow paired-end hg38 (chr22)

## 1. Summary

A single large **homozygous deletion of ~1 Mb** was detected on chromosome 22:

| quantity | value |
|---|---|
| chromosome | chr22 (GRCh38) |
| reported start (rounded to 100 kb) | **20,000,000** |
| reported end (rounded to 100 kb) | **21,000,000** |
| reported size | **1,000,000 bp** |
| exact breakpoints (evidence) | last retained base 20,000,000 / first retained base 21,000,001 |
| exact deleted interval | 20,000,001-21,000,000 (1,000,000 bp) |
| breakpoint evidence source | junction_reads |

The reported start/end are the task-required rounding of the two breakpoints to the
nearest 100 kb; they are **coordinate-reporting limits, not measurement uncertainty**
(Section 5).

## 2. Data and methods

* Input: `inputs\find.deletion.r1.fq.gz`, `inputs\find.deletion.r2.fq.gz` (paired by record order,
  150 bp reads, 332,123 pairs,
  ~2.7x over mappable chr22),
  reference `inputs\reference\GRCh38_chr22.fa.gz` (GRCh38 chr22).
* Alignment: Bowtie 2 2.5.5 `--very-sensitive` against chr22 only (index built from
  the supplied FASTA). Discordant alignments were retained.
* Signals used: (i) read depth in 100 kb bins normalized by mappable (non-N) bp,
  refined at 1 kb; (ii) junction/split reads located by exact k-mer/split alignment
  against the reference; (iii) discordant FR pairs whose mates map on opposite
  flanks of the candidate region.

## 3. QC highlights (see qc.json)

* 332,120/332,123
  pairs have both mates mapped; only 0
  pairs are fully unmapped.
* Library: FR orientation dominant
  (328,008 FR pairs);
  insert median 499 bp (MAD 29).
  Note: bowtie2 FLAG 0x2 "proper" counts (169,570)
  are capped by the default `-X 500` and understate true concordance; insert stats
  were recomputed directly from mate coordinates.
* Mean depth over mappable chr22: 2.68x;
  baseline per-bp depth 2.74x.
* Depth inside the called deletion: 0.046x
  = 1.7% of baseline -> the region is
  absent from the sample (homozygous/hemizygous), not merely reduced. The residual
  reads are low-MAPQ/paralogous mismaps expected when aligning whole-genome-ish
  reads to a single-chromosome reference.
* No other chromosome arm shows a half-depth (heterozygous) or zero-depth segment
  after mappability normalization.

## 4. Evidence for the deletion

### 4.1 Read depth (100 kb bins, mappability-normalized)

Consecutive 100 kb bins with ~0 coverage (baseline 2.74x):

| bin (1-based) | covered bp | mappable bp | depth |
|---|---|---|---|
| 20,000,001 - 20,100,000 | 0 | 100,000 | 0.000 |
| 20,100,001 - 20,200,000 | 887 | 100,000 | 0.009 |
| 20,200,001 - 20,300,000 | 5,032 | 100,000 | 0.050 |
| 20,300,001 - 20,400,000 | 24,462 | 100,000 | 0.245 |
| 20,400,001 - 20,500,000 | 150 | 100,000 | 0.002 |
| 20,500,001 - 20,600,000 | 0 | 100,000 | 0.000 |
| 20,600,001 - 20,700,000 | 5,841 | 100,000 | 0.058 |
| 20,700,001 - 20,800,000 | 6,134 | 100,000 | 0.061 |
| 20,800,001 - 20,900,000 | 0 | 100,000 | 0.000 |
| 20,900,001 - 21,000,000 | 7 | 100,000 | 0.000 |

At 1 kb resolution the drop is sharp: full-depth bins continue through
~20,000,000 and coverage is ~0 from ~20,000,001; on the right side
coverage resumes at ~21,000,001. This places both breakpoints to within
about 1 kb from depth alone.

### 4.2 Junction (split) reads -- base-resolution breakpoints

3 read(s) align exactly as left-flank sequence + right-flank
sequence across the junction (no mismatches at the join):

| read | last retained base (left flank) | first retained base (right flank) |
|---|---|---|
| READ_15241/2 | 20,000,000 | 21,000,001 |
| READ_152441/2 | 20,000,000 | 21,000,001 |
| READ_250584/2 | 20,000,000 | 21,000,001 |

All junction reads agree on the same join, giving the exact deleted interval
20,000,001-21,000,000.

### 4.3 Discordant spanning pairs

2 FR pair(s) map with one mate on each flank and an apparent
span of fragment_length + deletion_size:

| read | left mate end | right mate start | apparent span | implied fragment |
|---|---|---|---|---|
| READ_220881 | 19,999,795 | 21,000,059 | 1,000,563 | 563 |
| READ_256753 | 19,999,821 | 21,000,014 | 1,000,492 | 492 |

Implied fragment lengths are consistent with the library insert distribution
(median 499 bp), supporting the same breakpoints.

### 4.4 Zygosity

Inside-depth 0.046x vs baseline
2.74x rules out a heterozygous deletion
(which would show ~half depth); the deletion is called homozygous/hemizygous.

## 5. Evidence vs. precision limits

* **Evidence (measurement):** the depth edges localize each breakpoint to ~1 kb;
  the 3 junction read(s) fix the join at single-base
  resolution; the spanning pairs independently agree.
* **Reported precision (task convention):** breakpoints are rounded to the nearest
  100 kb in `output/deletion.tsv`
  (20,000,001 -> 20,000,000;
  21,000,000 -> 21,000,000), so the reported coordinates
  carry an intrinsic +/-50 kb rounding band that is unrelated to data quality.
  Here the rounding happens to land on the Mb grid, and reported size equals the
  exact size (1,000,000 bp).

## 6. Regions considered and rejected

* Recurrent long-span FR clusters with near-constant spans (e.g. mates ~6.2 Mb or
  ~2.4 Mb apart) occur in segmental-duplication-rich pericentromeric/q11 regions.
  Depth across each implied interval is normal, so these are paralogy/mismapping
  artifacts, not deletions.
* The raw-coverage dip at 10.5-13 Mb reflects 50 kb N gaps in the GRCh38 chr22
  reference; mappability-normalized depth there is normal.
* The zero-coverage 12.9-15.2 Mb block is the reference centromeric N gap, not a
  deletion (no mappable sequence exists there).

## 7. Limitations

* Shallow depth (~2.7x): only a handful of spanning/junction
  molecules are expected at breakpoints; the call rests on their exact agreement
  plus the depth profile.
* Single-chromosome reference: reads from other chromosomes are forced to map or
  drop out, inflating apparent alignment rate and creating the small residual
  coverage inside the deletion.
* Bowtie 2 does not perform spliced/split alignment, so junction reads are
  identified here by exact re-alignment of read sequences against the reference
  (Section 4.2), not from CIGAR strings.

## 8. Reproducibility

`python output/analysis.py` regenerates `work/` intermediates (if missing) and
rewrites `output/deletion.tsv`, `output/qc.json`, and this report.
