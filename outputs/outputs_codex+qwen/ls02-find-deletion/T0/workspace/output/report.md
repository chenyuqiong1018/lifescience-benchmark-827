# Large deletion on GRCh38 chr22 - analysis report

## Summary

A large **homozygous** deletion of **1,000,000 bp** was identified on **chr22**. Rounded to the nearest 100 kb as required, the event spans 1-based coordinates **20,000,000 to 21,000,000** (size = end - start = 1,000,000 bp). The raw breakpoint evidence places the junctions at **20,000,000** (last retained base) and **21,000,000** (last deleted base).

| field | value |
|---|---|
| chromosome | chr22 (GRCh38) |
| start_100kb (1-based) | 20,000,000 |
| end_100kb (1-based) | 21,000,000 |
| size_bp | 1,000,000 |
| zygosity estimate | homozygous (region depth ratio 0.006 of genomic median) |
| raw breakpoint method | split_reads |

## Data

- 332,123 paired-end read pairs x 2 x 150 bp (99,636,900 bp), mean Q R1=40.0, R2=40.0; pairs are matched by record order.
- Reference: chr22 (50,818,468 bp; 13,746,488 non-ACGT bases, mostly the p-arm/centromeric gaps).
- Effective depth ~2.652x (shallow). Reads were mapped with a custom numpy seed-and-vote aligner: R1 98.68% / R2 98.64% mapped, mean identity 149.9 and 149.9 matches out of 150.
- Concordant FR insert size: median 500 bp (MAD 34 bp) over 156,835 pairs.

## Methods

1. **Mapping** - repeat-masked 20-mer index of chr22 (k-mers with >16 copies treated as repetitive); seeds every 10 bp on both strands; diagonal voting; full-read verification (>=130/150 matches).
2. **Read depth** - mapped read starts counted in 100 kb and 10 kb bins, normalized by per-bin mappable read-start count (k-mer uniqueness); the deletion is the longest run of bins with ratio < 0.62x the genomic median.
3. **Discordant pairs** - FR pairs with insert above the robust threshold, both mates near-unique (>=140/150 matches); junction-spanning pairs bracket the breakpoints (left bound = max left-mate end, right bound = min right-mate start).
4. **Split reads** - reads with two strong diagonal clusters (or a partial alignment whose unmatched segment maps elsewhere) re-aligned segment-by-segment for nucleotide-resolution junctions; only split reads whose two junctions bracket the depth interval are accepted.

## Evidence for the deletion

1. **Read depth.** 10 consecutive 100 kb bins (20,000,000-21,000,000, 0-based) have median normalized ratio 0.006 versus the genomic median 0.983: a near-complete, homozygous loss. At 10 kb resolution the depletion starts at 20,000,000 and ends at 21,000,000 (0-based).
2. **Split reads.** 3 read(s) span the junction with near-perfect segment matches:
   - R1 read 139652 (strand 1): left junction 20,000,000, right junction 21,000,000; segment matches 44/44 and 106/106
   - R2 read 15240 (strand 0): left junction 20,000,000, right junction 21,000,000; segment matches 98/98 and 52/52
   - R2 read 152440 (strand 0): left junction 20,000,001, right junction 21,000,001; segment matches 100/100 and 50/50
3. **Discordant pairs.** 1 uniqueness-filtered pair(s) span the interval:
   - left mate starts at 19,999,645 (ends 19,999,795), right mate starts at 21,000,058; apparent insert 1,000,563 bp vs median 500 bp

All three signals agree: breakpoints at 20,000,000 / 21,000,000, deleted size 1,000,000 bp.

## Evidence vs. precision limits

- The **evidence** is the raw breakpoint estimate (method: split_reads): median junction position across split reads that bracket the depth interval; nucleotide-resolution evidence.
- The **reported precision** is coarser than the evidence: the task requires rounding each breakpoint to the nearest 100 kb, which caps coordinate precision at +/-50 kb. Here the raw breakpoints already fall on exact 100 kb grid positions, so rounding changes nothing, but in general the TSV coordinates must be read as +/-50 kb intervals around the raw values (preserved in `output/qc.json`, `breakpoints.raw_*`).
- Signal-specific resolution: read-depth bins localize the edges to ~10 kb (bin boundary); the spanning discordant pair brackets each junction to within the insert-size spread (few hundred bp); split reads give nucleotide-level junctions. The final call uses the finest available evidence (split reads).

## Specificity notes (segmental duplications on chr22)

chr22 contains large segmental duplications; copies of the deleted segment also reside near 18.4-18.9 Mb, 21.1-21.2 Mb, 21.4-21.5 Mb and 24.2 Mb. This causes two classes of artifacts that were explicitly checked and rejected:
- *Secondary depth dips.* Additional 100 kb windows with ratio < 0.3 (see `qc.json` artifact_audit) consist almost entirely of multi-copy k-mers (paralogous sequence); none shows junction split reads or spanning discordant pairs, consistent with reads from those loci being placed at paralogous positions rather than true deletion.
- *Spurious split clusters.* Split alignments whose junctions do not bracket the depth interval (e.g. apparent junctions inside the deleted segment or at ~1.5 kb / ~0.83 Mb offsets in duplicated sequence near 21.13 Mb / 21.17 Mb) were excluded from the call.

## Output files

- `output/deletion.tsv` - final call (chrom, start_100kb, end_100kb, size_bp, supporting_signals).
- `output/qc.json` - QC metrics, raw breakpoint evidence, precision notes and artifact audit.
- `output/analysis.py` - this analysis (self-contained, numpy only).
- `output/report.md` - this report.

