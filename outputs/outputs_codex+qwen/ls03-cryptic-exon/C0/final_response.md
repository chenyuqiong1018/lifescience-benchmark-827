All deliverables are complete and verified. 

## Summary

**Answer: the cryptic exon is in the protein-coding HGNC gene `GNG10` (ENSG00000242616), chr9:111,664,537–111,664,589 (53 bp), inside GNG10 intron 1.**

| Item | Value |
|---|---|
| Cryptic exon | `chr9:111664537-111664589` (53 bp, between GNG10 exon 1 and exon 2) |
| Left novel junction (exon 1 → cryptic exon) | `chr9:111661716-111664536` — **49 reads**, canonical GT–AG, absent from all Ensembl 112 transcripts |
| Right novel junction (cryptic exon → exon 2) | `chr9:111664590-111666814` — **49 reads**, canonical GT–AG, absent from all Ensembl 112 transcripts |
| Expression | 506 uniquely-mapped GNG10 reads; 21 reads span **both** junctions; exon body covered 48–55× (mean 51.5×); **0** reads span the constitutive intron 1 → the exon is included in ~100% of GNG10 transcripts |

**How it was done** (no aligner was available, so `output/analysis.py` implements a full pure-Python/NumPy pipeline):
1. Read-k-mer-restricted genome index (k=16, repeats capped), seed/diagonal-voting placement with base-level verification (≤3 mm), uniqueness filtering against segmental duplications.
2. Junction calling → GT–AG canonicalization → novelty tested as `(chrom, intron_start, intron_end)` tuples against **all introns of the supplied Ensembl 112 annotation** (1-based inclusive, per the reference README).
3. Cryptic-exon search: the GNG10 pair was the **only** well-supported novel-junction pair whose `left_intron_end+1 .. right_intron_start-1` interval overlaps no annotated exon and lies inside a protein-coding gene.
4. A 3-bp `AAG` microhomology at the exon-1/cryptic-exon boundary was resolved by re-aligning all reads to every motif-valid splice variant — only the reported coordinates are supported by reads spanning both junctions (0–1 mm).
5. Traps ruled out: ~270 other novel junctions are ±1–3 bp shifted, often opposite-strand copies of annotated introns from inverted segmental-duplication copies of highly expressed genes (RPS6, RPL7A, HSPA5, SET…); all their candidate intervals overlap annotated exons.

**Deliverables** (in `output/`): `cryptic_exon.tsv`, `junctions.tsv` (721 junctions, 272 novel), `analysis.py` (deterministic; re-runs end-to-end from `./inputs` in ~3.5 min), and `report.md`.