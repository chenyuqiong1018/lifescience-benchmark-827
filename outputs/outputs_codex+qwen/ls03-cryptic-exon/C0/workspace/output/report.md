# Cryptic-exon analysis report

## Question
Identify the protein-coding HGNC gene containing the highly expressed cryptic
exon supported by two novel splice junctions. Novelty is assessed against the
supplied Ensembl release-112 protein-coding annotation.

## Data
- `inputs/cryptic.exon.q1.fq.gz`: 200,414 single-end reads, 107 bp each.
- `inputs/reference/GRCh38_chr9.fa.gz`: chromosome chr9 (138,394,717 bp, Broad hg38/v0).
- `inputs/reference/ensembl112_protein_coding_exons.tsv.gz`: Ensembl 112
  protein-coding exon annotation (283,075 transcript introns used as the
  known-junction set).

## Method (implemented in `analysis.py`, pure Python + NumPy)
No splice-aware aligner was available, so the pipeline:

1. builds a genome index restricted to 16-mers present in the reads (k-mers
   with >20 genomic copies are treated as repetitive);
2. places every read (both orientations) by seed/diagonal voting, verifying
   alignments base-by-base (<=3 mismatches / 107 bp; introns 30-200,000 bp);
   reads with several equally good loci (segmental duplications) are
   discarded, so every counted read is uniquely placed;
3. canonicalises each unique junction cluster to the GT-AG splice motif
   (+/-3 bp window, strand-aware);
4. calls a junction novel iff `(chrom, intron_start, intron_end)` (1-based
   inclusive skipped interval) is absent from every supplied transcript;
5. searches pairs of novel junctions whose interval
   `left_intron_end+1 .. right_intron_start-1` overlaps no annotated exon and
   lies inside a protein-coding gene;
6. resolves microhomology-ambiguous boundaries by re-aligning all reads to
   every canonical splice-site variant of the candidate and choosing the
   variant whose spliced product is supported by reads spanning both
   junctions;
7. quantifies the host gene against full-length canonical and cryptic isoform
   references.

## Result

**The cryptic exon lies in the protein-coding HGNC gene GNG10
(ENSG00000242616), chromosome chr9.**

| item | value |
|---|---|
| Gene | GNG10 (protein-coding, + strand, transcript ENST00000374293) |
| Cryptic exon | chr9:111664537-111664589 (53 bp) |
| Location | GNG10 intron 1 (111661716-111666814), between exon 1 (111661605-111661715) and exon 2 (111666815-111666946) |
| Left novel junction (exon 1 -> cryptic exon) | chr9:111661716-111664536, 49 reads |
| Right novel junction (cryptic exon -> exon 2) | chr9:111664590-111666814, 49 reads |
| Reads spanning both novel junctions | 21 |
| Cryptic-exon body coverage | 48-55x (mean 51.5x) |
| GNG10 uniquely-mapped reads | 506 |
| Reads spanning constitutive intron 1 (canonical isoform) | 0 |

### Novelty against the supplied Ensembl 112 annotation
- Left junction `chr9:111661716-111664536`: `(chrom, 111661716, 111664536)` is absent
  from every supplied transcript (annotated intron 1 is 111661716-111666814; the
  cryptic 3' splice site 111664536 is new). **Novel.**
- Right junction `chr9:111664590-111666814`: `(chrom, 111664590, 111666814)` is absent
  from every supplied transcript (the cryptic 5' splice site 111664590 is new; the
  downstream acceptor 111666814 is the annotated one). **Novel.**
- Both junctions carry canonical GT-AG dinucleotides on the + strand.
- The interval 111664537-111664589 does not overlap any annotated exon of any supplied
  transcript.

### Expression evidence
GNG10 is among the most highly expressed genes in this library, and
every uniquely-mapped GNG10 read is consistent with the cryptic isoform
(exon 1 - cryptic exon - exon 2 - exon 3): **0 reads span constitutive intron
1**, i.e. the cryptic exon is included in ~100% of GNG10 transcripts.
Direct evidence for the exon body: 21 reads span both novel junctions and
thereby cover the entire 53 bp exon; together with the remaining
junction reads the exon body is covered at 48-55x
(mean 51.5x). Junction support: 49 left + 49 right
(double-spanning reads counted in both).

Boundary note: exon 1 ends in `...AAG` and the genomic cryptic-exon region
begins `AAGTTG...` (3-bp AAG microhomology), which makes raw split alignment
of the left junction ambiguous by 3 bp. Re-aligning reads to candidate
spliced transcripts shows the transcript keeps the AAG on the exon side, so
the cryptic exon starts at 111664537 (acceptor AG at 111664535-111664536); this is the
only variant supported by reads spanning both junctions with 0-1 mismatches.

### Ruled-out alternatives
`junctions.tsv` lists 721 detected junctions, 272 of them novel
against Ensembl 112. The other novel junctions are
+/-1-3 bp shifted copies of annotated introns, frequently supported on both
strands - the signature of reads from inverted segmental-duplication copies
of highly expressed genes (RPS6, RPL7A, HSPA5, SET, ...). Every interval
defined by those junction pairs overlaps an annotated exon (the shifted exon
of the duplicated copy), so none qualifies as a cryptic exon. The GNG10
pair is the only well-supported pair of novel junctions whose intervening
interval is unannotated, gene-internal and expressed.

*Generated by `output/analysis.py` (deterministic; re-run with
`python output/analysis.py`).*
