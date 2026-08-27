# Cryptic-exon analysis

## Result

A highly expressed **53-bp cryptic exon** was detected in the protein-coding HGNC
gene **GNG10** at **chr9:111,664,537-111,664,589**
(GRCh38, 1-based inclusive).

## Junction evidence

- Left novel intron: `chr9:111661716-111664536`;
  45 supporting reads.
- Right novel intron: `chr9:111664590-111666814`;
  42 supporting reads.
- 77 full-read alignments overlap the cryptic exon, and
  14 reads span both novel junctions.
- Both junctions use canonical `GT-AG` splice boundaries on the plus strand.

The interval follows the supplied rule exactly: left intron end + 1 =
111664536 + 1 = 111664537, and right intron start - 1 =
111664590 - 1 = 111664589.

## Novelty assessment

All adjacent exon pairs from the supplied `ensembl112_protein_coding_exons.tsv.gz`
were converted to 1-based inclusive introns. Neither complete intron tuple is present
in any supplied Ensembl Genes release 112 protein-coding transcript. Nearby annotated
junctions were resolved before novel calls, preventing sequencing-error boundary
shifts from being mislabeled as novel.

## Method

The discovery pass built a unique 21-mer index over supplied chr9 protein-coding
exons, completed candidate split alignments within the same gene span, and preferred
exact annotated boundaries or strand-consistent canonical boundaries. Novel donor-
and acceptor-side junctions were paired only when they enclosed a 15-500 bp interval
in one gene. The winning pair had discovery supports 20 and
21 and was unique.

For final counts, the flanking GNG10 exons and candidate cryptic exon were assembled
into a candidate transcript. All 200,414 FASTQ reads were aligned in
both orientations; alignments allowed at most 3 mismatches over
the complete 107-bp read, and junction counts required at least
8 aligned bases on each side. This produced 90
full-read candidate-transcript alignments.

## Reproduction

Run `python analysis.py` from any directory. It locates the repository from its own
path and deterministically rewrites the two TSV files and this report.
