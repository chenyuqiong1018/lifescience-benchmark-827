# Annotation-locked cryptic-exon report

## Result

The unique high-support event is a **53-bp cryptic exon in GNG10** at
`chr9:111664537-111664589` (GRCh38; 1-based inclusive).

| junction | 1-based inclusive intron | final read support | discovery support | motif |
|---|---:|---:|---:|---|
| left | chr9:111661716-111664536 | 45 | 20 | GT-AG |
| right | chr9:111664590-111666814 | 42 | 21 | GT-AG |

There are 77 complete-read alignments overlapping the exon; 14
span both novel junctions. The interval is exactly `111664536 + 1` through
`111664590 - 1`.

## Annotation and sequence validation

The supplied Ensembl Genes release 112 table was converted into complete transcript
exon chains and a set of 1-based inclusive introns. Neither selected intron tuple
occurs in any supplied protein-coding transcript. Both boundaries were independently
extracted from the supplied GRCh38 chr9 FASTA and have plus-strand canonical `GT-AG`
motifs. The event maps to the supplied protein-coding GNG10 gene model: the left
junction starts immediately after an annotated exon and the right junction ends
immediately before the next annotated exon.

## Computational workflow

An exon-sequence index anchored split reads to protein-coding gene spans. Candidate
boundaries were resolved in priority order: exact frozen-annotation junction,
strand-consistent canonical junction, then other split. Novel donor- and acceptor-side
junctions were paired only inside the same gene when enclosing 15-500 bp. The GNG10
pair was the only qualifying event and each side had at least 20 discovery reads.

The flanking exons and candidate exon were assembled into a sequence-validated local
transcript. All 200,414 107-bp reads were tested in both
orientations; final support requires a full-read alignment with at most
3 mismatches and at least 8 aligned bases on each side.
Online sequence or annotation services were not queried because the supplied frozen
reference is the required novelty authority.

Run `python analysis.py` to regenerate and self-validate all outputs.
