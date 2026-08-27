# Frozen splice reference

- `GRCh38_chr9.fa.gz`: Broad hg38/v0 chromosome 9.
- `ensembl112_protein_coding_exons.tsv.gz`: Ensembl Genes release 112 BioMart export of protein-coding transcript exon coordinates on 2026-08-17.

Report introns as 1-based inclusive skipped intervals. A junction is novel only when its `(chrom,intron_start,intron_end)` is absent from every supplied transcript. A cryptic exon interval between two ordered novel junctions is `left_intron_end+1` through `right_intron_start-1`.
