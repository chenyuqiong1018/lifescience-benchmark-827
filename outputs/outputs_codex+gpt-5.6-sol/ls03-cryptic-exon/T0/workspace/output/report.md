# Transcript-centric cryptic-exon analysis

## Call

The unique highly supported event is a **53-bp cryptic exon in GNG10** at
`chr9:111664537-111664589` (GRCh38, 1-based inclusive).

- Left novel intron: `chr9:111661716-111664536` with 45 reads.
- Right novel intron: `chr9:111664590-111666814` with 42 reads.
- 77 complete-read alignments overlap the exon; 14 span both junctions.

Both introns are absent from every transcript in the supplied Ensembl Genes release
112 protein-coding exon table and both have plus-strand `GT-AG` boundaries. The exon
coordinates follow the supplied construction rule: `111664536 + 1 =
111664537` and `111664590 - 1 = 111664589`.

## Workflow

Following the transcriptome-analysis workflow, the analysis first reconstructed all
frozen transcript exon chains, built their known intron catalogue, and indexed chr9
exon sequence. Split reads were resolved within a containing gene; exact Ensembl 112
boundaries were preferred before canonical unannotated boundaries. Novel donor-side
and acceptor-side junctions were paired only inside one protein-coding gene. The top
pair had discovery supports 20 and 21
and was the only qualifying paired candidate.

The GNG10 flanking exons plus the candidate exon were then assembled into a local
transcript. All 200,414 107-bp reads were aligned in
both orientations with at most 3 mismatches across the complete read;
junction support required at least 8 bases on each side. Online Ensembl
or UCSC calls were deliberately not used because the task requires novelty relative
to the supplied frozen release rather than a mutable current database.

Run `python analysis.py` to regenerate all outputs deterministically.
