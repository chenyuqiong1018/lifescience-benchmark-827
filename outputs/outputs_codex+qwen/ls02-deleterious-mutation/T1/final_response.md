Analysis complete. Here's the summary:

## Result: mosaic nonsense SNV in **STXBP1**

| Field | Value |
|---|---|
| Variant | **chr9:127,661,125 G>T** (GRCh38, 1-based) |
| Gene | **STXBP1** (ENSG00000136854), chr9q34.11 |
| Consequence | **stop_gained** — canonical transcript ENST00000373299.5 (STXBP1-201), exon 6/19, codon `GAA>TAA`, **p.Glu117Ter** (affects 17 STXBP1 transcripts) |
| Allele fraction | **0.1277** (6 ALT / 47 unique-molecule reads) |
| LoF intolerance | gnomAD v4: **pLI = 1.0**, o/e LoF = 0.038 (CI 0.017–0.099), lof_z = 7.23 |

## How it was established
- **Pipeline:** bwa mem 0.7.17 → samtools markdup → duplicate-filtered BAM → bcftools mpileup/call + a streaming mosaic scan of all 16.3M covered sites → custom strand/phase-aware GENCODE v47 annotation engine (2,998 protein-coding transcripts reconstructed from the reference).
- **Mosaicism:** AF 0.128 rejects germline het (binomial P = 8.9e-08 vs AF=0.5) and sequencing error (P = 7.6e-06 vs 1% error); the diploid caller misses the site entirely.
- **High confidence:** 6 independent ALT molecules (distinct starts), both strands, all baseQ 33, all MAPQ 60, no nearby variants, clean non-homopolymer context, counts identical with/without duplicates.
- **Uniqueness:** 5 stop-gained sites were found (GRIN3A, ALAD, COL27A1, C5, STXBP1) — the other four are either germline hets (AF ~0.5), strand-bias artifacts (COL27A1), or in LoF-tolerant genes; STXBP1 is the only highly LoF-intolerant hit (the other genes' gnomAD pLI ≤ 0.05).

## Deliverables (in `output/`)
- **variant.tsv** — `chr9 127661125 G T STXBP1 stop_gained 6 47 0.1277`
- **evidence.json** — read-level support, binomial stats, gnomAD constraint, rejected candidates, full provenance
- **analysis.py** — reproducible end-to-end pipeline (re-runs verified; reuses intermediates in `work/`)
- **report.md** — full methods, evidence tables, and version documentation

**Versions documented:** GRCh38 primary assembly chr9 (Broad GATK bundle, 138,394,717 bp), GENCODE v47 (Ensembl 113, 2024-07-19), gnomAD v4 constraint, samtools/htslib 1.19.2/1.19, bcftools 1.19, bwa 0.7.17-r1188; coordinates GRCh38 1-based.