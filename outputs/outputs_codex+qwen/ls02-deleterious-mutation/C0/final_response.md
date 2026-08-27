Task complete. The high-confidence mosaic nonsense SNV is:

## Called variant
| chrom | pos | ref | alt | gene | consequence | alt_reads | total_reads | allele_fraction |
|---|---|---|---|---|---|---|---|---|
| chr9 | 127661125 | G | T | **STXBP1** | stop_gained | 18 | 93 | **0.1935** |

**`output/variant.tsv`** contains exactly that row (columns chrom, pos, ref, alt, gene, consequence, alt_reads, total_reads, allele_fraction; AF on 0–1 scale).

## How it was found
1. **Alignment** — 431,396 single-end reads aligned to GRCh38 chr9 with minimap2 2.31 (`-ax sr`); 99.1% primary mapping rate.
2. **Pileup/calling** (`output/analysis.py`, pure-Python, reproducible) — strand-aware counts with MAPQ≥20/baseQ≥20, then mosaic-SNV filters (≥4 alt reads, AF 0.02–0.50, both strands, mean alt baseQ≥25, no read-end artifact) → 328 candidates.
3. **Consequence reconstruction** against **GENCODE v47** protein-coding transcripts (codon-level; GTF reading frames cross-checked against the reconstructed CDS — coordinate systems fully consistent). Two stop-gained sites emerged; selection required a mosaic AF plus a highly LoF-intolerant gene (gnomAD LOEUF < 0.35, pLI ≥ 0.9).

## Why this variant
- **Mosaicism:** AF 18/93 = 0.1935 (far below the 0.5 germline-het expectation), alt on both strands, all alt reads MAPQ 60/baseQ 33, distributed across read positions, independent recount matches, and dbSNP shows no record at the site.
- **LoF intolerance:** STXBP1 (gnomAD v4, GRCh38): obs 3 vs exp 78.5 LoF variants, **LOEUF = 0.099, pLI = 1.0** — among the most constrained human genes; heterozygous LoF causes DEE4 (OMIM 612164) via haploinsufficiency.
- **Consequence:** c.349G>T → GAA>TAA → **p.Glu117\*** in MANE Select ENST00000373299.5 (p.Glu103\* in shorter isoforms), stop-gained in all 17 STXBP1 transcripts.
- **Rejected alternative:** chr9:122477222 G>T in OR1J1 (p.Cys235\*) — AF 0.462 (germline-het-like) in an olfactory receptor gene with no LoF constraint (`no_exp_lof`). A relaxed-threshold scan confirmed no other nonsense site has meaningful support.

## Deliverables (all in `output/`)
- `variant.tsv` — called variant
- `evidence.json` — full evidence bundle (versions: GRCh38/GATK bundle, GENCODE v47, minimap2 2.31, gnomAD v4 constraint; filters, recount, selection decisions)
- `analysis.py` — end-to-end reproducible pipeline (regenerates all three output artifacts)
- `report.md` — methods, results, mosaicism/LoF-intolerance evidence, rejected candidate, limitations

Intermediates (SAM, candidate lists, cached gnomAD response) are kept under `work/` for traceability.