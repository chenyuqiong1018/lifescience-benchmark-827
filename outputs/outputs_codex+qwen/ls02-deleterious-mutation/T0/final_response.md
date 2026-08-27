Task complete. The analysis identified a single high-confidence mosaic nonsense SNV:

## Result

| Field | Value |
|---|---|
| Variant | **chr9:127661125 G>T** (GRCh38, 1-based) |
| Gene | **STXBP1** (ENSG00000136854.26) |
| Consequence | **stop_gained** — p.(Glu117Ter), codon GAA→TAA (MANE Select ENST00000373299.5; consistent across all 17 STXBP1 transcripts) |
| Allele fraction | **0.1935** (18 alt / 93 total Q≥20 reads; +14/−4 strands) |

## How it was found

Since no aligner binaries were installed, I built a complete pure-Python/numpy pipeline (`output/analysis.py`): 16-mer seed-and-verify mapping over protein-coding CDS±100 bp targets → ungapped placement verification with end-trimming → Q≥20 CDS pileup → 968 SNV candidates → transcript-aware codon reconstruction (26 stop-gained) → **allele co-occurrence filter**.

The co-occurrence step was decisive: the dataset contains aggressive decoys — most notably a FAM78A locus with ~50 co-occurring candidate alleles at 600–800× depth and near-zero "clean" support (paralog mismapping), plus germline-range (AF≈0.5–0.57) stop variants in C5/ALAD/GRIN3A. The STXBP1 variant is the only stop-gained candidate where **all 18 alt reads are clean** (carry no other candidate allele, mean BQ 33, most reads differ from GRCh38 at exactly the variant base), on both strands, at a clearly mosaic AF of ~0.19. STXBP1 is among the genome's most LoF-intolerant genes (gnomAD v2.1.1 pLI = 1.0, LOEUF ≈ 0.05; haploinsufficiency causes STXBP1 epileptic encephalopathy, OMIM 612164). A whole-chr9 paralog scan confirmed the locus is unique.

## Deliverables (all in `output/`)
- **variant.tsv** — `chr9 127661125 G T STXBP1 stop_gained 18 93 0.1935`
- **evidence.json** — reference/annotation versions (GRCh38 primary assembly chr9; GENCODE v47/Ensembl 113), pipeline parameters, read/coverage stats, all 26 stop-gained candidates with co-occurrence metrics, read-level evidence, and the LoF-constraint prior
- **analysis.py** — deterministic, re-runnable pipeline (~6 min, inputs-only)
- **report.md** — full methods, decoy-rejection table, evidence, and interpretation