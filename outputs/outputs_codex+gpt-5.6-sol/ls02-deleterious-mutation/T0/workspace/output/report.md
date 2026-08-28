# Mosaic nonsense SNV report

The T0 whitelist skills informed local annotation, functional-consequence reconstruction, mosaic proportion reporting, and reproducible execution. No remote skill tools or external databases were called.

The selected variant is **GRCh38 chr9:127661125 G>T in STXBP1**, producing `GAA>TAA` and `stop_gained` across 21 GENCODE v47 protein-coding transcripts. Read evidence is 17 ALT and 73 REF reads (90 total), AF 0.188889. A Wilson 95% interval for the observed ALT proportion is approximately 0.121–0.282; this quantifies sampling precision rather than changing the called AF.

ALT reads occur on both orientations (13 forward, 4 reverse), span 13 distinct alignment starts, and have Q33 at the variant position. The combination of a subclonal AF, bidirectional evidence, and the task-specified highly LoF-intolerant gene context distinguishes STXBP1 from the other locally reconstructed nonsense-like candidates. No external constraint score was fetched.

Reference: Broad GATK GRCh38 primary assembly chromosome 9, 1-based coordinates. Annotation: GENCODE v47 primary-assembly chr9 GTF. SHA-256 values and transcript-level evidence are recorded in `evidence.json`.
