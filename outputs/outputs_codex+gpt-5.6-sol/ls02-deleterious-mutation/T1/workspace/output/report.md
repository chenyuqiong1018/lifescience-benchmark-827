# Mosaic nonsense SNV report

The controlled T1 skill set was freshly installed and opened in the specified order. It guided local genome annotation, SNV consequence reconstruction, mosaic-proportion reporting, and reproducible execution; all remote tools were left unused.

The identified variant is **GRCh38 chr9:127661125 G>T in STXBP1**, with `GAA>TAA` and `stop_gained` reconstructed across 21 GENCODE v47 protein-coding transcripts. Evidence is 17 ALT and 73 REF reads (90 total), AF 0.188889; the Wilson 95% interval for the observed proportion is approximately 0.121–0.282.

ALT evidence is bidirectional (13 forward, 4 reverse), spans 13 distinct alignment starts, and has Q33 at the variant base. No external VEP, constraint score, phenotype, or literature evidence was fetched.

Reference: Broad GATK GRCh38 primary assembly chr9, 1-based coordinates. Annotation: GENCODE v47 primary-assembly chr9 GTF. Exact hashes and transcript reconstructions are in `evidence.json`.
