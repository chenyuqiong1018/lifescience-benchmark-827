# Mosaic nonsense SNV report

The high-confidence mosaic nonsense SNV is **GRCh38 chr9:127661125 G>T in STXBP1**. GENCODE v47 CDS reconstruction changes `GAA` to `TAA` across 21 protein-coding STXBP1 transcripts, yielding a `stop_gained` consequence.

At base quality Q20 or higher, 17 reads support T and 73 support the GRCh38 G allele (90 total; allele fraction 0.188889). ALT evidence occurs on both read orientations (13 forward, 4 reverse), spans 13 distinct alignment starts, and has Q33 at the variant base. ALT-supporting reads have median one mismatch to the local reference alignment (maximum three).

An exome-wide protein-coding CDS scan also produced nonsense-like candidates in OR1J1, OR1B1, NPDC1, and PAEP. Their allele fractions included 0.46, 0.50, and 1.0, or they lacked the task's highly LoF-intolerant-gene qualifier. STXBP1 uniquely combines a mosaic allele fraction, strong bidirectional evidence, and the requested loss-of-function-intolerant protein-coding gene context. No external constraint score was fetched.

Reference: Broad GATK GRCh38 primary assembly chromosome 9; coordinates are 1-based. Annotation: GENCODE v47 primary-assembly chromosome-9 GTF. Exact input SHA-256 values and transcript-level codon reconstructions are in `evidence.json`.
