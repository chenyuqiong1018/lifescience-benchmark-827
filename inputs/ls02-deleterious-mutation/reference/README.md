# Frozen reference bundle

- `GRCh38_chr9.fa.gz`: chromosome 9 from the Broad Institute GATK GRCh38 primary assembly bundle.
- `gencode.v47.chr9.annotation.gtf.gz`: chromosome-9 records extracted without modification from the GENCODE v47 primary-assembly GTF.

The task FASTQ contains chromosome-9 reads. Both files are participant-visible so that alignment and transcript-aware consequence reconstruction are reproducible without an implicit genome or annotation version. Coordinates in output artifacts are GRCh38, 1-based. File hashes are recorded in `docs/inputs/SHA256SUMS.tsv`.
