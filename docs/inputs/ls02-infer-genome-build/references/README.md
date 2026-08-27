# Diagnostic chromosome references

The three gzip FASTA files are exact chromosome 20 sequences retrieved on 2026-08-17: UCSC hg18, Broad hg19/v0 (GRCh37), and Broad hg38/v0 (GRCh38). The evaluator checks every VCF REF allele at its declared 1-based coordinate. Source URLs and SHA-256 values are recorded in `reference_manifest.json`. T2T/hs1 is listed in the prompt as a candidate; if its reference file is absent, it must be excluded explicitly rather than treated as a mismatch. The final call still requires one supplied reference to dominate reproducibly.
