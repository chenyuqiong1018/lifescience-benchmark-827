Completed the T0 large-deletion analysis with the selected whitelist skills `chromosome_analysis`, `genome_annotation`, and `code_execution_analysis`.

The supported call is a 1,000,000 bp deletion on chr22 with 100-kb-rounded boundaries at 20,000,000 and 21,000,000. Evidence comprises 10 consecutive callable 100-kb bins with zero uniquely mapped read depth, 2 correctly oriented FR pairs spanning the interval, and 4 exact junction-spanning reads. Precision and coordinate limits are documented in `workspace/output/report.md`.

Created:

- `workspace/output/deletion.tsv`
- `workspace/output/qc.json`
- `workspace/output/analysis.py`
- `workspace/output/report.md`
