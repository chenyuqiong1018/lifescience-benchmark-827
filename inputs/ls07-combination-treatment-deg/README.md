# Input contract: ls07-combination-treatment-deg

## Files and roles

- `counts_raw_unfiltered.csv`: raw integer gene-by-sample count matrix.
- `sample_layout.csv`: sample identifiers and condition labels; only the frozen six-sample contrast belongs in the model.
- `ensg_to_gene_name.tsv`: capsule-supplied Ensembl-to-display-name support mapping. DE statistics remain keyed by Ensembl ID.

## Provenance and redistribution boundary

Files came from the official BixBench v1.5 `bix-43` capsule retrieved on 2026-08-14. The dataset card is Apache-2.0; the capsule does not state independent redistribution terms for the original experimental data. Preserve upstream attribution and keep the files inside the benchmark workspace unless those terms have been checked.

## Schema, units, and missing values

Counts are non-negative integers. `sample_layout.csv` uses the supplied group strings verbatim. Missing/unusable mapping rows do not remove genes from the DE table. Independent-filtering missing statistics are represented as empty CSV cells or JSON `null`, never string `NaN` or invented zeroes.

## Integrity

See `../SHA256SUMS.tsv`; the per-run manifest is generated before execution.
