# Input contract: ls08-multiome-column-match

## Files and roles

- `multiome.match.atac.rna.q1.atac.tsv.gz`: ATAC-bin features by eight permuted population columns.
- `multiome.match.atac.rna.q1.rna.tsv.gz`: RNA-TPM features by eight population columns.

## Provenance and redistribution boundary

Retrieved on 2026-08-14 from `Genentech/compbiobench-data-v1`, repository revision `c673f0855fce09d320f1677f168f7864eec52c1a`, for task `multiome-match-atac-rna-q1`. Files retain the upstream dataset terms; this repository asserts no new license.

## Schema, units, and missing values

Both files are gzip-compressed tabular matrices. Population labels and feature identifiers are strings; matrix entries retain upstream units (ATAC-bin signal and RNA TPM). Missing numeric values remain missing and must not be imputed to zero unless a frozen preprocessing policy explicitly says so.

## Integrity and readiness boundary

See `../SHA256SUMS.tsv`. The inputs are complete; the remaining formal blocker is evaluator-side acceptance of the hidden permutation and normalization policy, not a missing participant-visible file.
