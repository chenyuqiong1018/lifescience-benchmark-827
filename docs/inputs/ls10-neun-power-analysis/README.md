# Input contract: ls10-neun-power-analysis

## File and role

- `NeuN_quantification.csv`: two labeled groups with observed NeuN measurements used for effect-size and power calculations.

## Provenance and redistribution boundary

Retained from the official BixBench v1.5 capsule for `bix-19` on 2026-08-14. The BixBench dataset card is Apache-2.0; the repository does not infer new rights over the underlying measurements.

## Schema, units, and missing values

Group labels and row order are frozen. Measurements retain their upstream quantitative unit. Calculations use sample standard deviations; missing observations remain missing and are not replaced by zero.

## Integrity

See `../SHA256SUMS.tsv`; each run records a task-local input manifest before execution.
