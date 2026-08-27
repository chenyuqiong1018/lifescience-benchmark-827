# Input contract: ls08-enhancer-promoter-integration

## Files and roles

- `ep.interactions.q1.hic.csv`: physical-contact evidence for candidate enhancer-promoter pairs.
- `ep.interactions.q1.expr.csv`: CRISPR-expression evidence for the same candidate keys.

## Provenance and redistribution boundary

Retrieved on 2026-08-14 from `Genentech/compbiobench-data-v1`, revision `c673f0855fce09d320f1677f168f7864eec52c1a`, for task `ep-interactions-q1`. Files retain upstream terms; no new license is asserted.

## Schema, units, and missing values

Join the two modalities only by the supplied pair key. Keep physical contact and perturbation evidence as separate measurements. Numeric scales and missing values are preserved as supplied; empty values remain empty/null.

## Integrity and readiness boundary

See `../SHA256SUMS.tsv`. The source question mentions eight candidates while the supplied files contain `EP1` through `EP7`; this discrepancy is deliberately disclosed and must not be silently repaired. The frozen integration rule, gold artifacts and accepted static checker define scoring for the supplied EP1–EP7 files.
