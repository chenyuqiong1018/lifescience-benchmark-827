# Input contract: ls10-treatment-response-model

## File and role

- `data.xlsx`: sheet `Sheet1` contains `Efficacy`, `Age`, `Gender`, `BMI`, and distractor covariates. Only the frozen outcome and three named predictors belong in the model.

## Provenance and redistribution boundary

Retained from the official BixBench v1.5 capsule for `bix-51` on 2026-08-14. The BixBench dataset card is Apache-2.0; the underlying workbook retains its upstream terms.

## Schema, units, and missing values

Workbook labels, categorical spellings, and empty cells are frozen. Complete-case handling applies only to the specified model variables. Missing values are not converted to zero or a new category.

## Integrity

See `../SHA256SUMS.tsv`; the per-run manifest freezes the copied workbook and this input contract.
