# Input contract: ls06-eno1-significance-audit

## Files and roles

- `Proteomic_data .xlsx`: analysis target. Sheet `Tumor vs Normal` contains both raw and adjusted p-value columns.
- `MeRIP_RNA_result.xlsx`: unrelated transcript/m6A decoy.

## Provenance and redistribution boundary

Both workbooks were retained from the official BixBench v1.5 capsule for `bix-37` on 2026-08-14. The BixBench dataset card is Apache-2.0, but the underlying experimental-data redistribution terms are not asserted by this repository. Preserve upstream attribution and confirm those terms before redistribution beyond the evaluation workspace.

## Schema, units, and missing values

Workbook sheet names, column labels, cell types, and empty cells are frozen. P-values are unitless and must remain in `[0,1]`; empty cells remain missing rather than zero.

## Integrity

See `../SHA256SUMS.tsv`. Each run creates a task-local `INPUT_MANIFEST.sha256.tsv` from the copied read-only files.
