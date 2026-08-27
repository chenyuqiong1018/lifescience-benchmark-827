# Input contract: ls06-eno1-effect-size

## Files and roles

- `Proteomic_data .xlsx`: analysis target. Use sheet `Tumor vs Normal`; relevant columns include `gene`, `Normal`, `Tumor`, `Ratio`, `FC`, `log2FC`, `p.value`, and `adj.Pval`.
- `MeRIP_RNA_result.xlsx`: unrelated transcript/m6A decoy. It must not be substituted for the proteomics workbook.

## Provenance and redistribution boundary

Both workbooks were retained from the official BixBench v1.5 capsule for `bix-37` on 2026-08-14. The BixBench dataset card is Apache-2.0, but that metadata does not establish new rights over the underlying experimental workbooks. Preserve upstream attribution and verify the original study's redistribution terms before republishing the workbooks outside the benchmark workspace.

## Schema, units, and missing values

Workbook cell types and sheet names are part of the frozen input. Numeric values are used as supplied; the task does not invent a new physical unit for normalized abundance. Empty workbook cells remain missing and must not be converted to zero.

## Integrity

Byte sizes and SHA-256 values are recorded in `../SHA256SUMS.tsv`. A run copies this directory read-only and creates `INPUT_MANIFEST.sha256.tsv` before the prompt is shown.
