# ENO1 adjusted-p-value significance audit

## Result

ENO1's adjusted p-value is `0.226`. At FDR `0.05`, ENO1 is **not significant** because the decision expression `0.226 <= 0.05` evaluates to `false`.

| Field | Value |
|---|---|
| Gene | `ENO1` |
| Adjusted p-value (`adj.Pval`) | `0.226` |
| FDR threshold | `0.05` |
| Significant | `false` |
| Source file | `Proteomic_data .xlsx` |
| Source sheet | `Tumor vs Normal` |

## Provenance and audit checks

- ENO1 appears exactly once at worksheet row 474 in the single visible target sheet.
- The raw `p.value` (`0.031`) and adjusted `adj.Pval` (`0.226`) are distinct columns and values. The raw value is not relabeled as adjusted.
- The worksheet's existing `is_sig` value is `True`, but it is not used for this answer because it conflicts with the requested direct rule `adj.Pval <= 0.05`.
- A T1-only MarkItDown conversion of the same target workbook preserved the ENO1 row and both p-value columns.
- The controlled code-execution handoff used the same threshold expression; the saved local script performs the authoritative executable calculation.

This result means ENO1 does not meet the specified FDR criterion; it is not proof that ENO1 has no biological effect. No additional hypothesis test was invented. External proteome and biomarker databases were not mixed into this workbook-derived significance decision. The unrelated RNA/m6A workbook was not opened or used.
