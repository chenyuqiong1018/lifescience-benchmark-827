# ENO1 tumor-versus-normal fold change

## Answer

ENO1 is **higher in tumor**. The tumor-versus-normal fold change is `4.80663981223069`, and the log2 fold change is `2.2650286985085`.

| Field | Value |
|---|---:|
| Gene | ENO1 |
| Tumor value | 350385456.451912 |
| Normal value | 72896133.2946858 |
| Fold change (`Tumor / Normal`) | 4.80663981223069 |
| log2 fold change | 2.2650286985085 |

Source: `Proteomic_data .xlsx`, sheet `Tumor vs Normal`.

## Validation and interpretation

- The designated workbook has 1 visible sheet, 3851 rows including the header, and 15 columns.
- ENO1 appears exactly once (worksheet row 474); all effect fields are present, numeric, non-formula cells, and visible.
- Direct calculations agree with the workbook's rounded `Ratio` and `FC` values (`4.81`) and rounded `log2FC` (`2.27`).
- A T1-only MarkItDown conversion of this same workbook preserved the ENO1 row and the same rounded values.
- This is a descriptive fold-change calculation. The supplied summary row is not used to invent a confidence interval or a new hypothesis test, and no claim of statistical significance is made here.
- No physical unit is assigned because the input does not specify one. External proteome annotations are not mixed into this file-derived calculation.

The unrelated RNA/m6A workbook was not opened or used.
