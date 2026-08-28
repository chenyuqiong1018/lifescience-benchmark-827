# ENO1 tumor-versus-normal effect size

## Result

ENO1 is **higher in tumor**. Using the supplied abundance values, the tumor-versus-normal fold change is `4.80663981223069` and the log2 fold change is `2.2650286985085`.

- Tumor value: `350385456.451912`
- Normal value: `72896133.2946858`
- Calculation: `Tumor / Normal`
- Source file: `Proteomic_data .xlsx`
- Source sheet: `Tumor vs Normal`

## Focused data checks

- The target workbook contains 1 visible analysis sheet with 3851 rows including the header and 15 columns.
- ENO1 occurs exactly once, at worksheet row 474.
- The ENO1 `gene`, `Normal`, `Tumor`, `Ratio`, `FC`, and `log2FC` cells are all present.
- The direct calculation agrees with the workbook's rounded `Ratio`/`FC` (`4.81`) and `log2FC` (`2.27`) values.

The workbook does not specify a physical unit for the abundance values, so none is invented. No inferential significance claim is made from this single requested effect-size calculation. The unrelated RNA/m6A workbook was not opened or used.
