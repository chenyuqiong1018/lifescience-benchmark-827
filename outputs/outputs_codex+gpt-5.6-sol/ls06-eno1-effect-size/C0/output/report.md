# ENO1 tumor-versus-normal effect size

ENO1 is higher in tumor than normal in the supplied proteomics results.

| Gene | Normal value | Tumor value | Tumor / normal fold change | log2 fold change |
|---|---:|---:|---:|---:|
| ENO1 | 72896133.2946858 | 350385456.451912 | 4.806639812231 | 2.265028698508 |

The calculation uses `Tumor / Normal` from `Proteomic_data .xlsx`, sheet `Tumor vs Normal`. The positive log2 fold change (2.265) corresponds to an approximately 4.81-fold increase in tumor. The workbook's `Ratio`, `FC`, and `log2FC` cells (4.81, 4.81, and 2.27) agree with the direct calculation after rounding to two decimals.

No physical unit is invented for the normalized abundance values. `MeRIP_RNA_result.xlsx` is an unrelated RNA/m6A workbook and was not used.
