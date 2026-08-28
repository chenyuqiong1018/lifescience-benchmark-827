# ENO1 tumor-versus-normal effect size (proteomics)

## Source
- File: `inputs/Proteomic_data .xlsx` (sheet `Tumor vs Normal`)
- The unrelated workbook `MeRIP_RNA_result.xlsx` (transcript/m6A decoy) was **not** used.

## ENO1 values
| Quantity | Value |
|---|---|
| Normal abundance | 72,896,133.29 |
| Tumor abundance | 350,385,456.45 |
| Fold change (Tumor / Normal) | 4.8066 |
| log2 fold change | 2.2650 |

## Direction
ENO1 is **upregulated in tumor**: tumor abundance is ~4.81-fold the normal
abundance (log2FC = +2.27).

## Validation
- Exactly one ENO1 row in the sheet (of 3850 data rows).
- Recomputed FC 4.8066 matches the sheet's recorded FC (4.81) after rounding.
- Recomputed log2FC 2.2650 matches the sheet's recorded log2FC (2.27) after rounding.
- Sheet-level statistics for context: p.value = 0.031, adj.Pval = 0.226, is_sig = True.
