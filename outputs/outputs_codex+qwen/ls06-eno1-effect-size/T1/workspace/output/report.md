# ENO1 Tumor-vs-Normal Effect Size (Proteomics)

## Result

| Field | Value |
|---|---|
| Gene | ENO1 (Alpha-enolase, UniProt P06733, gene_id 2023) |
| Tumor abundance | 350,385,456.451912 |
| Normal abundance | 72,896,133.2946858 |
| Fold change (Tumor / Normal) | **4.8066** |
| log2 fold change | **2.2650** |
| Source file | `Proteomic_data .xlsx` |
| Source sheet | `Tumor vs Normal` |

## Direction

**ENO1 is upregulated in tumor relative to normal tissue.** Tumor abundance is
~4.81x the normal abundance (fold change > 1, log2FC = +2.27 > 0), so the
effect direction is tumor-up / up in tumor.

## Method

1. Loaded `inputs/Proteomic_data .xlsx` (openpyxl, data-only values) and read
   sheet `Tumor vs Normal` (3,850 protein rows).
2. Located the row with `gene == "ENO1"` (row index 472, protein P06733).
3. Computed `fold_change = Tumor / Normal` and
   `log2_fold_change = log2(fold_change)` directly from the supplied `Normal`
   and `Tumor` abundance columns; no values were imputed or rescaled.

## Validation

- Workbook-reported columns for the same ENO1 row: `FC = 4.81`,
  `log2FC = 2.27`, `Ratio = 4.81` (stored rounded to 2 decimals in the sheet),
  `p.value = 0.031`, `adj.Pval = 0.226`. These agree with the values computed
  here from the raw abundances (4.8066 and 2.2650) within the sheet's rounding.
- The unrelated workbook `MeRIP_RNA_result.xlsx` (transcript/m6A decoy) was not
  used, per the input contract.

## Deliverables

- `output/eno1_effect.json` — machine-readable result with keys
  `gene, tumor_value, normal_value, fold_change, log2_fold_change, source_file, source_sheet`.
- `output/analysis.py` — reproducible script that regenerates the JSON.
- `output/report.md` — this report.
