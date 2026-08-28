# ENO1 Tumor vs Normal Effect Size

## Result

ENO1 (alpha-enolase, UniProt P06733) is **upregulated in tumor relative to
normal tissue** in the supplied proteomics dataset.

| Quantity | Value |
|---|---|
| Tumor abundance | 350,385,456.45 |
| Normal abundance | 72,896,133.29 |
| Fold change (Tumor / Normal) | 4.8066 |
| log2 fold change | +2.2650 |
| Workbook-reported FC (rounded) | 4.81 |
| Workbook-reported log2FC (rounded) | 2.27 |
| p.value / adj.Pval | 0.031 / 0.226 |

## Direction

Fold-change direction: **tumor > normal (up in tumor)**. The tumor abundance
is ~4.8-fold higher than the normal abundance (log2FC = +2.27, positive).

## Method

- Source workbook: `inputs/Proteomic_data .xlsx`, sheet `Tumor vs Normal`
  (3,850 protein rows; the ENO1 row was located by exact match on the `gene`
  column). The unrelated `MeRIP_RNA_result.xlsx` transcriptomics workbook was
  not used, per the input contract.
- `fold_change` was computed as `Tumor / Normal` = 350385456.451912 /
  72896133.2946858 = 4.8066, and `log2_fold_change` as
  `log2(fold_change)` = +2.2650, from the raw (unrounded) abundance values.
- These computed values agree with the workbook's own rounded `FC` (4.81) and
  `log2FC` (2.27) columns, confirming internal consistency.
- All computations are in `output/analysis.py`; machine-readable results are
  in `output/eno1_effect.json`.

## Notes

- Numeric values were used as supplied; no missing-value imputation was
  needed (the ENO1 row has no empty cells in the relevant columns).
- Abundance values are normalized protein abundances as reported in the
  workbook; no physical unit is asserted beyond the source data.
