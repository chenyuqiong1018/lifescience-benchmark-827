# ENO1 significance audit (FDR 0.05)

## Task

Retrieve ENO1's adjusted p-value from the supplied proteomics results and give a
threshold-calibrated interpretation at FDR 0.05.

## Source

- Source file: `inputs/Proteomic_data .xlsx`
- Source sheet: `Tumor vs Normal` (only sheet in the workbook)
- ENO1 row: protein `P06733`, gene_id 2023, log2FC = 2.27, FC = 4.81
- Comparison: Tumor vs Normal

## Extracted values

| Column | Value | Meaning |
|---|---|---|
| `p.value` | 0.031 | raw (unadjusted) p-value - context only |
| `adj.Pval` | 0.226 | adjusted p-value (multiple-testing corrected) - used for the decision |

The sheet carries both a raw p-value column (`p.value`) and an adjusted
p-value column (`adj.Pval`). Per the task rules, the **adjusted** value
(`adj.Pval` = 0.226) is used for the significance call. The raw p-value is
**not** relabeled or treated as an adjusted value.

## Decision at FDR 0.05

- Adjusted p-value: **0.226**
- FDR threshold: **0.05**
- 0.226 > 0.05 -> **NOT significant** after multiple-testing correction.

## Interpretation

Although ENO1's raw p-value (0.031) is below 0.05 and the sheet's
`is_sig` flag is `True`, the adjusted p-value (0.226) exceeds the FDR 0.05
threshold. Across the 3850 quantified proteins, only 0 pass
adj.Pval <= 0.05, so ENO1's raw-level evidence does not survive
multiple-testing correction. ENO1 (log2FC = 2.27, strongly up-regulated in
tumor) is therefore **not statistically significant at FDR 0.05** in this
comparison; its up-regulation is a descriptive fold-change observation rather
than an FDR-controlled hit. Note the sheet's `is_sig` flag appears to reflect
an unadjusted criterion and should not override the adjusted p-value.

## Reproducibility

Run `python output/analysis.py` from the workspace root. It regenerates
`output/eno1_significance.json` and this report from the frozen input files.
