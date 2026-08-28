# ENO1 significance audit (Tumor vs Normal proteomics)

## Question
Is ENO1 differentially abundant at FDR 0.05 in the supplied proteomics results, using the multiple-testing **adjusted** p-value (not the raw p-value)?

## Data and method
- Source file: `Proteomic_data .xlsx` (inputs/)
- Source sheet: `Tumor vs Normal` (the only sheet; 3,850 protein rows)
- Columns used: `gene` for selection, `adj.Pval` as the adjusted p-value, `p.value` recorded only as the raw (unadjusted) value for reference. `MeRIP_RNA_result.xlsx` is an unrelated transcript/m6A file and was not used.
- Rule applied: significant iff `adj.Pval < 0.05`.

## Result
| gene | protein | raw p.value | adj.Pval | FDR threshold | significant |
|------|---------|-------------|----------|---------------|-------------|
| ENO1 | P06733  | 0.031       | 0.226    | 0.05          | No          |

ENO1 (alpha-enolase, P06733) shows higher tumor abundance (FC ~4.81, log2FC ~2.27), and its raw p-value (0.031) is nominally below 0.05. However, after multiple-testing correction its adjusted p-value is **0.226**, well above the 0.05 FDR cutoff, so **ENO1 is not significant at FDR 0.05**.

## Calibration notes
- The workbook's `is_sig` flag marks ENO1 `True`; it is derived from the **raw** p-value (p < 0.05), not the adjusted one. Across all 3,850 rows, none has `adj.Pval < 0.05` (minimum observed adj.Pval ~0.127), while 708 rows have raw p < 0.05. This confirms `is_sig` tracks raw p-values and must not be read as an FDR call.
- The raw p-value was never relabeled as adjusted; only the `adj.Pval` column feeds the significance decision in `output/eno1_significance.json`.

## Reproduction
`powershell
python output/analysis.py   # regenerates output/eno1_significance.json
`

## Machine-readable output
`output/eno1_significance.json` with keys `gene, adjusted_p_value, fdr_threshold, significant, source_file, source_sheet`.
