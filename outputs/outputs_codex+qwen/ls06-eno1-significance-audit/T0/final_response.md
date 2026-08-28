All deliverables are complete and validated.

## Result: ENO1 is **not significant** at FDR 0.05

| Field | Value |
|---|---|
| Adjusted p-value (`adj.Pval`) | **0.226** |
| FDR threshold | 0.05 |
| Significant | **false** (0.226 > 0.05) |
| Raw p-value (`p.value`, context only) | 0.031 |
| Source | `inputs/Proteomic_data .xlsx`, sheet `Tumor vs Normal` |

**Key interpretation:** The raw p-value (0.031) falls below 0.05 and the sheet's `is_sig` flag is `True`, but that flag reflects an unadjusted criterion. The adjusted p-value of 0.226 does not survive multiple-testing correction across the 3,850 quantified proteins (only 118 pass adj.Pval ≤ 0.05). ENO1's strong up-regulation (log2FC = 2.27, FC = 4.81) remains a descriptive fold-change observation, not an FDR-controlled hit. Per the task rules, the raw value was never relabeled as adjusted.

**Artifacts written (all validated):**
- `output/eno1_significance.json` — exact schema `gene, adjusted_p_value, fdr_threshold, significant, source_file, source_sheet`; `significant=false`
- `output/analysis.py` — reproducible script generating the JSON and report from the frozen inputs
- `output/report.md` — threshold-calibrated interpretation with the raw-vs-adjusted distinction documented

**Validation:** one production run plus one independent cross-check (openpyxl re-read of the workbook without pandas) confirming the JSON's adjusted value matches `adj.Pval` exactly, differs from the raw p-value, has correct keys/types, and satisfies `significant == (adjusted_p_value ≤ fdr_threshold)` — all assertions passed.