All deliverables are complete and validated.

**Result — ENO1 (from `Proteomic_data .xlsx`, sheet `Tumor vs Normal`):**
- Tumor abundance: 350,385,456.45 · Normal abundance: 72,896,133.29
- **Fold change (Tumor/Normal) = 4.8066**, **log2FC = +2.2650**
- **Direction: ENO1 is upregulated in tumor** (FC > 1, positive log2FC)

**Validation:** the sheet's own pre-computed columns for the ENO1 row (FC = 4.81, log2FC = 2.27, stored rounded to 2 decimals) agree with the values recomputed from the raw abundances. The decoy workbook `MeRIP_RNA_result.xlsx` was not used, per the input contract.

**Artifacts written:**
- `output/eno1_effect.json` — exact keys `gene, tumor_value, normal_value, fold_change, log2_fold_change, source_file, source_sheet` (schema assertion passed)
- `output/analysis.py` — reproducible script that regenerates the JSON
- `output/report.md` — method, direction statement, and cross-check details