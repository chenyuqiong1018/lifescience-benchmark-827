All deliverables are complete and validated.

**Result:** ENO1 is **upregulated in tumor vs normal** — tumor abundance (350,385,456.45) is **4.8066×** the normal abundance (72,896,133.29), giving **log2 fold change = +2.2650**. These computed values agree with the workbook's own rounded `FC` (4.81) and `log2FC` (2.27) columns.

**Details:**
- Source: `inputs/Proteomic_data .xlsx`, sheet `Tumor vs Normal` (ENO1 found by exact match on the `gene` column; protein P06733/alpha-enolase). The decoy `MeRIP_RNA_result.xlsx` was not used.
- `output/eno1_effect.json` — single record with exactly the requested keys: `gene, tumor_value, normal_value, fold_change, log2_fold_change, source_file, source_sheet` (schema assertion passed: key order, types, FC/log2FC consistency, positive direction).
- `output/analysis.py` — reads the proteomics workbook, computes FC = Tumor/Normal and log2FC from raw values, cross-checks against the workbook's FC/log2FC columns, and writes the JSON.
- `output/report.md` — states the direction (tumor > normal, upregulated) and documents method, cross-check, and provenance.