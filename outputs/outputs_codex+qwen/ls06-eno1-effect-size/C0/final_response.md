All deliverables are complete and validated.

## Results — ENO1 tumor-vs-normal effect size

| Quantity | Value |
|---|---|
| Normal abundance | 72,896,133.29 |
| Tumor abundance | 350,385,456.45 |
| **Fold change (Tumor/Normal)** | **4.8066** |
| **log2 fold change** | **+2.2650** |

**Direction: ENO1 is upregulated in tumor** (~4.81-fold higher in tumor than normal; sheet context: p.value = 0.031, adj.Pval = 0.226).

## Deliverables (all in `output/`)
- **`eno1_effect.json`** — exact schema `gene, tumor_value, normal_value, fold_change, log2_fold_change, source_file, source_sheet`; source is `Proteomic_data .xlsx`, sheet `Tumor vs Normal`.
- **`analysis.py`** — reads the proteomics workbook's `Tumor vs Normal` sheet, isolates the single ENO1 row (1 of 3,850 rows), computes FC and log2FC, and asserts consistency against the sheet's own rounded FC/log2FC columns.
- **`report.md`** — values, direction statement, and validation notes.

## Validation
- Primary run's recomputed FC (4.8066) and log2FC (2.2650) match the workbook's recorded FC = 4.81 / log2FC = 2.27 after rounding.
- Independent recompute via a separate openpyxl code path reproduced identical values; JSON key order, types, and direction sign all passed schema assertions.
- The decoy workbook `MeRIP_RNA_result.xlsx` was never opened or substituted.