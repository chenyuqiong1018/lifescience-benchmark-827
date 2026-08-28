Done. All deliverables are in `output/` and schema-validated.

**Result: ENO1 is NOT significant at FDR 0.05.**

- From `Proteomic_data .xlsx`, sheet `Tumor vs Normal` (3,850 protein rows), the unique ENO1 row (P06733, Alpha-enolase) has:
  - **adjusted p-value (`adj.Pval`) = 0.226** — the reported statistic
  - raw p-value (`p.value`) = 0.031 — reference only, deliberately not relabeled as adjusted (guards in the script verify adj.Pval ≠ raw and adj.Pval ≥ raw, as BH requires)
- 0.226 > 0.05 → `significant: false`. Direction of effect (tumor ↑, FC 4.81) does not survive correction; notably **no** protein in the sheet passes FDR 0.05 (min adj.Pval ≈ 0.127). The report also flags that the sheet's own `is_sig=True` for ENO1 is driven by the raw p-value, not the adjusted one.

**Deliverables:**
- `output/eno1_significance.json` — exact keys `gene, adjusted_p_value, fdr_threshold, significant, source_file, source_sheet` (schema assertion passed)
- `output/analysis.py` — extraction script with mislabeling guards
- `output/report.md` — threshold-calibrated interpretation + provenance
- `INPUT_MANIFEST.sha256.tsv` — task-local integrity manifest per the input README contract