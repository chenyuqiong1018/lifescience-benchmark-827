# ENO1 adjusted-significance audit

## Result

ENO1 has adjusted p-value `0.226` from the `adj.Pval` column. At an FDR threshold of `0.05`, it is **not significant** because `0.226 > 0.05`.

- Source file: `Proteomic_data .xlsx`
- Source sheet: `Tumor vs Normal`
- Gene: `ENO1`
- Adjusted p-value: `0.226`
- FDR threshold: `0.05`
- Significant: `false`

## Audit note

The ENO1 raw `p.value` is `0.031`, but it is not relabeled or used as the adjusted p-value. The result is recomputed directly as `adj.Pval <= 0.05`; any pre-existing worksheet significance flag is not substituted for this threshold-calibrated decision.

The unrelated RNA/m6A workbook was not opened or used.
