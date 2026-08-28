# ENO1 FDR significance audit

## Threshold-calibrated result

ENO1's adjusted p-value is `0.226`. It **does not pass** the prespecified FDR `0.05` threshold, so `significant` is `false`.

- Gene: `ENO1`
- Adjusted column used: `adj.Pval`
- Adjusted p-value: `0.226`
- Decision rule: `adjusted_p_value <= 0.05`
- Source: `Proteomic_data .xlsx`, sheet `Tumor vs Normal`

## Statistical interpretation

The raw `p.value` in the ENO1 row is `0.031`, which is below 0.05, but it is not the requested multiplicity-adjusted value. The adjusted value `0.226` is above 0.05 and therefore does not meet the requested FDR criterion. This threshold result is not proof of no biological effect; it means only that ENO1 is not significant under the specified adjusted-p-value rule.

ENO1 occurs exactly once (worksheet row 474) in the target sheet. No new hypothesis test was invented, and no worksheet flag was substituted for the explicit FDR comparison. The unrelated RNA/m6A workbook was not opened or used.
