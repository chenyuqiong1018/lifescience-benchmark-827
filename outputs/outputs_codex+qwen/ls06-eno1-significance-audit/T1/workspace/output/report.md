# ENO1 Significance Audit — Tumor vs Normal Proteomics

## Question

Is ENO1 differentially abundant in the supplied proteomics results after
multiple-testing correction, at FDR 0.05?

## Data source

| Item | Value |
|---|---|
| Source file | `Proteomic_data .xlsx` |
| Sheet | `Tumor vs Normal` |
| Rows / unique genes | 3,850 rows / 3,746 unique gene symbols |
| ENO1 row | protein `P06733` (Alpha-enolase), gene_id 2023 |

The sheet carries two p-value columns: `p.value` (raw, uncorrected) and
`adj.Pval` (Benjamini–Hochberg-style FDR-adjusted). Per the input contract,
`MeRIP_RNA_result.xlsx` is an unrelated transcript/m6A decoy and was not used.

## Result

| Quantity | Value |
|---|---|
| ENO1 raw p-value (`p.value`) | 0.031 (reference only, **not** the reported statistic) |
| ENO1 adjusted p-value (`adj.Pval`) | **0.226** |
| FDR threshold | 0.05 |
| Significant at FDR 0.05? | **No** |

Direction of effect (context only): Tumor abundance ≈ 3.50e8 vs Normal
≈ 7.29e7 (FC 4.81, log2FC 2.27) — ENO1 is higher in tumor, but this effect
does not survive multiple-testing correction.

## Interpretation (threshold-calibrated)

- The adjusted p-value for ENO1 is 0.226, i.e. 0.226 > 0.05. At an FDR
  threshold of 0.05, ENO1 is **not** a statistically significant hit. The raw
  p-value of 0.031 would be "significant" only if no multiple-testing
  correction were applied; with 3,850 proteins tested simultaneously, that
  uncorrected comparison is not valid, and this report deliberately does not
  relabel the raw p-value as adjusted.
- The adjustment inflates ENO1's p-value ~7.3-fold (0.031 -> 0.226), the
  expected behavior of BH/FDR correction when many hypotheses are tested and
  few are truly small.
- Whole-sheet context: the smallest `adj.Pval` in the sheet is ≈ 0.127, so
  **zero** proteins pass FDR 0.05. ENO1's non-significance is consistent with
  the global evidence pattern, not an isolated artifact.
- Caution about the sheet's own `is_sig` flag: it is `True` for ENO1 even
  though adj.Pval = 0.226. The flag is driven by the raw p-value (p < 0.05),
  not by the adjusted p-value, and should not be read as an FDR-controlled
  significance call.

## Method

`output/analysis.py` reads the `Tumor vs Normal` sheet with pandas, selects
the unique row with gene symbol `ENO1` (case/whitespace-insensitive exact
match), takes `adj.Pval` as the reported statistic, and applies sanity guards:
the value must lie in [0, 1], must differ from the raw `p.value`, and must be
>= the raw p-value (BH-adjusted p-values cannot be smaller than their raw
counterparts). It then compares against 0.05 and writes
`output/eno1_significance.json`.

## Machine-readable artifact

`output/eno1_significance.json`:

```json
{
  "gene": "ENO1",
  "adjusted_p_value": 0.226,
  "fdr_threshold": 0.05,
  "significant": false,
  "source_file": "Proteomic_data .xlsx",
  "source_sheet": "Tumor vs Normal"
}
```

## Provenance

Workbooks retained from the official BixBench v1.5 capsule for `bix-37`
(2026-08-14); dataset card Apache-2.0, underlying experimental-data
redistribution terms not asserted — preserve upstream attribution before any
redistribution beyond this workspace.
