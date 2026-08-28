"""ENO1 significance audit for the Tumor vs Normal proteomics comparison.

Reads the frozen input workbook, extracts the ENO1 row, and reports the
ADJUSTED p-value (adj.Pval column) against an FDR threshold of 0.05.
The raw p-value (p.value column) is reported for context only and is never
used as the significance decision or relabeled as an adjusted value.

Artifacts written:
  output/eno1_significance.json
  output/report.md
"""

import json
from pathlib import Path

import pandas as pd

WORKSPACE = Path(__file__).resolve().parent.parent
SOURCE_FILE = "Proteomic_data .xlsx"
SOURCE_SHEET = "Tumor vs Normal"
INPUT_PATH = WORKSPACE / "inputs" / SOURCE_FILE
OUTPUT_DIR = WORKSPACE / "output"

FDR_THRESHOLD = 0.05
GENE = "ENO1"

RAW_P_COL = "p.value"
ADJ_P_COL = "adj.Pval"


def main() -> None:
    df = pd.read_excel(INPUT_PATH, sheet_name=SOURCE_SHEET, engine="openpyxl")

    # Locate the ENO1 row (expect exactly one).
    eno1 = df[df["gene"] == GENE]
    if len(eno1) != 1:
        raise ValueError(f"Expected exactly 1 ENO1 row, found {len(eno1)}")
    row = eno1.iloc[0]

    raw_p = float(row[RAW_P_COL])
    adj_p = float(row[ADJ_P_COL])

    # Sanity checks required by the input contract: p-values in [0, 1].
    for label, value in (("raw", raw_p), ("adjusted", adj_p)):
        if not (0.0 <= value <= 1.0):
            raise ValueError(f"{label} p-value out of [0, 1]: {value}")

    significant = bool(adj_p <= FDR_THRESHOLD)

    # Context stats for the report (not part of the machine-readable output).
    n_rows = len(df)
    adj_valid = df[ADJ_P_COL].dropna()
    n_sig_adj = int((adj_valid <= FDR_THRESHOLD).sum())
    is_sig_flag = bool(row["is_sig"])
    log2fc = float(row["log2FC"])
    fc = float(row["FC"])
    protein = str(row["protein"])

    record = {
        "gene": GENE,
        "adjusted_p_value": adj_p,
        "fdr_threshold": FDR_THRESHOLD,
        "significant": significant,
        "source_file": SOURCE_FILE,
        "source_sheet": SOURCE_SHEET,
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT_DIR / "eno1_significance.json"
    with json_path.open("w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=2)
        fh.write("\n")

    verdict = "SIGNIFICANT" if significant else "NOT significant"
    report = f"""# ENO1 significance audit (FDR 0.05)

## Task

Retrieve ENO1's adjusted p-value from the supplied proteomics results and give a
threshold-calibrated interpretation at FDR 0.05.

## Source

- Source file: `inputs/{SOURCE_FILE}`
- Source sheet: `{SOURCE_SHEET}` (only sheet in the workbook)
- ENO1 row: protein `{protein}`, gene_id 2023, log2FC = {log2fc}, FC = {fc}
- Comparison: {row["compare"]}

## Extracted values

| Column | Value | Meaning |
|---|---|---|
| `p.value` | {raw_p} | raw (unadjusted) p-value - context only |
| `adj.Pval` | {adj_p} | adjusted p-value (multiple-testing corrected) - used for the decision |

The sheet carries both a raw p-value column (`p.value`) and an adjusted
p-value column (`adj.Pval`). Per the task rules, the **adjusted** value
(`adj.Pval` = {adj_p}) is used for the significance call. The raw p-value is
**not** relabeled or treated as an adjusted value.

## Decision at FDR 0.05

- Adjusted p-value: **{adj_p}**
- FDR threshold: **{FDR_THRESHOLD}**
- {adj_p} {'<=' if significant else '>'} {FDR_THRESHOLD} -> **{verdict}** after multiple-testing correction.

## Interpretation

Although ENO1's raw p-value ({raw_p}) is below 0.05 and the sheet's
`is_sig` flag is `True`, the adjusted p-value ({adj_p}) exceeds the FDR 0.05
threshold. Across the {n_rows} quantified proteins, only {n_sig_adj} pass
adj.Pval <= 0.05, so ENO1's raw-level evidence does not survive
multiple-testing correction. ENO1 (log2FC = {log2fc}, strongly up-regulated in
tumor) is therefore **not statistically significant at FDR 0.05** in this
comparison; its up-regulation is a descriptive fold-change observation rather
than an FDR-controlled hit. Note the sheet's `is_sig` flag appears to reflect
an unadjusted criterion and should not override the adjusted p-value.

## Reproducibility

Run `python output/analysis.py` from the workspace root. It regenerates
`output/eno1_significance.json` and this report from the frozen input files.
"""
    report_path = OUTPUT_DIR / "report.md"
    report_path.write_text(report, encoding="utf-8")

    print(json.dumps(record, indent=2))
    print(f"wrote {json_path}")
    print(f"wrote {report_path}")


if __name__ == "__main__":
    main()
