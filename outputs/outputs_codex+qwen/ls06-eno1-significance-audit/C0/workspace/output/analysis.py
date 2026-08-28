#!/usr/bin/env python3
"""Retrieve ENO1's adjusted p-value from the proteomics workbook and
interpret it against an FDR 0.05 threshold.

Primary method: read sheet 'Tumor vs Normal' of 'Proteomic_data .xlsx',
select the ENO1 row, and take the column 'adj.Pval' (the multiple-testing
corrected p-value). The raw p-value column 'p.value' is intentionally NOT
used for the significance call and is never relabeled as adjusted.
Writes output/eno1_significance.json.
"""
import json
from pathlib import Path

import pandas as pd

SOURCE_FILE = "Proteomic_data .xlsx"
SOURCE_SHEET = "Tumor vs Normal"
RAW_P_COL = "p.value"
ADJ_P_COL = "adj.Pval"
GENE = "ENO1"
FDR_THRESHOLD = 0.05

ROOT = Path(__file__).resolve().parent
INPUT = ROOT.parent / "inputs" / SOURCE_FILE


def main() -> None:
    df = pd.read_excel(INPUT, sheet_name=SOURCE_SHEET)

    required = {"gene", RAW_P_COL, ADJ_P_COL}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing expected columns: {sorted(missing)}")

    rows = df[df["gene"] == GENE]
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one {GENE} row, found {len(rows)}")

    raw_p = float(rows[RAW_P_COL].iloc[0])
    adj_p = float(rows[ADJ_P_COL].iloc[0])

    # Sanity: p-values must be valid probabilities and distinct columns.
    for name, v in ((RAW_P_COL, raw_p), (ADJ_P_COL, adj_p)):
        if not 0.0 <= v <= 1.0:
            raise ValueError(f"{name} out of [0,1]: {v}")
    if raw_p == adj_p == 0.0:
        raise ValueError("Suspicious all-zero p-values")

    result = {
        "gene": GENE,
        "adjusted_p_value": adj_p,
        "fdr_threshold": FDR_THRESHOLD,
        "significant": bool(adj_p < FDR_THRESHOLD),
        "source_file": SOURCE_FILE,
        "source_sheet": SOURCE_SHEET,
    }

    out = ROOT / "eno1_significance.json"
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    print(f"gene:                {GENE}")
    print(f"raw p.value (info):  {raw_p}")
    print(f"adjusted p (adj.Pval): {adj_p}")
    print(f"FDR threshold:       {FDR_THRESHOLD}")
    print(f"significant:         {result['significant']}")
    print(f"wrote: {out}")


if __name__ == "__main__":
    main()
