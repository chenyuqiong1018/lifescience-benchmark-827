"""ENO1 significance audit.

Extracts ENO1's multiple-testing-corrected p-value (column `adj.Pval`) from
`inputs/Proteomic_data .xlsx`, sheet `Tumor vs Normal`, and evaluates it
against an FDR threshold of 0.05. The raw p-value column (`p.value`) is read
only for cross-checking and is never substituted for the adjusted value.
"""

import json
from pathlib import Path

import pandas as pd

WORKSPACE = Path(__file__).resolve().parent.parent
INPUT_FILE = WORKSPACE / "inputs" / "Proteomic_data .xlsx"
SOURCE_FILE_NAME = "Proteomic_data .xlsx"
SOURCE_SHEET = "Tumor vs Normal"
GENE = "ENO1"
FDR_THRESHOLD = 0.05
OUTPUT_JSON = Path(__file__).resolve().parent / "eno1_significance.json"


def main() -> None:
    df = pd.read_excel(INPUT_FILE, sheet_name=SOURCE_SHEET)

    # Locate the ENO1 row (exact gene-symbol match, case/whitespace tolerant).
    mask = df["gene"].astype(str).str.strip().str.upper().eq(GENE)
    hits = df[mask]
    if len(hits) != 1:
        raise ValueError(f"Expected exactly 1 ENO1 row, found {len(hits)}")
    row = hits.iloc[0]

    raw_p = float(row["p.value"])       # raw (uncorrected) p-value, reference only
    adj_p = float(row["adj.Pval"])      # BH/FDR-adjusted p-value: the reported value

    # Guards against mislabeling a raw p-value as adjusted.
    if not (0.0 <= adj_p <= 1.0):
        raise ValueError(f"adj.Pval out of [0, 1]: {adj_p}")
    if adj_p == raw_p:
        raise ValueError("adj.Pval identical to raw p.value; check column mapping")
    if adj_p < raw_p:
        # Benjamini-Hochberg adjusted p-values are >= raw p-values.
        raise ValueError("adj.Pval smaller than raw p.value; check column mapping")

    significant = adj_p <= FDR_THRESHOLD

    result = {
        "gene": GENE,
        "adjusted_p_value": adj_p,
        "fdr_threshold": FDR_THRESHOLD,
        "significant": significant,
        "source_file": SOURCE_FILE_NAME,
        "source_sheet": SOURCE_SHEET,
    }
    OUTPUT_JSON.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    # Context stats for the report (not written into the JSON artifact).
    n_genes = int(df["gene"].nunique())
    n_sig_adj = int((df["adj.Pval"] <= FDR_THRESHOLD).sum())
    is_sig_flag = bool(row["is_sig"])
    print(json.dumps(result, indent=2))
    print(f"rows={len(df)} unique_genes={n_genes}")
    print(f"ENO1 protein={row['protein']} raw_p={raw_p} adj_p={adj_p} sheet_is_sig={is_sig_flag}")
    print(f"genes passing adj.Pval <= {FDR_THRESHOLD}: {n_sig_adj}")


if __name__ == "__main__":
    main()
