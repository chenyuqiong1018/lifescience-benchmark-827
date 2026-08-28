"""Compute ENO1 tumor-vs-normal fold change from the proteomics workbook.

Primary method: read sheet "Tumor vs Normal" from inputs/Proteomic_data .xlsx,
take the ENO1 row, and compute fold_change = Tumor / Normal and
log2_fold_change = log2(fold_change). The sheet's own rounded FC / log2FC
columns are used only as an internal consistency check.

The MeRIP_RNA_result.xlsx workbook is an unrelated transcript/m6A decoy and is
intentionally never opened here.
"""

import json
import math
from pathlib import Path

import pandas as pd

WORKSPACE = Path(__file__).resolve().parent.parent
SOURCE_FILE = "Proteomic_data .xlsx"
SOURCE_SHEET = "Tumor vs Normal"
SOURCE_PATH = WORKSPACE / "inputs" / SOURCE_FILE
OUT_DIR = WORKSPACE / "output"


def main() -> None:
    df = pd.read_excel(SOURCE_PATH, sheet_name=SOURCE_SHEET)

    required = ["gene", "Normal", "Tumor", "FC", "log2FC"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing expected columns: {missing}")

    eno1 = df[df["gene"].astype(str).str.strip() == "ENO1"]
    if len(eno1) != 1:
        raise SystemExit(f"Expected exactly 1 ENO1 row, found {len(eno1)}")
    row = eno1.iloc[0]

    normal = float(row["Normal"])
    tumor = float(row["Tumor"])
    if normal <= 0 or tumor <= 0:
        raise SystemExit("Non-positive abundance value; fold change undefined")

    fold_change = tumor / normal
    log2_fold_change = math.log2(fold_change)

    # Internal consistency check against the sheet's own rounded columns.
    sheet_fc = float(row["FC"])
    sheet_log2fc = float(row["log2FC"])
    assert abs(round(fold_change, 2) - sheet_fc) < 0.011, (fold_change, sheet_fc)
    assert abs(round(log2_fold_change, 2) - sheet_log2fc) < 0.011, (
        log2_fold_change,
        sheet_log2fc,
    )

    direction = (
        "upregulated in tumor"
        if fold_change > 1
        else ("downregulated in tumor" if fold_change < 1 else "unchanged")
    )

    result = {
        "gene": "ENO1",
        "tumor_value": tumor,
        "normal_value": normal,
        "fold_change": fold_change,
        "log2_fold_change": log2_fold_change,
        "source_file": SOURCE_FILE,
        "source_sheet": SOURCE_SHEET,
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_DIR / "eno1_effect.json", "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=2)
        fh.write("\n")

    report = f"""# ENO1 tumor-versus-normal effect size (proteomics)

## Source
- File: `inputs/{SOURCE_FILE}` (sheet `{SOURCE_SHEET}`)
- The unrelated workbook `MeRIP_RNA_result.xlsx` (transcript/m6A decoy) was **not** used.

## ENO1 values
| Quantity | Value |
|---|---|
| Normal abundance | {normal:,.2f} |
| Tumor abundance | {tumor:,.2f} |
| Fold change (Tumor / Normal) | {fold_change:.4f} |
| log2 fold change | {log2_fold_change:.4f} |

## Direction
ENO1 is **{direction}**: tumor abundance is ~{fold_change:.2f}-fold the normal
abundance (log2FC = +{log2_fold_change:.2f}).

## Validation
- Exactly one ENO1 row in the sheet (of {len(df)} data rows).
- Recomputed FC {fold_change:.4f} matches the sheet's recorded FC ({sheet_fc}) after rounding.
- Recomputed log2FC {log2_fold_change:.4f} matches the sheet's recorded log2FC ({sheet_log2fc}) after rounding.
- Sheet-level statistics for context: p.value = {row['p.value']}, adj.Pval = {row['adj.Pval']}, is_sig = {row['is_sig']}.
"""
    with open(OUT_DIR / "report.md", "w", encoding="utf-8") as fh:
        fh.write(report)

    print(json.dumps(result, indent=2))
    print("Direction:", direction)


if __name__ == "__main__":
    main()
