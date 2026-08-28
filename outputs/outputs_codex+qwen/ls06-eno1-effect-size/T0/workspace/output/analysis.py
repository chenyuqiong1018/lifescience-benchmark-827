"""Compute ENO1 tumor-vs-normal effect size from the proteomics workbook.

Reads only the analysis-target workbook `Proteomic_data .xlsx`
(sheet `Tumor vs Normal`); the unrelated `MeRIP_RNA_result.xlsx`
transcriptomics workbook is never used.
"""

import json
import math
from pathlib import Path

import openpyxl

WORKSPACE = Path(__file__).resolve().parent.parent
SOURCE_FILE = "Proteomic_data .xlsx"
SOURCE_SHEET = "Tumor vs Normal"
GENE = "ENO1"
OUT_JSON = Path(__file__).resolve().parent / "eno1_effect.json"


def main() -> None:
    wb = openpyxl.load_workbook(WORKSPACE / "inputs" / SOURCE_FILE, read_only=True, data_only=True)
    ws = wb[SOURCE_SHEET]
    rows = ws.iter_rows(values_only=True)
    header = [str(h) for h in next(rows)]
    idx = {name: i for i, name in enumerate(header)}
    required = ["gene", "Normal", "Tumor", "FC", "log2FC"]
    missing = [c for c in required if c not in idx]
    if missing:
        raise RuntimeError(f"Missing expected columns: {missing}")

    eno1 = None
    for row in rows:
        if row[idx["gene"]] == GENE:
            eno1 = row
            break
    if eno1 is None:
        raise RuntimeError(f"{GENE} not found in sheet {SOURCE_SHEET!r}")

    normal_value = float(eno1[idx["Normal"]])
    tumor_value = float(eno1[idx["Tumor"]])

    # Compute fold change (tumor relative to normal) and log2 fold change.
    fold_change = tumor_value / normal_value
    log2_fold_change = math.log2(fold_change)

    # Cross-check against the workbook's own rounded FC / log2FC columns.
    wb_fc = float(eno1[idx["FC"]])
    wb_log2fc = float(eno1[idx["log2FC"]])
    assert abs(fold_change - wb_fc) < 0.01, (fold_change, wb_fc)
    assert abs(log2_fold_change - wb_log2fc) < 0.01, (log2_fold_change, wb_log2fc)

    result = {
        "gene": GENE,
        "tumor_value": tumor_value,
        "normal_value": normal_value,
        "fold_change": round(fold_change, 4),
        "log2_fold_change": round(log2_fold_change, 4),
        "source_file": SOURCE_FILE,
        "source_sheet": SOURCE_SHEET,
    }
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    direction = "upregulated in tumor vs normal" if fold_change > 1 else (
        "downregulated in tumor vs normal" if fold_change < 1 else "no change")
    print(json.dumps(result, indent=2))
    print(f"Direction: ENO1 is {direction} (FC = {fold_change:.4f}, "
          f"log2FC = {log2_fold_change:+.4f})")


if __name__ == "__main__":
    main()
