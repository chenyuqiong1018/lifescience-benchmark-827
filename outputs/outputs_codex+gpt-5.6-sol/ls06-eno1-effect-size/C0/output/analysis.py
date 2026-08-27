from __future__ import annotations

import json
import math
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


OUTPUT_DIR = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[5]
SOURCE_FILE = "Proteomic_data .xlsx"
SOURCE_SHEET = "Tumor vs Normal"
SOURCE_PATH = REPO / "inputs" / "ls06-eno1-effect-size" / SOURCE_FILE

MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_DOC = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
REL_PKG = "http://schemas.openxmlformats.org/package/2006/relationships"


def column_number(cell_ref: str) -> int:
    letters = re.match(r"[A-Z]+", cell_ref).group(0)
    number = 0
    for letter in letters:
        number = number * 26 + ord(letter) - 64
    return number


with zipfile.ZipFile(SOURCE_PATH) as archive:
    shared_strings: list[str] = []
    if "xl/sharedStrings.xml" in archive.namelist():
        shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        for item in shared_root.findall(f"{{{MAIN}}}si"):
            shared_strings.append("".join(node.text or "" for node in item.iter(f"{{{MAIN}}}t")))

    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    target_rel = None
    for sheet in workbook.findall(f".//{{{MAIN}}}sheet"):
        if sheet.attrib["name"] == SOURCE_SHEET:
            target_rel = sheet.attrib[f"{{{REL_DOC}}}id"]
            break
    if target_rel is None:
        raise ValueError(f"Sheet not found: {SOURCE_SHEET}")

    rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    target_xml = None
    for relation in rels.findall(f"{{{REL_PKG}}}Relationship"):
        if relation.attrib["Id"] == target_rel:
            target_xml = relation.attrib["Target"]
            break
    if target_xml is None:
        raise ValueError("Worksheet relationship not found")
    sheet_member = target_xml.lstrip("/")
    if not sheet_member.startswith("xl/"):
        sheet_member = "xl/" + sheet_member
    sheet_root = ET.fromstring(archive.read(sheet_member))


def cell_value(cell: ET.Element) -> object:
    kind = cell.attrib.get("t")
    if kind == "inlineStr":
        return "".join(node.text or "" for node in cell.iter(f"{{{MAIN}}}t"))
    value_node = cell.find(f"{{{MAIN}}}v")
    if value_node is None or value_node.text is None:
        return None
    raw = value_node.text
    if kind == "s":
        return shared_strings[int(raw)]
    if kind == "b":
        return raw == "1"
    try:
        return float(raw)
    except ValueError:
        return raw


rows: list[dict[int, object]] = []
for row in sheet_root.findall(f".//{{{MAIN}}}row"):
    values: dict[int, object] = {}
    for cell in row.findall(f"{{{MAIN}}}c"):
        values[column_number(cell.attrib["r"])] = cell_value(cell)
    rows.append(values)
if not rows:
    raise ValueError("Worksheet has no rows")

headers = {str(value): column for column, value in rows[0].items() if value is not None}
required = {"gene", "Normal", "Tumor", "Ratio", "FC", "log2FC"}
if not required.issubset(headers):
    raise ValueError(f"Missing required columns: {sorted(required - set(headers))}")

matches = [row for row in rows[1:] if row.get(headers["gene"]) == "ENO1"]
if len(matches) != 1:
    raise ValueError(f"Expected exactly one ENO1 row; found {len(matches)}")
eno1 = matches[0]
normal_value = float(eno1[headers["Normal"]])
tumor_value = float(eno1[headers["Tumor"]])
if normal_value <= 0 or tumor_value <= 0:
    raise ValueError("Fold change and log2 fold change require positive abundances")

fold_change = tumor_value / normal_value
log2_fold_change = math.log2(fold_change)

# The supplied workbook also contains rounded summary columns. Treat these as an
# internal consistency check, not as a substitute for calculating from abundances.
if round(fold_change, 2) != float(eno1[headers["Ratio"]]):
    raise ValueError("Calculated fold change does not match the workbook Ratio column")
if round(fold_change, 2) != float(eno1[headers["FC"]]):
    raise ValueError("Calculated fold change does not match the workbook FC column")
if round(log2_fold_change, 2) != float(eno1[headers["log2FC"]]):
    raise ValueError("Calculated log2 fold change does not match the workbook log2FC column")

result = {
    "gene": "ENO1",
    "tumor_value": tumor_value,
    "normal_value": normal_value,
    "fold_change": fold_change,
    "log2_fold_change": log2_fold_change,
    "source_file": SOURCE_FILE,
    "source_sheet": SOURCE_SHEET,
}
(OUTPUT_DIR / "eno1_effect.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

report = f"""# ENO1 tumor-versus-normal effect size

ENO1 is higher in tumor than normal in the supplied proteomics results.

| Gene | Normal value | Tumor value | Tumor / normal fold change | log2 fold change |
|---|---:|---:|---:|---:|
| ENO1 | {normal_value:.7f} | {tumor_value:.6f} | {fold_change:.12f} | {log2_fold_change:.12f} |

The calculation uses `Tumor / Normal` from `{SOURCE_FILE}`, sheet `{SOURCE_SHEET}`. The positive log2 fold change ({log2_fold_change:.3f}) corresponds to an approximately {fold_change:.2f}-fold increase in tumor. The workbook's `Ratio`, `FC`, and `log2FC` cells (4.81, 4.81, and 2.27) agree with the direct calculation after rounding to two decimals.

No physical unit is invented for the normalized abundance values. `MeRIP_RNA_result.xlsx` is an unrelated RNA/m6A workbook and was not used.
"""
(OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")

print(json.dumps(result, indent=2))
