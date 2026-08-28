"""Audit ENO1 adjusted significance from the designated proteomics workbook."""

from __future__ import annotations

import json
import math
import re
import zipfile
from pathlib import Path, PurePosixPath
from xml.etree import ElementTree as ET


GENE = "ENO1"
SOURCE_FILE = "Proteomic_data .xlsx"
SOURCE_SHEET = "Tumor vs Normal"
FDR_THRESHOLD = 0.05
MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
CELL_REF = re.compile(r"^([A-Z]+)([0-9]+)$")


def column_index(reference: str) -> int:
    match = CELL_REF.fullmatch(reference)
    if not match:
        raise ValueError(f"Invalid Excel cell reference: {reference!r}")
    value = 0
    for letter in match.group(1):
        value = value * 26 + ord(letter) - ord("A") + 1
    return value - 1


def shared_strings(archive: zipfile.ZipFile) -> list[str]:
    name = "xl/sharedStrings.xml"
    if name not in archive.namelist():
        return []
    root = ET.fromstring(archive.read(name))
    return [
        "".join(node.text or "" for node in item.iter(f"{{{MAIN_NS}}}t"))
        for item in root.findall(f"{{{MAIN_NS}}}si")
    ]


def worksheet_path(archive: zipfile.ZipFile, sheet_name: str) -> str:
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    sheets = workbook.find(f"{{{MAIN_NS}}}sheets")
    if sheets is None:
        raise ValueError("Workbook has no sheets collection")
    matches = [sheet for sheet in sheets if sheet.attrib.get("name") == sheet_name]
    if len(matches) != 1:
        raise ValueError(f"Expected one sheet named {sheet_name!r}; found {len(matches)}")
    relationship_id = matches[0].attrib.get(f"{{{REL_NS}}}id")
    relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    targets = [
        rel.attrib["Target"]
        for rel in relationships.findall(f"{{{PKG_REL_NS}}}Relationship")
        if rel.attrib.get("Id") == relationship_id
    ]
    if len(targets) != 1:
        raise ValueError("Could not resolve target worksheet relationship")
    target = PurePosixPath(targets[0].lstrip("/"))
    if target.parts and target.parts[0] == "xl":
        return str(target)
    return str(PurePosixPath("xl") / target)


def cell_value(cell: ET.Element, strings: list[str]) -> object:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return "".join(node.text or "" for node in cell.iter(f"{{{MAIN_NS}}}t"))
    value_node = cell.find(f"{{{MAIN_NS}}}v")
    if value_node is None or value_node.text is None:
        return None
    text = value_node.text
    if cell_type == "s":
        return strings[int(text)]
    if cell_type == "b":
        return text == "1"
    if cell_type in {"str", "e"}:
        return text
    number = float(text)
    return int(number) if number.is_integer() else number


def read_rows(workbook_path: Path) -> list[dict[int, object]]:
    with zipfile.ZipFile(workbook_path) as archive:
        strings = shared_strings(archive)
        sheet = ET.fromstring(archive.read(worksheet_path(archive, SOURCE_SHEET)))
        rows: list[dict[int, object]] = []
        for row in sheet.findall(f".//{{{MAIN_NS}}}sheetData/{{{MAIN_NS}}}row"):
            parsed: dict[int, object] = {}
            for cell in row.findall(f"{{{MAIN_NS}}}c"):
                parsed[column_index(cell.attrib["r"])] = cell_value(cell, strings)
            rows.append(parsed)
        return rows


def analyze(workbook_path: Path) -> tuple[dict[str, object], dict[str, object]]:
    rows = read_rows(workbook_path)
    if not rows:
        raise ValueError("Target sheet is empty")
    header = {value: index for index, value in rows[0].items() if isinstance(value, str)}
    required = ("gene", "p.value", "adj.Pval")
    missing = [name for name in required if name not in header]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    if header["p.value"] == header["adj.Pval"]:
        raise ValueError("Raw and adjusted p-values do not resolve to distinct columns")

    matches = [row for row in rows[1:] if row.get(header["gene"]) == GENE]
    if len(matches) != 1:
        raise ValueError(f"Expected one ENO1 row; found {len(matches)}")
    row = matches[0]
    raw_p = row.get(header["p.value"])
    adjusted_p = row.get(header["adj.Pval"])
    if isinstance(adjusted_p, bool) or not isinstance(adjusted_p, (int, float)):
        raise ValueError(f"ENO1 adj.Pval is not numeric: {adjusted_p!r}")
    adjusted_p = float(adjusted_p)
    if not math.isfinite(adjusted_p) or not 0 <= adjusted_p <= 1:
        raise ValueError(f"ENO1 adj.Pval is outside [0,1]: {adjusted_p!r}")
    if isinstance(raw_p, bool) or not isinstance(raw_p, (int, float)):
        raise ValueError(f"ENO1 raw p.value is not numeric: {raw_p!r}")
    raw_p = float(raw_p)
    if not math.isfinite(raw_p) or not 0 <= raw_p <= 1:
        raise ValueError(f"ENO1 raw p.value is outside [0,1]: {raw_p!r}")

    result: dict[str, object] = {
        "gene": GENE,
        "adjusted_p_value": adjusted_p,
        "fdr_threshold": FDR_THRESHOLD,
        "significant": adjusted_p <= FDR_THRESHOLD,
        "source_file": SOURCE_FILE,
        "source_sheet": SOURCE_SHEET,
    }
    audit = {
        "raw_p_value": raw_p,
        "data_rows": len(rows) - 1,
        "raw_column": "p.value",
        "adjusted_column": "adj.Pval",
    }
    return result, audit


def report(result: dict[str, object], audit: dict[str, object]) -> str:
    outcome = "significant" if result["significant"] else "not significant"
    return f"""# ENO1 adjusted-significance audit

## Result

ENO1 has adjusted p-value `{result['adjusted_p_value']}` from the `adj.Pval` column. At an FDR threshold of `{result['fdr_threshold']}`, it is **{outcome}** because `{result['adjusted_p_value']} > {result['fdr_threshold']}`.

- Source file: `{result['source_file']}`
- Source sheet: `{result['source_sheet']}`
- Gene: `{result['gene']}`
- Adjusted p-value: `{result['adjusted_p_value']}`
- FDR threshold: `{result['fdr_threshold']}`
- Significant: `{str(result['significant']).lower()}`

## Audit note

The ENO1 raw `p.value` is `{audit['raw_p_value']}`, but it is not relabeled or used as the adjusted p-value. The result is recomputed directly as `adj.Pval <= 0.05`; any pre-existing worksheet significance flag is not substituted for this threshold-calibrated decision.

The unrelated RNA/m6A workbook was not opened or used.
"""


def main() -> None:
    output_dir = Path(__file__).resolve().parent
    repo_root = Path(__file__).resolve().parents[5]
    workbook_path = repo_root / "inputs" / "ls06-eno1-significance-audit" / SOURCE_FILE
    result, audit = analyze(workbook_path)
    (output_dir / "eno1_significance.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (output_dir / "report.md").write_text(report(result, audit), encoding="utf-8")


if __name__ == "__main__":
    main()
