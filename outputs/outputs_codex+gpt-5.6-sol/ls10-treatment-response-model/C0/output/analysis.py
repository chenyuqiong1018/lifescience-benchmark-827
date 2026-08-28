#!/usr/bin/env python3
"""Fit the prespecified treatment-response logistic regression."""

from __future__ import annotations

import csv
import json
import math
import re
import zipfile
from collections import Counter
from pathlib import Path, PurePosixPath
from xml.etree import ElementTree as ET

import numpy as np
import statsmodels.api as sm


HERE = Path(__file__).resolve().parent
TASK_INPUT = HERE.parents[4] / "inputs" / "ls10-treatment-response-model"
MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
DOC_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
MODEL_COLUMNS = ["Efficacy", "BMI", "Age", "Gender"]
OUTCOME_CODING = {"PR": 1, "SD": 0, "PD": 0}
GENDER_CODING = {"Female": 0, "Male": 1}


def abort(message: str) -> None:
    raise SystemExit(f"ABORT: {message}")


def column_number(cell_reference: str) -> int:
    match = re.match(r"([A-Z]+)", cell_reference)
    if not match:
        abort(f"invalid cell reference {cell_reference!r}")
    value = 0
    for letter in match.group(1):
        value = value * 26 + ord(letter) - ord("A") + 1
    return value - 1


def read_sheet(path: Path, sheet_name: str) -> list[list[object | None]]:
    """Read cached XLSX cell values without modifying the workbook."""
    with zipfile.ZipFile(path) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        relation_targets = {
            item.attrib["Id"]: item.attrib["Target"]
            for item in relationships.findall(f"{{{PKG_REL_NS}}}Relationship")
        }
        selected = None
        sheets = workbook.find(f"{{{MAIN_NS}}}sheets")
        if sheets is None:
            abort("workbook has no sheets collection")
        for sheet in sheets:
            if sheet.attrib.get("name") == sheet_name:
                selected = sheet
                break
        if selected is None:
            abort(f"missing required sheet {sheet_name!r}")
        relationship_id = selected.attrib.get(f"{{{DOC_REL_NS}}}id")
        if relationship_id not in relation_targets:
            abort("sheet relationship is unresolved")
        target = relation_targets[relationship_id]
        sheet_path = (
            target.lstrip("/")
            if target.startswith("/")
            else str(PurePosixPath("xl") / target)
        )

        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in shared.findall(f"{{{MAIN_NS}}}si"):
                shared_strings.append(
                    "".join(node.text or "" for node in item.iter(f"{{{MAIN_NS}}}t"))
                )

        worksheet = ET.fromstring(archive.read(sheet_path))
        parsed_rows: list[list[object | None]] = []
        for row in worksheet.findall(f".//{{{MAIN_NS}}}row"):
            parsed: dict[int, object | None] = {}
            for cell in row.findall(f"{{{MAIN_NS}}}c"):
                index = column_number(cell.attrib["r"])
                cell_type = cell.attrib.get("t")
                value_node = cell.find(f"{{{MAIN_NS}}}v")
                if cell_type == "inlineStr":
                    value: object | None = "".join(
                        node.text or "" for node in cell.iter(f"{{{MAIN_NS}}}t")
                    )
                elif value_node is None:
                    value = None
                elif cell_type == "s":
                    value = shared_strings[int(value_node.text)]
                elif cell_type == "b":
                    value = value_node.text == "1"
                else:
                    raw = value_node.text or ""
                    try:
                        numeric = float(raw)
                        value = int(numeric) if numeric.is_integer() else numeric
                    except ValueError:
                        value = raw
                parsed[index] = value
            if parsed:
                width = max(parsed) + 1
                parsed_rows.append([parsed.get(i) for i in range(width)])
        return parsed_rows


def present(value: object | None) -> bool:
    return value is not None and str(value).strip() != ""


def main() -> None:
    contract = (TASK_INPUT / "README.md").read_text(encoding="utf-8")
    if "Only the frozen outcome and three named predictors belong in the model" not in contract:
        abort("input contract does not freeze the model variables")
    matrix = read_sheet(TASK_INPUT / "data.xlsx", "Sheet1")
    if not matrix:
        abort("Sheet1 is empty")
    header = {str(value).strip(): index for index, value in enumerate(matrix[0]) if present(value)}
    missing_columns = [name for name in MODEL_COLUMNS if name not in header]
    if missing_columns:
        abort(f"missing model columns: {missing_columns}")

    input_rows: list[dict[str, object | None]] = []
    for raw_row in matrix[1:]:
        row = {
            name: raw_row[header[name]] if header[name] < len(raw_row) else None
            for name in MODEL_COLUMNS
        }
        if any(present(value) for value in row.values()):
            input_rows.append(row)
    complete = [row for row in input_rows if all(present(row[name]) for name in MODEL_COLUMNS)]
    dropped = len(input_rows) - len(complete)
    if not complete:
        abort("no complete cases for the specified model")

    unexpected_outcomes = sorted(
        {str(row["Efficacy"]).strip() for row in complete} - set(OUTCOME_CODING)
    )
    unexpected_genders = sorted(
        {str(row["Gender"]).strip() for row in complete} - set(GENDER_CODING)
    )
    if unexpected_outcomes:
        abort(f"unrecognized Efficacy labels: {unexpected_outcomes}")
    if unexpected_genders:
        abort(f"unrecognized Gender labels: {unexpected_genders}")

    y = np.asarray(
        [OUTCOME_CODING[str(row["Efficacy"]).strip()] for row in complete], dtype=float
    )
    bmi = np.asarray([float(row["BMI"]) for row in complete], dtype=float)
    age = np.asarray([float(row["Age"]) for row in complete], dtype=float)
    gender_male = np.asarray(
        [GENDER_CODING[str(row["Gender"]).strip()] for row in complete], dtype=float
    )
    design = np.column_stack((np.ones(len(complete)), bmi, age, gender_male))
    if np.linalg.matrix_rank(design) != design.shape[1]:
        abort("design matrix is rank deficient")

    fit = sm.Logit(y, design).fit(method="newton", maxiter=100, disp=False)
    if not bool(fit.mle_retvals.get("converged", False)):
        abort("logistic regression did not converge")
    terms = ["Intercept", "BMI", "Age", "Gender_Male"]
    rows = []
    for index, term in enumerate(terms):
        estimate = float(fit.params[index])
        standard_error = float(fit.bse[index])
        z_value = float(fit.tvalues[index])
        p_value = float(fit.pvalues[index])
        rows.append(
            {
                "term": term,
                "estimate": estimate,
                "std_error": standard_error,
                "z": z_value,
                "p_value": p_value,
                "odds_ratio": math.exp(estimate),
            }
        )

    fields = ["term", "estimate", "std_error", "z", "p_value", "odds_ratio"]
    with (HERE / "model_coefficients.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {key: row[key] if key == "term" else format(float(row[key]), ".15g") for key in fields}
            )

    outcome_counts = Counter(str(row["Efficacy"]).strip() for row in complete)
    metadata = {
        "workbook": "data.xlsx",
        "sheet": "Sheet1",
        "outcome": "Efficacy",
        "outcome_coding": OUTCOME_CODING,
        "modeled_event": "PR (treatment response)",
        "gender_reference_level": "Female",
        "gender_coding": GENDER_CODING,
        "predictors": ["BMI", "Age", "Gender_Male"],
        "complete_case_variables": MODEL_COLUMNS,
        "n_input_rows": len(input_rows),
        "n_complete_cases": len(complete),
        "n_dropped_missing": dropped,
        "outcome_counts_complete_cases": dict(outcome_counts),
        "model": "binary logistic regression",
        "link": "logit",
        "converged": True,
    }
    (HERE / "model_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8", newline="\n"
    )

    age_row = rows[terms.index("Age")]
    report = f"""# Treatment-response logistic regression

## Model definition

The binary modeled event is treatment response: **Efficacy = PR was coded 1**, while **SD and PD were coded 0**. Gender used **Female as the reference level** and `Gender_Male = 1` for Male. The logit model included only BMI, Age, and Gender_Male, plus an intercept; distractor workbook covariates were excluded.

Complete-case filtering was restricted to Efficacy, BMI, Age, and Gender. All {len(input_rows)} data rows were complete, so the fitted sample was n = {len(complete)} with {int(y.sum())} modeled responses and {len(complete) - int(y.sum())} non-responses.

## Requested age result

Holding BMI and gender fixed, the **Age log-odds coefficient is {age_row['estimate']:.6f}** (SE = {age_row['std_error']:.6f}, z = {age_row['z']:.6f}, **two-sided p = {age_row['p_value']:.6g}**). Its odds ratio is {age_row['odds_ratio']:.6f} per one-unit increase in age.

The model converged using maximum-likelihood binary logistic regression with a logit link. Coefficients, standard errors, Wald z statistics, two-sided p-values, and exponentiated coefficients are provided in `model_coefficients.csv`.
"""
    (HERE / "report.md").write_text(report, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
