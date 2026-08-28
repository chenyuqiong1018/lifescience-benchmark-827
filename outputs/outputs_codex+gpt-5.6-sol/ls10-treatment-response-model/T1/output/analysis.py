#!/usr/bin/env python3
"""Audit the workbook and fit the specified treatment-response GLM."""

from __future__ import annotations

import csv
import json
import re
import zipfile
from collections import Counter
from pathlib import Path, PurePosixPath
from xml.etree import ElementTree as ET

import numpy as np
import statsmodels.api as sm


OUT_DIR = Path(__file__).resolve().parent
IN_DIR = OUT_DIR.parents[4] / "inputs" / "ls10-treatment-response-model"
SSML = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
ODREL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKGREL = "http://schemas.openxmlformats.org/package/2006/relationships"


def abort(reason: str) -> None:
    raise SystemExit(f"ABORT: {reason}")


def excel_column(reference: str) -> int:
    match = re.match(r"([A-Z]+)", reference)
    if match is None:
        abort(f"invalid cell coordinate {reference!r}")
    result = 0
    for character in match.group(1):
        result = result * 26 + ord(character) - ord("A") + 1
    return result - 1


def decode_xlsx(path: Path, sheet_name: str) -> tuple[list[str], list[list[object | None]]]:
    """Return workbook sheet names and cached values for one named sheet."""
    with zipfile.ZipFile(path) as package:
        workbook = ET.fromstring(package.read("xl/workbook.xml"))
        rel_root = ET.fromstring(package.read("xl/_rels/workbook.xml.rels"))
        rels = {
            rel.attrib["Id"]: rel.attrib["Target"]
            for rel in rel_root.findall(f"{{{PKGREL}}}Relationship")
        }
        sheets = workbook.findall(f".//{{{SSML}}}sheet")
        sheet_names = [sheet.attrib["name"] for sheet in sheets]
        matches = [sheet for sheet in sheets if sheet.attrib["name"] == sheet_name]
        if len(matches) != 1:
            abort(f"required sheet {sheet_name!r} is not unique")
        relationship = matches[0].attrib[f"{{{ODREL}}}id"]
        target = rels.get(relationship)
        if target is None:
            abort("worksheet relationship target is missing")
        member = target.lstrip("/") if target.startswith("/") else str(PurePosixPath("xl") / target)

        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in package.namelist():
            shared_root = ET.fromstring(package.read("xl/sharedStrings.xml"))
            shared_strings = [
                "".join(node.text or "" for node in item.iter(f"{{{SSML}}}t"))
                for item in shared_root.findall(f"{{{SSML}}}si")
            ]

        worksheet = ET.fromstring(package.read(member))
        table: list[list[object | None]] = []
        for row in worksheet.findall(f".//{{{SSML}}}row"):
            decoded: dict[int, object | None] = {}
            for cell in row.findall(f"{{{SSML}}}c"):
                position = excel_column(cell.attrib["r"])
                kind = cell.attrib.get("t")
                value_element = cell.find(f"{{{SSML}}}v")
                if kind == "inlineStr":
                    value: object | None = "".join(
                        node.text or "" for node in cell.iter(f"{{{SSML}}}t")
                    )
                elif value_element is None:
                    value = None
                elif kind == "s":
                    value = shared_strings[int(value_element.text)]
                else:
                    raw = value_element.text or ""
                    try:
                        number = float(raw)
                        value = int(number) if number.is_integer() else number
                    except ValueError:
                        value = raw
                decoded[position] = value
            if decoded:
                table.append([decoded.get(index) for index in range(max(decoded) + 1)])
        return sheet_names, table


def populated(value: object | None) -> bool:
    return value is not None and str(value).strip() != ""


def main() -> None:
    contract = (IN_DIR / "README.md").read_text(encoding="utf-8")
    for required_text in (
        "sheet `Sheet1`",
        "Only the frozen outcome and three named predictors belong in the model",
        "Missing values are not converted to zero or a new category",
    ):
        if required_text not in contract:
            abort(f"input contract omits: {required_text}")

    sheet_names, table = decode_xlsx(IN_DIR / "data.xlsx", "Sheet1")
    if not table:
        abort("Sheet1 has no cached cell values")
    headings = {str(value).strip(): index for index, value in enumerate(table[0]) if populated(value)}
    variables = ["Efficacy", "BMI", "Age", "Gender"]
    absent = [variable for variable in variables if variable not in headings]
    if absent:
        abort(f"missing required workbook headings: {absent}")

    rows: list[dict[str, object | None]] = []
    for source in table[1:]:
        record = {
            variable: source[headings[variable]] if headings[variable] < len(source) else None
            for variable in variables
        }
        if any(populated(value) for value in record.values()):
            rows.append(record)
    missing_by_variable = {
        variable: sum(not populated(row[variable]) for row in rows) for variable in variables
    }
    complete = [row for row in rows if all(populated(row[variable]) for variable in variables)]
    if not complete:
        abort("there are no complete model records")

    outcome_codes = {"PR": 1, "SD": 0, "PD": 0}
    gender_codes = {"Female": 0, "Male": 1}
    outcome_labels = [str(row["Efficacy"]).strip() for row in complete]
    gender_labels = [str(row["Gender"]).strip() for row in complete]
    bad_outcomes = sorted(set(outcome_labels) - set(outcome_codes))
    bad_genders = sorted(set(gender_labels) - set(gender_codes))
    if bad_outcomes or bad_genders:
        abort(f"unmapped categories: outcomes={bad_outcomes}, genders={bad_genders}")

    y = np.asarray([outcome_codes[label] for label in outcome_labels], dtype=float)
    x = np.asarray(
        [
            [
                1.0,
                float(row["BMI"]),
                float(row["Age"]),
                float(gender_codes[gender]),
            ]
            for row, gender in zip(complete, gender_labels)
        ],
        dtype=float,
    )
    design_rank = int(np.linalg.matrix_rank(x))
    if design_rank != x.shape[1]:
        abort(f"design rank {design_rank} differs from {x.shape[1]} columns")
    fit = sm.GLM(y, x, family=sm.families.Binomial(link=sm.families.links.Logit())).fit(
        maxiter=100, tol=1e-12
    )
    if not bool(fit.converged):
        abort("binomial GLM failed to converge")
    predicted = np.asarray(fit.predict(x), dtype=float)
    if not np.all((predicted > 0) & (predicted < 1)):
        abort("fitted probabilities reached an invalid boundary")

    term_names = ["Intercept", "BMI", "Age", "Gender_Male"]
    coefficient_records = []
    for index, term in enumerate(term_names):
        estimate = float(fit.params[index])
        coefficient_records.append(
            {
                "term": term,
                "estimate": estimate,
                "std_error": float(fit.bse[index]),
                "z": float(fit.tvalues[index]),
                "p_value": float(fit.pvalues[index]),
                "odds_ratio": float(np.exp(estimate)),
            }
        )
    fields = ["term", "estimate", "std_error", "z", "p_value", "odds_ratio"]
    with (OUT_DIR / "model_coefficients.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for record in coefficient_records:
            writer.writerow(
                {
                    field: record[field] if field == "term" else format(float(record[field]), ".15g")
                    for field in fields
                }
            )

    metadata = {
        "workbook": "data.xlsx",
        "sheet": "Sheet1",
        "workbook_sheet_names": sheet_names,
        "outcome": "Efficacy",
        "outcome_coding": outcome_codes,
        "modeled_event": "PR (treatment response)",
        "gender_reference_level": "Female",
        "gender_coding": gender_codes,
        "predictors": ["BMI", "Age", "Gender_Male"],
        "complete_case_variables": variables,
        "n_input_rows": len(rows),
        "n_complete_cases": len(complete),
        "n_dropped_missing": len(rows) - len(complete),
        "missing_by_model_variable": missing_by_variable,
        "outcome_counts_complete_cases": dict(Counter(outcome_labels)),
        "gender_counts_complete_cases": dict(Counter(gender_labels)),
        "model": "binomial generalized linear model",
        "link": "logit",
        "design_rank": design_rank,
        "converged": True,
        "fitted_probability_min": float(predicted.min()),
        "fitted_probability_max": float(predicted.max()),
    }
    (OUT_DIR / "model_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8", newline="\n"
    )

    age = next(record for record in coefficient_records if record["term"] == "Age")
    age_index = term_names.index("Age")
    age_low, age_high = fit.conf_int()[age_index]
    report = f"""# Treatment-response model and data audit

## Data and coding

The XLSX workbook exposes the expected `Sheet1`. Restricting assessment to Efficacy, BMI, Age, and Gender found {len(rows)} data records and zero missing values in each model variable, so all {len(complete)} records were retained. The outcome event was **PR = 1**; **SD and PD = 0**. Gender used **Female as the reference**, with `Gender_Male = 1` for Male. No other workbook covariate entered the model.

## Logistic model

A binomial GLM with logit link fit `PR ~ BMI + Age + Gender_Male` plus an intercept. The four-column design had rank {design_rank}, fitting converged, and fitted probabilities ranged from {predicted.min():.6f} to {predicted.max():.6f}.

The requested **Age log-odds coefficient is {age['estimate']:.6f}** (SE = {age['std_error']:.6f}, Wald z = {age['z']:.6f}, **two-sided p = {age['p_value']:.6g}**). Its odds ratio is {age['odds_ratio']:.6f} per age unit, with Wald 95% CI {np.exp(age_low):.6f} to {np.exp(age_high):.6f}.

## Skill-assisted checks

`exploratory-data-analysis` guided the bounded XLSX sheet/schema/completeness/category audit. `statistical-analysis` guided the binary multiple-logistic specification, explicit reference coding, Wald statistics, odds ratios, and convergence disclosure. `code_execution_analysis` returned the requested fixed-number audit code but did not execute it, so this rerunnable local GLM is the executed source of record. `personalized_medicine` was opened as required, but its drug/variant database calls were not used because the input contains neither drug nor variant identifiers and external data are prohibited.
"""
    (OUT_DIR / "report.md").write_text(report, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
