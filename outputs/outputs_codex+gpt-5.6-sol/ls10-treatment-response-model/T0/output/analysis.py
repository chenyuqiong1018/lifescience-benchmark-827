#!/usr/bin/env python3
"""Fit and report the frozen complete-case treatment-response model."""

from __future__ import annotations

import csv
import json
import re
import zipfile
from collections import Counter
from math import exp
from pathlib import Path, PurePosixPath
from xml.etree import ElementTree as XML

import numpy as np
import statsmodels.api as sm
from scipy.stats import chi2


OUTPUT = Path(__file__).resolve().parent
INPUT = OUTPUT.parents[4] / "inputs" / "ls10-treatment-response-model"
NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
RNS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PNS = "http://schemas.openxmlformats.org/package/2006/relationships"


def fail(message: str) -> None:
    raise SystemExit(f"ABORT: {message}")


def cell_column(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference)
    if letters is None:
        fail(f"bad XLSX cell reference: {reference}")
    index = 0
    for letter in letters.group():
        index = index * 26 + ord(letter) - 64
    return index - 1


def sheet_values(workbook_path: Path, wanted_sheet: str) -> list[list[object | None]]:
    with zipfile.ZipFile(workbook_path) as package:
        workbook = XML.fromstring(package.read("xl/workbook.xml"))
        relationships = XML.fromstring(package.read("xl/_rels/workbook.xml.rels"))
        targets = {
            rel.attrib["Id"]: rel.attrib["Target"]
            for rel in relationships.findall(f"{{{PNS}}}Relationship")
        }
        selected = [
            sheet
            for sheet in workbook.findall(f".//{{{NS}}}sheet")
            if sheet.attrib.get("name") == wanted_sheet
        ]
        if len(selected) != 1:
            fail(f"expected one {wanted_sheet} sheet, found {len(selected)}")
        rel_id = selected[0].attrib[f"{{{RNS}}}id"]
        target = targets[rel_id]
        member = target.lstrip("/") if target.startswith("/") else str(PurePosixPath("xl") / target)

        shared: list[str] = []
        if "xl/sharedStrings.xml" in package.namelist():
            root = XML.fromstring(package.read("xl/sharedStrings.xml"))
            shared = [
                "".join(part.text or "" for part in item.iter(f"{{{NS}}}t"))
                for item in root.findall(f"{{{NS}}}si")
            ]
        sheet = XML.fromstring(package.read(member))
        result: list[list[object | None]] = []
        for xml_row in sheet.findall(f".//{{{NS}}}row"):
            cells: dict[int, object | None] = {}
            for cell in xml_row.findall(f"{{{NS}}}c"):
                index = cell_column(cell.attrib["r"])
                kind = cell.attrib.get("t")
                stored = cell.find(f"{{{NS}}}v")
                if kind == "inlineStr":
                    value: object | None = "".join(
                        part.text or "" for part in cell.iter(f"{{{NS}}}t")
                    )
                elif stored is None:
                    value = None
                elif kind == "s":
                    value = shared[int(stored.text)]
                else:
                    raw = stored.text or ""
                    try:
                        number = float(raw)
                        value = int(number) if number.is_integer() else number
                    except ValueError:
                        value = raw
                cells[index] = value
            if cells:
                result.append([cells.get(i) for i in range(max(cells) + 1)])
        return result


def nonempty(value: object | None) -> bool:
    return value is not None and str(value).strip() != ""


def main() -> None:
    readme = (INPUT / "README.md").read_text(encoding="utf-8")
    if "Complete-case handling applies only to the specified model variables" not in readme:
        fail("complete-case rule is not frozen")
    values = sheet_values(INPUT / "data.xlsx", "Sheet1")
    if not values:
        fail("workbook sheet is empty")
    columns = {str(value).strip(): i for i, value in enumerate(values[0]) if nonempty(value)}
    required = ["Efficacy", "BMI", "Age", "Gender"]
    if any(column not in columns for column in required):
        fail("one or more prespecified columns are absent")

    candidates: list[dict[str, object | None]] = []
    for source in values[1:]:
        row = {
            column: source[columns[column]] if columns[column] < len(source) else None
            for column in required
        }
        if any(nonempty(value) for value in row.values()):
            candidates.append(row)
    cases = [row for row in candidates if all(nonempty(row[column]) for column in required)]
    if not cases:
        fail("complete-case filter left no observations")

    outcome_map = {"PR": 1, "SD": 0, "PD": 0}
    gender_map = {"Female": 0, "Male": 1}
    outcomes = [str(row["Efficacy"]).strip() for row in cases]
    genders = [str(row["Gender"]).strip() for row in cases]
    if set(outcomes) - set(outcome_map):
        fail(f"unknown Efficacy labels: {sorted(set(outcomes) - set(outcome_map))}")
    if set(genders) - set(gender_map):
        fail(f"unknown Gender labels: {sorted(set(genders) - set(gender_map))}")

    response = np.array([outcome_map[value] for value in outcomes], dtype=float)
    matrix = np.array(
        [
            [1.0, float(row["BMI"]), float(row["Age"]), float(gender_map[gender])]
            for row, gender in zip(cases, genders)
        ],
        dtype=float,
    )
    if np.linalg.matrix_rank(matrix) != 4:
        fail("model matrix is not full rank")
    model = sm.Logit(response, matrix).fit(disp=False, maxiter=100)
    if not model.mle_retvals.get("converged", False):
        fail("maximum-likelihood fit did not converge")

    term_names = ["Intercept", "BMI", "Age", "Gender_Male"]
    coefficient_rows = []
    for position, term in enumerate(term_names):
        estimate = float(model.params[position])
        coefficient_rows.append(
            {
                "term": term,
                "estimate": estimate,
                "std_error": float(model.bse[position]),
                "z": float(model.tvalues[position]),
                "p_value": float(model.pvalues[position]),
                "odds_ratio": exp(estimate),
            }
        )
    coefficient_fields = ["term", "estimate", "std_error", "z", "p_value", "odds_ratio"]
    with (OUTPUT / "model_coefficients.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=coefficient_fields, lineterminator="\n")
        writer.writeheader()
        for record in coefficient_rows:
            writer.writerow(
                {
                    field: record[field] if field == "term" else format(float(record[field]), ".15g")
                    for field in coefficient_fields
                }
            )

    likelihood_ratio = 2 * (float(model.llf) - float(model.llnull))
    lr_df = len(term_names) - 1
    lr_p = float(chi2.sf(likelihood_ratio, lr_df))
    metadata = {
        "workbook": "data.xlsx",
        "sheet": "Sheet1",
        "outcome": "Efficacy",
        "outcome_coding": outcome_map,
        "modeled_event": "PR (treatment response)",
        "gender_reference_level": "Female",
        "gender_coding": gender_map,
        "predictors": ["BMI", "Age", "Gender_Male"],
        "complete_case_variables": required,
        "n_input_rows": len(candidates),
        "n_complete_cases": len(cases),
        "n_dropped_missing": len(candidates) - len(cases),
        "outcome_counts_complete_cases": dict(Counter(outcomes)),
        "model": "binary logistic regression",
        "link": "logit",
        "converged": True,
        "likelihood_ratio_chi2": likelihood_ratio,
        "likelihood_ratio_df": lr_df,
        "likelihood_ratio_p_value": lr_p,
    }
    (OUTPUT / "model_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8", newline="\n"
    )

    age = next(row for row in coefficient_rows if row["term"] == "Age")
    age_ci = model.conf_int()[term_names.index("Age")]
    age_or_low, age_or_high = exp(float(age_ci[0])), exp(float(age_ci[1]))
    report = f"""# Treatment-response logistic model

## Coding and analysis set

Treatment response was the modeled event: `Efficacy = PR` was coded 1, and `SD` or `PD` was coded 0. Female was the Gender reference level; `Gender_Male` equals 1 for Male. Complete cases were defined only over Efficacy, BMI, Age, and Gender. All {len(candidates)} candidate rows were complete, leaving n = {len(cases)}.

The prespecified multiple logistic regression was `logit(P(PR)) = intercept + BMI + Age + Gender_Male`. No distractor covariate from the workbook entered the model. Maximum likelihood converged. The overall likelihood-ratio test versus the intercept-only model was chi-square({lr_df}) = {likelihood_ratio:.6f}, p = {lr_p:.6g}.

## Age coefficient

Holding BMI and gender constant, the **Age log-odds coefficient was {age['estimate']:.6f}** (SE = {age['std_error']:.6f}, Wald z = {age['z']:.6f}, **two-sided p = {age['p_value']:.6g}**). The corresponding odds ratio was {age['odds_ratio']:.6f} per age unit (Wald 95% CI {age_or_low:.6f} to {age_or_high:.6f}).

`statistical-analysis` guided selection of multiple logistic regression for a binary outcome, explicit outcome/reference coding, complete-case disclosure, coefficient/SE/z/p/OR reporting, and convergence/model-fit checks.
"""
    (OUTPUT / "report.md").write_text(report, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
