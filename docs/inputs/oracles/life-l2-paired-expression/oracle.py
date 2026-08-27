#!/usr/bin/env python3
"""Deterministic grader for life-l2-paired-expression."""

from __future__ import annotations

import argparse
import csv
import json
import math
import struct
import sys
from pathlib import Path

TASK_ID = "life-l2-paired-expression"
TOLERANCE = 1e-6
EXPECTED_COLUMNS = [
    "subject_id",
    "pre_expression",
    "post_expression",
    "paired_change",
    "direction",
]


def direction(value: float) -> str:
    if value > TOLERANCE:
        return "up"
    if value < -TOLERANCE:
        return "down"
    return "unchanged"


def load_input(path: Path) -> dict[str, tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    result: dict[str, tuple[float, float]] = {}
    for row in rows:
        subject = row["subject_id"].strip()
        if not subject or subject in result:
            raise ValueError("input contains an empty or duplicate subject_id")
        result[subject] = (float(row["pre_expression"]), float(row["post_expression"]))
    return result


def finite_float(value: object) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("non-finite numeric value")
    return number


def valid_png(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size < 33:
        return False
    data = path.read_bytes()
    if not data.startswith(b"\x89PNG\r\n\x1a\n") or data[12:16] != b"IHDR":
        return False
    width, height = struct.unpack(">II", data[16:24])
    return width > 0 and height > 0 and data[-12:] == b"\x00\x00\x00\x00IEND\xaeB\x60\x82"


def grade(workspace: Path) -> dict[str, object]:
    input_path = workspace / "inputs" / "paired_expression.csv"
    output = workspace / "output"
    results_path = output / "paired_expression_results.csv"
    summary_path = output / "summary.json"
    chart_path = output / "paired_expression.png"

    criteria = {
        "required_files_and_formats": {"points": 0, "max_points": 10},
        "subject_identity_and_pairing": {"points": 0, "max_points": 15},
        "values_preserved": {"points": 0, "max_points": 10},
        "paired_changes": {"points": 0, "max_points": 20},
        "direction_labels": {"points": 0, "max_points": 10},
        "summary": {"points": 0, "max_points": 10},
        "png_export": {"points": 0, "max_points": 5},
    }
    failures: list[str] = []

    try:
        expected = load_input(input_path)
    except Exception as exc:
        return {
            "task_id": TASK_ID,
            "deterministic_score": 0,
            "criteria": criteria,
            "hardgate_pass": False,
            "failure_codes": ["INPUT_INVALID"],
            "error": str(exc),
        }

    rows: list[dict[str, str]] = []
    summary: dict[str, object] | None = None
    csv_parseable = json_parseable = False

    try:
        with results_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != EXPECTED_COLUMNS:
                failures.append("CSV_COLUMNS_INVALID")
            else:
                rows = list(reader)
                csv_parseable = True
    except Exception:
        failures.append("RESULTS_MISSING_OR_INVALID")

    try:
        with summary_path.open(encoding="utf-8") as handle:
            loaded = json.load(handle)
        if isinstance(loaded, dict):
            summary = loaded
            json_parseable = True
        else:
            failures.append("SUMMARY_NOT_OBJECT")
    except Exception:
        failures.append("SUMMARY_MISSING_OR_INVALID")

    png_ok = valid_png(chart_path)
    if not png_ok:
        failures.append("PNG_MISSING_OR_INVALID")

    if csv_parseable and json_parseable and png_ok:
        criteria["required_files_and_formats"]["points"] = 10

    parsed: dict[str, tuple[float, float, float, str]] = {}
    row_values_ok = True
    for row in rows:
        try:
            subject = row["subject_id"].strip()
            if not subject or subject in parsed:
                raise ValueError("empty or duplicate subject")
            parsed[subject] = (
                finite_float(row["pre_expression"]),
                finite_float(row["post_expression"]),
                finite_float(row["paired_change"]),
                row["direction"].strip().lower(),
            )
        except Exception:
            row_values_ok = False
            break

    identities_ok = row_values_ok and set(parsed) == set(expected) and len(parsed) == len(expected)
    if identities_ok:
        criteria["subject_identity_and_pairing"]["points"] = 15
    else:
        failures.append("PAIRING_OR_IDENTITY_MISMATCH")

    values_ok = identities_ok and all(
        math.isclose(parsed[s][0], expected[s][0], abs_tol=TOLERANCE)
        and math.isclose(parsed[s][1], expected[s][1], abs_tol=TOLERANCE)
        for s in expected
    )
    if values_ok:
        criteria["values_preserved"]["points"] = 10
    else:
        failures.append("SOURCE_VALUES_CHANGED")

    changes_ok = values_ok and all(
        math.isclose(parsed[s][2], expected[s][1] - expected[s][0], abs_tol=TOLERANCE)
        for s in expected
    )
    if changes_ok:
        criteria["paired_changes"]["points"] = 20
    else:
        failures.append("PAIRED_CHANGE_INCORRECT")

    directions_ok = changes_ok and all(parsed[s][3] == direction(parsed[s][2]) for s in expected)
    if directions_ok:
        criteria["direction_labels"]["points"] = 10
    else:
        failures.append("DIRECTION_INCORRECT")

    expected_changes = [post - pre for pre, post in expected.values()]
    expected_directions = [direction(value) for value in expected_changes]
    summary_ok = False
    if summary is not None:
        try:
            summary_ok = (
                summary.get("task_id") == TASK_ID
                and int(summary["pair_count"]) == len(expected)
                and math.isclose(finite_float(summary["mean_paired_change"]), sum(expected_changes) / len(expected_changes), abs_tol=TOLERANCE)
                and int(summary["up_count"]) == expected_directions.count("up")
                and int(summary["down_count"]) == expected_directions.count("down")
                and int(summary["unchanged_count"]) == expected_directions.count("unchanged")
            )
        except Exception:
            summary_ok = False
    if summary_ok:
        criteria["summary"]["points"] = 10
    else:
        failures.append("SUMMARY_INCORRECT")

    if png_ok:
        criteria["png_export"]["points"] = 5

    score = sum(int(item["points"]) for item in criteria.values())
    failures = list(dict.fromkeys(failures))
    return {
        "task_id": TASK_ID,
        "deterministic_score": score,
        "criteria": criteria,
        "hardgate_pass": score == 80,
        "failure_codes": failures,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    args = parser.parse_args()
    result = grade(args.workspace.resolve())
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if result["hardgate_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())

