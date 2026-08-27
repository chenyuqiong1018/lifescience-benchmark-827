from __future__ import annotations

import csv
import json
import math
import re
from decimal import Decimal, InvalidOperation
from pathlib import Path

ACCEPTED = True

# Immutable truth recomputed from the authorized MatrixMarket inputs.  The raw
# matrices contain 6,295 and 5,004 cells.  A conservative horizontal-cell
# annotation gives about 179 and 4 cells; marker-only sensitivity analyses
# (LHX1 and ONECUT3) independently show the same sample-2 collapse.
RAW_TOTALS = {1: 6295, 2: 5004}
HORIZONTAL_ALIASES = ("horizontal",)


def _sample_index(value):
    text = str(value).strip().casefold()
    compact = re.sub(r"[^a-z0-9]+", "", text)
    if text in {"1", "sample 1", "sample_1", "sample-1", "q1.1"} or compact in {"s1", "sample1", "q11"}:
        return 1
    if text in {"2", "sample 2", "sample_2", "sample-2", "q1.2"} or compact in {"s2", "sample2", "q12"}:
        return 2
    return None


def _horizontal(value):
    compact = re.sub(r"[^a-z0-9]+", "", str(value).casefold())
    return any(alias in compact for alias in HORIZONTAL_ALIASES) or compact in {"hc", "hcs"}


def _integer(value):
    try:
        number = Decimal(str(value).strip())
        if not number.is_finite() or number != number.to_integral_value() or number < 0:
            return None
        return int(number)
    except (InvalidOperation, ValueError):
        return None


def _fraction(value):
    raw = str(value).strip()
    try:
        dec = Decimal(raw)
        number = float(dec)
    except (InvalidOperation, ValueError, OverflowError):
        return None, None
    if not math.isfinite(number) or not 0 <= number <= 1:
        return None, None
    if number in {0.0, 1.0} and "." not in raw and "e" not in raw.casefold():
        tolerance = 1e-9
    else:
        exponent = dec.as_tuple().exponent
        tolerance = max(1e-9, 0.500001 * (10.0 ** exponent)) if exponent < 0 else 1e-9
    return number, tolerance


def _column(fieldnames, aliases):
    normalized = {re.sub(r"[^a-z0-9]+", "", name.casefold()): name for name in (fieldnames or [])}
    for alias in aliases:
        hit = normalized.get(re.sub(r"[^a-z0-9]+", "", alias.casefold()))
        if hit:
            return hit
    return None


def _read_composition(path):
    result = {"readable": False, "valid_rows": False, "rows": [], "totals": {}, "horizontal": {}}
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            sample_col = _column(reader.fieldnames, ("sample", "sample_id"))
            type_col = _column(reader.fieldnames, ("cell_type", "celltype", "population"))
            count_col = _column(reader.fieldnames, ("n_cells", "ncells", "count"))
            fraction_col = _column(reader.fieldnames, ("fraction", "proportion", "frac"))
            if not all((sample_col, type_col, count_col, fraction_col)):
                return result
            rows = []
            for raw in reader:
                sample = _sample_index(raw.get(sample_col, ""))
                count = _integer(raw.get(count_col, ""))
                fraction, tolerance = _fraction(raw.get(fraction_col, ""))
                cell_type = str(raw.get(type_col, "")).strip()
                if sample is None or count is None or fraction is None or not cell_type:
                    return result
                rows.append({"sample": sample, "cell_type": cell_type, "count": count,
                             "fraction": fraction, "tolerance": tolerance})
    except (OSError, csv.Error, UnicodeError):
        return result
    result["readable"] = True
    if not rows or {row["sample"] for row in rows} != {1, 2}:
        return result
    totals = {sample: sum(row["count"] for row in rows if row["sample"] == sample) for sample in (1, 2)}
    if not all(totals.values()):
        return result
    arithmetic = all(
        abs(row["fraction"] - row["count"] / totals[row["sample"]]) <= row["tolerance"] + 1e-12
        for row in rows
    )
    sums = all(
        abs(sum(row["fraction"] for row in rows if row["sample"] == sample) - 1.0)
        <= sum(row["tolerance"] for row in rows if row["sample"] == sample) + 1e-9
        for sample in (1, 2)
    )
    horizontal = {}
    for sample in (1, 2):
        selected = [row for row in rows if row["sample"] == sample and _horizontal(row["cell_type"])]
        horizontal[sample] = {
            "present": bool(selected),
            "count": sum(row["count"] for row in selected),
            "fraction": sum(row["count"] for row in selected) / totals[sample],
        }
    result.update(valid_rows=arithmetic and sums, rows=rows, totals=totals, horizontal=horizontal)
    return result


def _values_for_keys(obj, keys):
    values = []
    if isinstance(obj, dict):
        for key, value in obj.items():
            norm = re.sub(r"[^a-z0-9]+", "", str(key).casefold())
            if norm in keys and isinstance(value, (str, int, float)):
                values.append(value)
            values.extend(_values_for_keys(value, keys))
    elif isinstance(obj, list):
        for value in obj:
            values.extend(_values_for_keys(value, keys))
    return values


def _read_call(path):
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return {"readable": False, "horizontal": False, "sample2": False, "sample1": False}
    type_keys = {"celltype", "depletedcelltype", "population", "depletedpopulation", "cellpopulation"}
    sample_keys = {"sample", "depletedsample", "depletedin", "sampledepleted", "direction"}
    types = _values_for_keys(data, type_keys)
    samples = _values_for_keys(data, sample_keys)
    return {
        "readable": True,
        "horizontal": any(_horizontal(value) for value in types),
        "sample2": any(_sample_index(value) == 2 for value in samples),
        "sample1": any(_sample_index(value) == 1 for value in samples),
    }


def _affirmed_conclusion(text):
    for sentence in re.split(r"[.!?\n]+", text.casefold()):
        has_identity = "horizontal" in sentence or re.search(r"\bhc\b", sentence)
        has_sample2 = bool(re.search(r"sample[ _-]*2|second sample|q1\.2", sentence))
        has_depletion = bool(re.search(r"deplet|severe|loss|reduc|rarer|collapse", sentence))
        negated = bool(re.search(r"\b(?:not|no|without|fails? to)\b.{0,30}(?:deplet|reduc|loss|collapse)", sentence))
        if has_identity and has_sample2 and has_depletion and not negated:
            return True
    return False


def _report_checks(path):
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        text = ""
    lower = text.casefold()
    markers = {marker for marker in ("lhx1", "prox1", "onecut1", "onecut2", "onecut3", "calb1") if marker in lower}
    marker_context = bool(re.search(r"marker|express|annotat|identity|evidence", lower))
    marker_negated = bool(re.search(r"(?:no|without|lack(?:ing|s|ed)?)\s+(?:horizontal[- ]cell\s+)?marker", lower))
    uncertainty = bool(re.search(r"uncertain|limitation|sensitivity|ambig|low[- ]count|sampling|qc|quality control", lower))
    return {"conclusion": _affirmed_conclusion(text),
            "marker_evidence": len(markers) >= 2 and marker_context and not marker_negated,
            "uncertainty": uncertainty}


def check(workspace: Path):
    workspace = Path(workspace)
    output = workspace / "output"
    composition = _read_composition(output / "composition.csv")
    call = _read_call(output / "depleted_call.json")
    report = _report_checks(output / "report.md")

    totals = composition.get("totals", {})
    totals_grounded = (
        4400 <= totals.get(1, -1) <= RAW_TOTALS[1]
        and 3500 <= totals.get(2, -1) <= RAW_TOTALS[2]
    )
    h1 = composition.get("horizontal", {}).get(1, {})
    h2 = composition.get("horizontal", {}).get(2, {})
    sample1_supported = (
        h1.get("present", False) and 100 <= h1.get("count", -1) <= 260
        and 0.015 <= h1.get("fraction", -1) <= 0.045
    )
    sample2_supported = (
        h2.get("present", False) and 0 <= h2.get("count", -1) <= 15
        and 0 <= h2.get("fraction", -1) <= 0.003
    )
    identity_agreement = h1.get("present", False) and h2.get("present", False) and call["horizontal"]
    severe_ratio = (
        sample1_supported and sample2_supported and h1["fraction"] > 0
        and h2["fraction"] / h1["fraction"] <= 0.15
        and h1["fraction"] - h2["fraction"] >= 0.015
    )
    call_direction = call["sample2"] and not call["sample1"]

    gates = {
        "VALID_COMPOSITION_PROVENANCE": composition["readable"] and composition["valid_rows"] and totals_grounded,
        "SUPPORTED_HORIZONTAL_IDENTITY": sample1_supported and sample2_supported and identity_agreement,
        "SEVERE_SAMPLE2_DEPLETION": severe_ratio and call_direction,
    }
    core = ((12 if sample1_supported else 0) + (12 if sample2_supported else 0)
            + (6 if composition["valid_rows"] else 0) + (4 if totals_grounded else 0)
            + (6 if identity_agreement else 0))
    direction = (8 if severe_ratio else 0) + (5 if call_direction else 0) + (2 if severe_ratio and call_direction else 0)
    summary = (2 if report["conclusion"] else 0) + (2 if report["marker_evidence"] else 0) + (1 if report["uncertainty"] else 0)

    criteria = {
        "fatal_gates": gates,
        "composition_readable": composition["readable"],
        "count_fraction_arithmetic": composition["valid_rows"],
        "grounded_sample_totals": totals_grounded,
        "sample1_horizontal_supported": sample1_supported,
        "sample2_horizontal_supported": sample2_supported,
        "cross_file_horizontal_identity": identity_agreement,
        "severe_fraction_depletion": severe_ratio,
        "depleted_call_explicitly_sample2": call_direction,
        "report_affirmed_conclusion": report["conclusion"],
        "report_marker_evidence": report["marker_evidence"],
        "report_uncertainty": report["uncertainty"],
    }
    failure_codes = []
    if not gates["VALID_COMPOSITION_PROVENANCE"]:
        failure_codes.append("INVALID_COMPOSITION_PROVENANCE")
    if not gates["SUPPORTED_HORIZONTAL_IDENTITY"]:
        failure_codes.append("UNSUPPORTED_HORIZONTAL_IDENTITY")
    if not gates["SEVERE_SAMPLE2_DEPLETION"]:
        failure_codes.append("SAMPLE2_DEPLETION_NOT_ESTABLISHED")
    if summary < 5:
        failure_codes.append("REPORT_SCIENTIFIC_SUMMARY_INCOMPLETE")
    return {"core_science": core, "direction": direction, "summary": summary,
            "hardgate_pass": all(gates.values()), "criteria": criteria,
            "failure_codes": failure_codes}
