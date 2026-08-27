from __future__ import annotations

import csv
import hashlib
import math
import re
from collections import defaultdict
from functools import lru_cache
from pathlib import Path


ACCEPTED = True

INPUT_CSV = (
    Path(__file__).resolve().parents[2]
    / "inputs"
    / "ls03-genome-coordinates"
    / "single_cell_dynamics_question.csv"
)
EXPECTED_INPUT_SHA256 = "fc7a5e166dcdba18e05275a1cd5d338b476153bc035fb718a4a945e14ac97646"
CONTACT_THRESHOLD_NM = 260.0
FATAL_GATES = (
    "AUTHORIZED_TRUTH_AVAILABLE",
    "CORE_RESULTS_GROUNDED",
    "CAUSAL_INTERPRETATION_SAFE",
)


def _norm(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def _finite(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@lru_cache(maxsize=1)
def _truth() -> dict:
    if not INPUT_CSV.is_file():
        raise FileNotFoundError(f"authorized input missing: {INPUT_CSV}")
    observed_sha = _sha256(INPUT_CSV)
    if observed_sha != EXPECTED_INPUT_SHA256:
        raise ValueError(
            f"authorized input SHA-256 mismatch: {observed_sha}; "
            f"expected {EXPECTED_INPUT_SHA256}"
        )

    series: dict[int, list[tuple[int, float, int, float]]] = defaultdict(list)
    with INPUT_CSV.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "cell_id", "time", "enh_x", "enh_y", "enh_z",
            "prom_x", "prom_y", "prom_z", "transcription",
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("authorized input schema is incomplete")
        for row in reader:
            cell = int(row["cell_id"])
            time = int(float(row["time"]))
            distance = math.sqrt(
                (float(row["enh_x"]) - float(row["prom_x"])) ** 2
                + (float(row["enh_y"]) - float(row["prom_y"])) ** 2
                + (float(row["enh_z"]) - float(row["prom_z"])) ** 2
            )
            contact = int(distance <= CONTACT_THRESHOLD_NM)
            transcription = float(row["transcription"])
            series[cell].append((time, distance, contact, transcription))

    cell_metrics = {}
    for cell, values in series.items():
        values.sort(key=lambda item: item[0])
        n = len(values)
        cell_metrics[cell] = {
            "mean_distance_nm": sum(item[1] for item in values) / n,
            "contact_fraction": sum(item[2] for item in values) / n,
            "transcription_fraction": sum(item[3] for item in values) / n,
            "n_observations": n,
        }
    return {
        "sha256": observed_sha,
        "series": dict(series),
        "cell_metrics": cell_metrics,
        "n_rows": sum(len(values) for values in series.values()),
    }


def _pearson(x: list[float], y: list[float]) -> float:
    n = len(x)
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    sum_xx = sum((value - mean_x) ** 2 for value in x)
    sum_yy = sum((value - mean_y) ** 2 for value in y)
    if sum_xx <= 0 or sum_yy <= 0:
        return math.nan
    return sum(
        (left - mean_x) * (right - mean_y) for left, right in zip(x, y)
    ) / math.sqrt(sum_xx * sum_yy)


@lru_cache(maxsize=512)
def _lag_truth(lag: int) -> dict[str, float | int]:
    contacts: list[float] = []
    distances: list[float] = []
    outcomes: list[float] = []
    for values in _truth()["series"].values():
        transcription_at = {item[0]: item[3] for item in values}
        for time, distance, contact, _ in values:
            target = time + lag
            if target in transcription_at:
                contacts.append(float(contact))
                distances.append(distance)
                outcomes.append(transcription_at[target])
    contacted = [outcome for contact, outcome in zip(contacts, outcomes) if contact == 1]
    not_contacted = [
        outcome for contact, outcome in zip(contacts, outcomes) if contact == 0
    ]
    risk_difference = (
        sum(contacted) / len(contacted) - sum(not_contacted) / len(not_contacted)
        if contacted and not_contacted
        else math.nan
    )
    return {
        "n_observations": len(outcomes),
        "contact_pearson": _pearson(contacts, outcomes),
        "distance_pearson": _pearson(distances, outcomes),
        "proximity_pearson": _pearson([-value for value in distances], outcomes),
        "contact_risk_difference": risk_difference,
    }


def _read_csv(path: Path, max_rows: int = 1_000_000) -> tuple[list[dict], list[str]]:
    if not path.is_file() or path.stat().st_size > 64 * 1024 * 1024:
        return [], []
    try:
        with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            rows = []
            for index, row in enumerate(reader):
                if index >= max_rows:
                    return [], fields
                rows.append(row)
            return rows, fields
    except (OSError, csv.Error):
        return [], []


def _field(fields: list[str], aliases: set[str]) -> str | None:
    mapping = {_norm(name): name for name in fields}
    for alias in aliases:
        if alias in mapping:
            return mapping[alias]
    return None


def _fraction(value: object, header: str) -> float | None:
    number = _finite(value)
    if number is None:
        return None
    normalized = _norm(header)
    if "percent" in normalized or "pct" in normalized:
        return number / 100.0
    if 1 < number <= 100:
        return number / 100.0
    return number


def _parse_cell_metrics(path: Path) -> tuple[dict[int, dict[str, float]], dict]:
    rows, fields = _read_csv(path)
    cell_col = _field(fields, {"cellid", "cell", "id"})
    distance_col = _field(
        fields,
        {
            "meandistance", "meandistancenm", "mean3ddistance",
            "mean3ddistancenm", "avgdistance", "averagedistance",
            "averagedistancenm",
        },
    )
    contact_col = _field(
        fields,
        {
            "contactfraction", "fractioncontact", "contactrate",
            "contactproportion", "pctcontact", "contactpercent",
        },
    )
    transcription_col = _field(
        fields,
        {
            "transcriptionfraction", "fractiontranscription",
            "transcriptionrate", "transcriptionproportion",
            "pcttranscription", "transcriptionpercent", "activefraction",
        },
    )
    parsed: dict[int, dict[str, float]] = {}
    duplicates = 0
    if all((cell_col, distance_col, contact_col, transcription_col)):
        for row in rows:
            try:
                cell_float = float(row[cell_col])
                cell = int(cell_float)
            except (TypeError, ValueError):
                continue
            if cell_float != cell or cell in parsed:
                duplicates += 1
                continue
            distance = _finite(row.get(distance_col))
            contact = _fraction(row.get(contact_col), contact_col)
            transcription = _fraction(row.get(transcription_col), transcription_col)
            if None not in (distance, contact, transcription):
                parsed[cell] = {
                    "mean_distance_nm": distance,
                    "contact_fraction": contact,
                    "transcription_fraction": transcription,
                }
    return parsed, {
        "rows": len(rows),
        "duplicates_or_invalid_ids": duplicates,
        "recognized_columns": all(
            (cell_col, distance_col, contact_col, transcription_col)
        ),
    }


def _parse_lags(path: Path) -> tuple[list[dict[str, float | int]], dict]:
    rows, fields = _read_csv(path)
    lag_col = _field(fields, {"lag", "lagtime", "timelag", "lagsteps", "lagstep"})
    association_col = _field(
        fields,
        {
            "association", "correlation", "coefficient", "effect",
            "associationvalue", "corr", "r",
        },
    )
    n_col = _field(
        fields,
        {"nobservations", "nobservation", "nobs", "npairs", "n", "samplesize"},
    )
    parsed = []
    seen = set()
    duplicates = 0
    if all((lag_col, association_col, n_col)):
        for row in rows:
            lag_number = _finite(row.get(lag_col))
            association = _finite(row.get(association_col))
            n_number = _finite(row.get(n_col))
            if None in (lag_number, association, n_number):
                continue
            lag = int(lag_number)
            n_observations = int(n_number)
            if (
                lag_number != lag
                or abs(lag) >= 250
                or n_number != n_observations
                or lag in seen
            ):
                duplicates += 1
                continue
            seen.add(lag)
            parsed.append(
                {
                    "lag": lag,
                    "association": association,
                    "n_observations": n_observations,
                }
            )
    return parsed, {
        "rows": len(rows),
        "duplicates_or_invalid_lags": duplicates,
        "recognized_columns": all((lag_col, association_col, n_col)),
    }


def _cell_assessment(candidate: dict[int, dict[str, float]], expected: dict) -> dict:
    overlap = sorted(set(candidate).intersection(expected))
    total = len(expected)
    tolerances = {
        "mean_distance_nm": 0.05,
        "contact_fraction": 0.0011,
        "transcription_fraction": 0.0011,
    }
    correct = {}
    for metric, tolerance in tolerances.items():
        correct[metric] = sum(
            abs(candidate[cell][metric] - expected[cell][metric]) <= tolerance
            for cell in overlap
        )
    return {
        "expected_cells": total,
        "matched_cells": len(overlap),
        "coverage": len(overlap) / total if total else 0.0,
        "distance_accuracy": correct["mean_distance_nm"] / total if total else 0.0,
        "contact_accuracy": correct["contact_fraction"] / total if total else 0.0,
        "transcription_accuracy": (
            correct["transcription_fraction"] / total if total else 0.0
        ),
    }


def _lag_assessment(rows: list[dict[str, float | int]]) -> dict:
    methods = (
        "contact_pearson",
        "distance_pearson",
        "proximity_pearson",
        "contact_risk_difference",
    )
    fits = []
    for orientation in (1, -1):
        for method in methods:
            comparisons = []
            for row in rows:
                truth = _lag_truth(orientation * int(row["lag"]))
                expected = float(truth[method])
                error = abs(float(row["association"]) - expected)
                comparisons.append((row, expected, error, truth))
            rmse = (
                math.sqrt(sum(item[2] ** 2 for item in comparisons) / len(comparisons))
                if comparisons
                else math.inf
            )
            fits.append((rmse, orientation, method, comparisons))
    rmse, orientation, method, comparisons = min(fits, key=lambda item: item[0])
    tolerance = 0.0025
    association_correct = [
        item for item in comparisons if item[2] <= tolerance
    ]
    informative = [
        item for item in comparisons if abs(item[1]) >= 0.02
    ]
    informative_correct = [
        item for item in informative if item[2] <= tolerance
    ]
    n_correct = [
        item
        for item in comparisons
        if int(item[0]["n_observations"]) == int(item[3]["n_observations"])
    ]
    peak_correct = False
    if comparisons:
        candidate_peak = max(
            comparisons, key=lambda item: abs(float(item[0]["association"]))
        )
        expected_peak = max(comparisons, key=lambda item: abs(item[1]))
        peak_correct = (
            candidate_peak[0]["lag"] == expected_peak[0]["lag"]
            and candidate_peak[2] <= tolerance
        )
    association_fraction = (
        len(association_correct) / len(comparisons) if comparisons else 0.0
    )
    informative_fraction = (
        len(informative_correct) / len(informative) if informative else 0.0
    )
    n_fraction = len(n_correct) / len(comparisons) if comparisons else 0.0
    return {
        "rows_evaluated": len(comparisons),
        "best_method": method,
        "best_orientation": (
            "association_source_at_t_vs_transcription_at_t_plus_lag"
            if orientation == 1
            else "opposite_lag_sign_convention"
        ),
        "rmse": rmse if math.isfinite(rmse) else None,
        "association_accuracy": association_fraction,
        "informative_association_accuracy": informative_fraction,
        "n_observations_accuracy": n_fraction,
        "has_zero_lag": any(int(row["lag"]) == 0 for row in rows),
        "has_one_step_lag": any(abs(int(row["lag"])) == 1 for row in rows),
        "has_informative_truth": bool(informative),
        "peak_direction_correct": peak_correct,
    }


def _sentences(text: str) -> list[str]:
    return [
        piece.strip()
        for piece in re.split(r"(?<=[.!?;])\s+|\n+", text.lower())
        if piece.strip()
    ]


def _has_causal_disclaimer(text: str) -> bool:
    causal = re.compile(r"caus|cause|directionalit|cause.and.effect")
    limitation = re.compile(
        r"cannot|can't|can not|does not|doesn't|do not|don't|not possible|"
        r"not establish|not prove|not infer|not support|unable|insufficient|"
        r"no evidence|preclud|unknown|unresolved|cannot conclude|not permit"
    )
    inferential = re.compile(
        r"observational|association|correlation|infer|establish|prove|conclude|"
        r"determine|support"
    )
    for sentence in _sentences(text):
        if causal.search(sentence) and limitation.search(sentence):
            return True
        if (
            "observational" in sentence
            and limitation.search(sentence)
            and inferential.search(sentence)
        ):
            return True
    return False


def _unsafe_causal_claims(text: str) -> list[str]:
    biological_left = re.compile(
        r"contact|proximity|distance|enhancer.{0,30}promoter|promoter.{0,30}enhancer"
    )
    biological_right = re.compile(r"transcription|expression|rna")
    causal_claim = re.compile(
        r"\bcaus(?:e|es|ed|ing)\b|\bdriv(?:e|es|en|ing)\b|"
        r"\binduc(?:e|es|ed|ing)\b|\bdetermin(?:e|es|ed|ing)\b|"
        r"\bleads?\s+to\b|\bresults?\s+in\b|"
        r"\bcausal\s+(?:effect|relationship|mechanism)\b|"
        r"\b(?:prove|proves|proved|demonstrates?|establishes?)\b.{0,35}\bcaus"
    )
    negated_or_hedged = re.compile(
        r"cannot|can't|can not|does not|doesn't|do not|don't|did not|"
        r"not\s+(?:cause|caused|causal|establish|prove|support|shown)|"
        r"no evidence|insufficient|unknown|unresolved|whether|"
        r"\bmay\b|\bmight\b|\bcould\b|hypothes|consistent with|"
        r"would require|fails? to|question"
    )
    unsafe = []
    clauses = re.split(
        r"[.;\n]+|\b(?:but|however|although|though|yet|nevertheless)\b",
        text.lower(),
    )
    for clause in clauses:
        if (
            biological_left.search(clause)
            and biological_right.search(clause)
            and causal_claim.search(clause)
            and not negated_or_hedged.search(clause)
        ):
            unsafe.append(clause.strip()[:240])
    return unsafe


def _report_evidence(text: str, lag: dict) -> dict[str, bool]:
    lower = text.lower()
    numbers = []
    for match in re.finditer(r"(?<![a-z0-9])(-?\d+(?:\.\d+)?)\s*(%)?", lower):
        value = float(match.group(1))
        numbers.append(value / 100.0 if match.group(2) else value)

    def near(target: float, tolerance: float) -> bool:
        return any(abs(value - target) <= tolerance for value in numbers)

    return {
        "states_260_nm_threshold": bool(
            re.search(r"\b260(?:\.0+)?\s*(?:nm|nanomet)", lower)
        ),
        "states_grounded_mean_distance": near(526.0669983528302, 0.6),
        "states_grounded_contact_fraction": (
            near(0.05336, 0.001) or near(5.336, 0.15)
        ),
        "states_grounded_transcription_fraction": (
            near(0.15390666666666666, 0.0015) or near(15.3906666667, 0.2)
        ),
        "summarizes_temporal_association": bool(
            re.search(r"\b(?:lag|lead|tempor|time.step|one.step)\b", lower)
            and re.search(r"associat|correlat|contact|distance|proximity", lower)
            and lag.get("peak_direction_correct", False)
        ),
    }


def check(workspace: Path) -> dict:
    workspace = Path(workspace)
    failures: list[str] = []
    criteria: dict[str, object] = {}

    truth_available = True
    truth_error = ""
    try:
        truth = _truth()
    except Exception as exc:
        truth_available = False
        truth_error = str(exc)
        truth = {"cell_metrics": {}, "sha256": None}

    criteria["fatal_gate_AUTHORIZED_TRUTH_AVAILABLE"] = {
        "fatal": True,
        "passed": truth_available,
        "input_sha256": truth.get("sha256"),
        "error": truth_error or None,
    }

    cell_rows, cell_parse = _parse_cell_metrics(
        workspace / "output" / "cell_metrics.csv"
    )
    lag_rows, lag_parse = _parse_lags(workspace / "output" / "lag_analysis.csv")
    cell = _cell_assessment(cell_rows, truth["cell_metrics"])
    lag = _lag_assessment(lag_rows) if truth_available else {
        "rows_evaluated": 0,
        "best_method": None,
        "best_orientation": None,
        "rmse": None,
        "association_accuracy": 0.0,
        "informative_association_accuracy": 0.0,
        "n_observations_accuracy": 0.0,
        "has_zero_lag": False,
        "has_one_step_lag": False,
        "has_informative_truth": False,
        "peak_direction_correct": False,
    }
    criteria["cell_metrics_recomputed"] = {**cell_parse, **cell}
    criteria["lag_analysis_recomputed"] = {**lag_parse, **lag}

    cell_grounded = (
        cell["coverage"] >= 0.90
        and cell["distance_accuracy"] >= 0.90
        and cell["contact_accuracy"] >= 0.90
        and cell["transcription_accuracy"] >= 0.90
    )
    lag_grounded = (
        lag["rows_evaluated"] >= 3
        and lag["has_zero_lag"]
        and lag["has_one_step_lag"]
        and lag["has_informative_truth"]
        and lag["association_accuracy"] >= 0.80
        and lag["informative_association_accuracy"] >= 0.80
        and lag["n_observations_accuracy"] >= 0.90
        and lag["peak_direction_correct"]
        and lag["rmse"] is not None
        and lag["rmse"] <= 0.005
    )
    core_grounded = truth_available and cell_grounded and lag_grounded
    criteria["fatal_gate_CORE_RESULTS_GROUNDED"] = {
        "fatal": True,
        "passed": core_grounded,
        "cell_metrics_grounded": cell_grounded,
        "lag_analysis_grounded": lag_grounded,
    }

    report_path = workspace / "output" / "report.md"
    try:
        report = (
            report_path.read_text(encoding="utf-8", errors="replace")
            if report_path.is_file() and report_path.stat().st_size <= 2 * 1024 * 1024
            else ""
        )
    except OSError:
        report = ""
    disclaimer = _has_causal_disclaimer(report)
    unsafe_claims = _unsafe_causal_claims(report)
    causal_safe = disclaimer and not unsafe_claims
    criteria["fatal_gate_CAUSAL_INTERPRETATION_SAFE"] = {
        "fatal": True,
        "passed": causal_safe,
        "observational_causal_limitation_present": disclaimer,
        "unsupported_causal_claims": unsafe_claims,
    }

    distance_score = round(12 * cell["distance_accuracy"])
    contact_score = round(10 * cell["contact_accuracy"])
    transcription_score = round(8 * cell["transcription_accuracy"])
    association_ratio = (
        0.5 * lag["association_accuracy"]
        + 0.5 * lag["informative_association_accuracy"]
    )
    association_score = round(8 * association_ratio)
    n_score = round(2 * lag["n_observations_accuracy"])
    core_science = max(
        0,
        min(
            40,
            distance_score + contact_score + transcription_score
            + association_score + n_score,
        ),
    )
    criteria["core_science_scoring"] = {
        "distance_3d": distance_score,
        "contact_at_inclusive_260_nm": contact_score,
        "transcription_dynamics": transcription_score,
        "lagged_association": association_score,
        "lag_pair_counts": n_score,
    }

    direction = (
        (7 if lag["peak_direction_correct"] and lag["has_informative_truth"] else 0)
        + (
            3
            if re.search(r"\b(?:lag|lead|tempor|time.step|one.step)\b", report.lower())
            and re.search(r"associat|correlat", report.lower())
            else 0
        )
        + (5 if causal_safe else 0)
    )
    direction = max(0, min(15, direction))
    criteria["direction_scoring"] = {
        "data_supported_lag_direction": (
            lag["peak_direction_correct"] and lag["has_informative_truth"]
        ),
        "temporal_association_interpreted": bool(
            re.search(r"\b(?:lag|lead|tempor|time.step|one.step)\b", report.lower())
            and re.search(r"associat|correlat", report.lower())
        ),
        "association_separated_from_causation": causal_safe,
    }

    report_checks = _report_evidence(report, lag)
    summary = sum(report_checks.values())
    criteria["summary_scoring"] = report_checks

    gates = {
        "AUTHORIZED_TRUTH_AVAILABLE": truth_available,
        "CORE_RESULTS_GROUNDED": core_grounded,
        "CAUSAL_INTERPRETATION_SAFE": causal_safe,
    }
    for name, passed in gates.items():
        if not passed:
            failures.append(f"FATAL_{name}")
    if cell["distance_accuracy"] < 1:
        failures.append("CELL_DISTANCE_MISMATCH")
    if cell["contact_accuracy"] < 1:
        failures.append("CONTACT_THRESHOLD_METRICS_MISMATCH")
    if cell["transcription_accuracy"] < 1:
        failures.append("TRANSCRIPTION_METRICS_MISMATCH")
    if lag["association_accuracy"] < 1 or not lag["peak_direction_correct"]:
        failures.append("LAG_ASSOCIATION_MISMATCH")
    if lag["n_observations_accuracy"] < 1:
        failures.append("LAG_N_OBSERVATIONS_MISMATCH")
    if not all(report_checks.values()):
        failures.append("REPORT_SUMMARY_INCOMPLETE")

    return {
        "core_science": int(core_science),
        "direction": int(direction),
        "summary": int(summary),
        "hardgate_pass": all(gates.values()),
        "criteria": criteria,
        "failure_codes": failures,
    }
