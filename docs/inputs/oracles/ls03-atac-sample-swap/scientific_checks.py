from __future__ import annotations

import csv
import gzip
import json
import math
import re
from functools import lru_cache
from pathlib import Path

import numpy as np

ACCEPTED = True

INPUT_MATRIX = (
    Path(__file__).resolve().parents[2]
    / "inputs"
    / "ls03-atac-sample-swap"
    / "sample.swap.atac.q1.tsv.gz"
)


def _norm(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def _number(value: object) -> float | None:
    try:
        ans = float(str(value).strip())
        return ans if math.isfinite(ans) else None
    except (TypeError, ValueError):
        return None


def _flatten(value: object) -> str:
    if isinstance(value, dict):
        return " ".join(_flatten(v) for v in value.values())
    if isinstance(value, (list, tuple)):
        return " ".join(_flatten(v) for v in value)
    return "" if value is None else str(value)


@lru_cache(maxsize=1)
def _reference() -> tuple[list[str], dict[str, np.ndarray]]:
    """Recompute three standard whole-matrix similarities from the authorized input."""
    if not INPUT_MATRIX.is_file():
        raise FileNotFoundError(INPUT_MATRIX)
    with gzip.open(INPUT_MATRIX, "rt", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    organs = header[3:]
    n = 0
    sums = sums2 = cross = raw_sums = raw2 = raw_cross = None
    with gzip.open(INPUT_MATRIX, "rt", encoding="utf-8", newline="") as handle:
        next(handle)
        chunk: list[np.ndarray] = []
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")[3:]
            if len(fields) != len(organs):
                continue
            chunk.append(np.asarray(fields, dtype=np.float64))
            if len(chunk) < 50000:
                continue
            n, sums, sums2, cross, raw_sums, raw2, raw_cross = _accumulate(
                np.vstack(chunk), n, sums, sums2, cross, raw_sums, raw2, raw_cross
            )
            chunk = []
        if chunk:
            n, sums, sums2, cross, raw_sums, raw2, raw_cross = _accumulate(
                np.vstack(chunk), n, sums, sums2, cross, raw_sums, raw2, raw_cross
            )
    if n == 0 or sums is None:
        raise ValueError("authorized ATAC matrix contains no usable rows")
    log_var = sums2 - np.outer(sums, sums) / n
    raw_var = raw2 - np.outer(raw_sums, raw_sums) / n
    log_corr = (cross - np.outer(sums, sums) / n) / np.sqrt(
        np.outer(np.diag(log_var), np.diag(log_var))
    )
    raw_corr = (raw_cross - np.outer(raw_sums, raw_sums) / n) / np.sqrt(
        np.outer(np.diag(raw_var), np.diag(raw_var))
    )
    cosine = raw_cross / np.sqrt(np.outer(np.diag(raw2), np.diag(raw2)))
    return organs, {"log1p_pearson": log_corr, "raw_pearson": raw_corr, "cosine": cosine}


def _accumulate(x, n, sums, sums2, cross, raw_sums, raw2, raw_cross):
    y = np.log1p(x)
    ys, ys2, yc = y.sum(0), y.T @ y, y.T @ y
    # ys2 is a matrix so its diagonal supplies per-column sums of squares.
    rs, r2, rc = x.sum(0), x.T @ x, x.T @ x
    if sums is None:
        return len(x), ys, ys2, yc, rs, r2, rc
    return n + len(x), sums + ys, sums2 + ys2, cross + yc, raw_sums + rs, raw2 + r2, raw_cross + rc


def _read_pairs(path: Path, organs: list[str]) -> list[tuple[str, str, float, bool]]:
    """Accept long/nearest-neighbour or square similarity layouts; never execute code."""
    if not path.is_file():
        return []
    try:
        with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
            rows = list(csv.reader(handle))
    except Exception:
        return []
    if len(rows) < 2:
        return []
    known = {_norm(o): o for o in organs}
    header = rows[0]
    organ_cols = [(j, known[_norm(h)]) for j, h in enumerate(header) if _norm(h) in known]
    pairs: list[tuple[str, str, float, bool]] = []
    # Square or sparse-square matrix: first column contains row organ labels.
    if len(organ_cols) >= 2:
        for row in rows[1:]:
            if not row:
                continue
            left = known.get(_norm(row[0]))
            if left is None:
                continue
            for j, right in organ_cols:
                value = _number(row[j]) if j < len(row) else None
                if value is not None and left != right:
                    pairs.append((left, right, value, False))
        return pairs
    names = {_norm(h): i for i, h in enumerate(header)}
    left_aliases = ("organa", "samplea", "tissuea", "organ1", "sample1", "tissue1", "row", "source")
    right_aliases = ("organb", "sampleb", "tissueb", "organ2", "sample2", "tissue2", "nearestorgan", "nearesttissue", "target")
    value_aliases = ("similarity", "correlation", "pearson", "pearsonr", "cosine", "score", "value", "distance", "dissimilarity")
    li = next((names[x] for x in left_aliases if x in names), None)
    ri = next((names[x] for x in right_aliases if x in names), None)
    vk = next((x for x in value_aliases if x in names), None)
    if li is None or ri is None or vk is None:
        return []
    vi, is_distance = names[vk], vk in {"distance", "dissimilarity"}
    for row in rows[1:]:
        if max(li, ri, vi) >= len(row):
            continue
        left, right = known.get(_norm(row[li])), known.get(_norm(row[ri]))
        value = _number(row[vi])
        if left and right and left != right and value is not None:
            pairs.append((left, right, value, is_distance))
    return pairs


def _similarity_evidence(path: Path):
    try:
        organs, methods = _reference()
    except Exception:
        return False, False, None, None, 0, 0, False
    pairs = _read_pairs(path, organs)
    index = {o: i for i, o in enumerate(organs)}
    unique = {tuple(sorted((a, b))) for a, b, _, _ in pairs}
    represented = {x for pair in unique for x in pair}
    breadth = len(unique) >= 8 and len(represented) >= 8
    best_name, best_mae, best_fraction = None, math.inf, 0.0
    for name, matrix in methods.items():
        errors = []
        for a, b, value, is_distance in pairs:
            expected = float(matrix[index[a], index[b]])
            expected = 1.0 - expected if is_distance else expected
            scaled = value / 100.0 if abs(value) > 1.0 and abs(value) <= 100.0 else value
            errors.append(abs(scaled - expected))
        if errors:
            mae = sum(errors) / len(errors)
            fraction = sum(e <= 0.05 for e in errors) / len(errors)
            if mae < best_mae:
                best_name, best_mae, best_fraction = name, mae, fraction
    numeric = bool(pairs) and best_mae <= 0.025 and best_fraction >= 0.85
    return numeric, breadth, best_name, (round(best_mae, 6) if math.isfinite(best_mae) else None), len(unique), len(represented), True


def check(workspace: Path):
    output = Path(workspace) / "output"
    try:
        call = json.loads((output / "swap_call.json").read_text(encoding="utf-8"))
        call_valid = isinstance(call, dict)
    except Exception:
        call, call_valid = {}, False
    unset = {"", "none", "null", "na", "n/a", "unknown", "undetermined", "notapplicable"}
    pair_unset = all(_norm(call.get(k)) in {_norm(v) for v in unset} for k in ("organ_a", "organ_b"))
    negative_call = call_valid and call.get("swap_detected") is False
    direction_ok = negative_call and pair_unset

    numeric, breadth, method, mae, pair_count, organ_count, input_ok = _similarity_evidence(
        output / "sample_similarity.csv"
    )
    grounded_similarity = input_ok and numeric and breadth
    report_path = output / "report.md"
    report = report_path.read_text(encoding="utf-8", errors="replace") if report_path.is_file() else ""
    text = (_flatten(call.get("evidence")) + " " + report).lower()
    ambiguity = bool(re.search(
        r"no\s+(?:unique|definitive|identifiable)|not\s+uniquely|cannot\s+(?:determine|identify|infer|assign)|"
        r"insufficient|ambig|uncertain|does\s+not\s+(?:prove|establish|identify)|"
        r"(?:not|no evidence to)\s+support(?:ed|ing)?\s+(?:a\s+|any\s+)?(?:unique\s+)?swap", text
    ))
    limitation = bool(re.search(
        r"replicat|reference|ground.?truth|verified\s+atlas|external|single\s+(?:sample|profile)|"
        r"one\s+(?:sample|profile).{0,25}(?:organ|tissue)|(?:similarity|correlation).{0,20}alone|"
        r"label\s+identity|identity\s+(?:reference|metadata)", text
    ))
    rationale_ok = ambiguity and limitation
    report_ok = bool(report.strip()) and ambiguity and bool(re.search(
        r"no\s+(?:unique|definitive)|swap\s+(?:was\s+)?not\s+detected|"
        r"cannot\s+(?:determine|identify|assign)|do(?:es)?\s+not\s+support.{0,25}swap|"
        r"not\s+enough\s+evidence.{0,25}swap", report.lower()
    ))

    core = (22 if numeric and input_ok else 0) + (10 if breadth else 0) + (8 if rationale_ok else 0)
    failures = []
    if not call_valid:
        failures.append("MISSING_OR_INVALID_SWAP_CALL")
    if not direction_ok:
        failures.append("UNSUPPORTED_SWAP_DIRECTION")
    if not breadth:
        failures.append("SIMILARITY_TABLE_INSUFFICIENT")
    if not input_ok:
        failures.append("INPUT_GROUNDING_FAILED")
    elif not numeric:
        failures.append("SIMILARITY_NOT_GROUNDED")
    if not rationale_ok:
        failures.append("UNCERTAINTY_NOT_JUSTIFIED")
    if not report_ok:
        failures.append("REPORT_MISMATCH")
    criteria = {
        "authorized_input_recomputed": input_ok,
        "negative_nonunique_call": negative_call,
        "organ_pair_unset": pair_unset,
        "similarity_numeric_agreement": numeric,
        "similarity_breadth": breadth,
        "matched_method": method,
        "mean_absolute_error": mae,
        "unique_pairs": pair_count,
        "organs_represented": organ_count,
        "uncertainty_identifiability_rationale": rationale_ok,
        "report_consistent": report_ok,
    }
    return {
        "core_science": core,
        "direction": 15 if direction_ok else 0,
        "summary": 5 if report_ok else 0,
        "hardgate_pass": direction_ok and grounded_similarity and rationale_ok,
        "criteria": criteria,
        "failure_codes": failures,
    }
