from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import statistics
from functools import lru_cache
from pathlib import Path

ACCEPTED = True

INPUT = Path(__file__).resolve().parents[2] / "inputs" / "ls07-combination-treatment-deg"
COUNT_SHA256 = "c233e3e1c68efea27818d97cc45dd1cceeb0bbb9e8693951b6964cd1a8f6ac17"
LAYOUT_SHA256 = "cb0212290186c6c5ecbd4d06f18e7fd0ad93d101ab0f07e7171e566e5cfa2584"
CONTROL = ("3_1", "3_2", "3_3")
TREATMENT = ("9_1", "9_2", "9_3")
EXPECTED_PASSING = 677

ALIASES = {
    "gene_id": ("gene_id", "ensembl_id", "ensg"),
    "gene_name": ("gene_name", "symbol"),
    "baseMean": ("baseMean", "base_mean", "mean_normalized_count"),
    "log2FoldChange": ("log2FoldChange", "log2_fold_change", "log2fc"),
    "pvalue": ("pvalue", "p_value"),
    "padj": ("padj", "adjusted_p_value", "fdr"),
    "pass": ("pass", "significant", "is_significant"),
}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _number(value):
    if value is None or str(value).strip().lower() in {"", "na", "nan", "null", "none"}:
        return math.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _boolean(value):
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "y"}:
        return True
    if text in {"false", "0", "no", "n", ""}:
        return False
    return None


def _pick(row, key):
    for name in ALIASES[key]:
        if name in row:
            return row[name]
    return None


@lru_cache(maxsize=1)
def _truth():
    counts_path, layout_path = INPUT / "counts_raw_unfiltered.csv", INPUT / "sample_layout.csv"
    integrity = (_sha256(counts_path) == COUNT_SHA256 and _sha256(layout_path) == LAYOUT_SHA256)
    raw = {}
    with counts_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        id_col = reader.fieldnames[0]
        for row in reader:
            raw[row[id_col]] = tuple(int(row[s]) for s in CONTROL + TREATMENT)
    usable = [v for v in raw.values() if all(x > 0 for x in v)]
    geometric = [math.exp(sum(math.log(x) for x in v) / 6.0) for v in usable]
    size_factors = tuple(statistics.median(v[j] / g for v, g in zip(usable, geometric)) for j in range(6))
    truth = {}
    for gene, values in raw.items():
        normalized = tuple(v / sf for v, sf in zip(values, size_factors))
        base = sum(normalized) / 6.0
        ctrl, trt = sum(normalized[:3]) / 3.0, sum(normalized[3:]) / 3.0
        sign = 1 if trt > ctrl else (-1 if trt < ctrl else 0)
        truth[gene] = (base, sign)
    return integrity, truth


def _read_table(path: Path):
    try:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))
    except Exception:
        return []


def _summary_count(data):
    for key in ("n_passing", "passing_genes", "significant_gene_count", "count"):
        try:
            return int(data[key])
        except (KeyError, TypeError, ValueError):
            pass
    return -1


def check(workspace: Path):
    out = Path(workspace) / "output"
    rows = _read_table(out / "differential_expression.csv")
    try:
        summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    except Exception:
        summary = {}
    try:
        report = (out / "report.md").read_text(encoding="utf-8", errors="replace")
    except Exception:
        report = ""

    try:
        integrity, truth = _truth()
    except Exception:
        integrity, truth = False, {}

    ids = [str(_pick(r, "gene_id") or "").strip() for r in rows]
    ids_valid = bool(rows) and all(ids) and len(ids) == len(set(ids)) and all(g in truth for g in ids)
    base_checked = base_matched = 0
    direction_checked = direction_matched = 0
    predicate_checked = predicate_matched = passing = 0
    missing_preserved = True
    for gene, row in zip(ids, rows):
        bm, lfc = _number(_pick(row, "baseMean")), _number(_pick(row, "log2FoldChange"))
        pval, padj = _number(_pick(row, "pvalue")), _number(_pick(row, "padj"))
        reported = _boolean(_pick(row, "pass"))
        if gene in truth and math.isfinite(bm):
            base_checked += 1
            expected_base, observed_sign = truth[gene]
            base_matched += int(math.isclose(bm, expected_base, rel_tol=2e-4, abs_tol=2e-3))
            if math.isfinite(lfc) and abs(lfc) > 0.5 and expected_base > 10 and observed_sign:
                direction_checked += 1
                direction_matched += int((lfc > 0) == (observed_sign > 0))
        expected_pass = (math.isfinite(bm) and math.isfinite(lfc) and math.isfinite(padj)
                         and padj < 0.05 and abs(lfc) > 0.5 and bm > 10)
        if reported is not None:
            predicate_checked += 1
            predicate_matched += int(reported == expected_pass)
            passing += int(reported)
        if not math.isfinite(padj) and reported is True:
            missing_preserved = False
        if str(_pick(row, "padj") or "").strip().lower() in {"nan", "inf", "-inf"}:
            missing_preserved = False
        if math.isfinite(padj) and (padj < 0 or padj > 1):
            missing_preserved = False
        if math.isfinite(pval) and (pval < 0 or pval > 1):
            missing_preserved = False

    base_rate = base_matched / base_checked if base_checked else 0.0
    direction_rate = direction_matched / direction_checked if direction_checked else 0.0
    predicate_rate = predicate_matched / predicate_checked if predicate_checked else 0.0
    summary_count = _summary_count(summary)
    text = (report + "\n" + json.dumps(summary, ensure_ascii=False)).lower()
    numerator = ("cisplatin_ic50_cbd_ic50" in text or
                 ("combination" in text and "cisplatin" in text and "cbd" in text))
    comparison = bool(re.search(r"(?:vs\.?|versus|compared\s+to|relative\s+to|reference|denominator|control)[^\n]{0,45}dmso|dmso[^\n]{0,45}(?:control|reference|denominator)", text))
    negated = bool(re.search(r"(?:not|doesn.t|didn.t|isn.t)[^\n]{0,45}(?:dmso|control|reference|contrast)|dmso[^\n]{0,35}(?:not|isn.t)[^\n]{0,20}(?:control|reference|denominator)", text))
    contrast_declared = numerator and comparison and not negated

    # Exactly three explicitly named fatal scientific gates.
    fatal_gates = {
        "INPUT_LINKED_QUANTITATIVE_GROUNDING": integrity and ids_valid and base_checked >= 100 and base_rate >= 0.95,
        "FROZEN_CONTRAST_DIRECTION": contrast_declared and direction_checked >= 50 and direction_rate >= 0.90,
        "OFFICIAL_ENDPOINT_AND_PREDICATE": (predicate_checked == len(rows) and predicate_rate >= 0.995
                                             and passing == EXPECTED_PASSING
                                             and summary_count == EXPECTED_PASSING),
    }

    core = 0
    core += 8 if integrity and ids_valid else 0
    core += round(12 * min(1.0, base_rate / 0.95)) if base_checked >= 100 else 0
    core += round(10 * predicate_rate) if predicate_checked == len(rows) and rows else 0
    core += 10 if passing == EXPECTED_PASSING else 0
    direction = 0
    if direction_checked >= 50:
        direction += round(10 * min(1.0, direction_rate / 0.90))
    direction += 5 if contrast_declared and direction_checked >= 50 and direction_rate >= 0.90 else 0
    summary_score = 5 if (summary_count == passing == EXPECTED_PASSING and predicate_rate >= 0.995
                          and missing_preserved) else 0

    failures = []
    if not fatal_gates["INPUT_LINKED_QUANTITATIVE_GROUNDING"]:
        failures.append("FATAL_INPUT_LINKED_QUANTITATIVE_GROUNDING")
    if not fatal_gates["FROZEN_CONTRAST_DIRECTION"]:
        failures.append("FATAL_FROZEN_CONTRAST_DIRECTION")
    if not fatal_gates["OFFICIAL_ENDPOINT_AND_PREDICATE"]:
        failures.append("FATAL_OFFICIAL_ENDPOINT_AND_PREDICATE")
    if not missing_preserved:
        failures.append("INVALID_MISSING_OR_PROBABILITY_VALUES")

    criteria = {
        "fatal_gates": fatal_gates,
        "authorized_input_integrity": integrity,
        "grounded_unique_gene_ids": ids_valid,
        "base_mean_anchor_rows": base_checked,
        "base_mean_match_rate": round(base_rate, 6),
        "direction_informative_rows": direction_checked,
        "direction_sign_match_rate": round(direction_rate, 6),
        "contrast_declared_without_negation": contrast_declared,
        "pass_predicate_match_rate": round(predicate_rate, 6),
        "passing_rows": passing,
        "official_passing_count": EXPECTED_PASSING,
        "summary_count": summary_count,
        "missing_values_preserved": missing_preserved,
    }
    return {"core_science": max(0, min(40, core)), "direction": max(0, min(15, direction)),
            "summary": max(0, min(5, summary_score)), "hardgate_pass": all(fatal_gates.values()),
            "criteria": criteria, "failure_codes": failures}
