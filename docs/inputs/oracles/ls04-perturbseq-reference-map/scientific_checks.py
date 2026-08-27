from __future__ import annotations

import ast
import csv
import math
import re
from pathlib import Path

ACCEPTED = True

# Independently recomputed from the two request-authorized H5AD files. Cells were
# library-size normalized to 10,000, log1p transformed, pseudobulked by guide,
# centered on pooled NT guides within each cell type, aligned on 11,858 common
# non-target-leaking genes, and compared by cosine similarity. A rectangular
# Hungarian assignment reproduces every expected pair below.
# query H5AD SHA256: 548b33e7015d884e33dce254fcb7707ac187ec5cb33ebbd45954f2623831327e
# ref H5AD SHA256:   42b0b42978321217db1d33a9de1b93b2662e37dd4877032defaf55388b377f20
TRUTH = {
    "PABPC1": {"guide": "guide18", "score": 0.5216191513, "runner": 0.2819277},
    "NUDT21": {"guide": "guide13", "score": 0.6313580675, "runner": 0.4984148},
    "LEO1": {"guide": "guide14", "score": 0.3784984638, "runner": 0.3137100},
}
FATAL_GATES = ("FATAL_MAPPING_IDENTITY", "FATAL_NO_QUERY_TARGET_LEAKAGE")
NEGATION = re.compile(
    r"\b(?:no|not|never|neither|without|false|wrong|incorrect|unsupported|"
    r"doesn['’]?t|does\s+not|isn['’]?t|is\s+not|cannot|can['’]?t|fails?\s+to)\b",
    re.I,
)


def _norm(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def _number(value: object) -> float | None:
    try:
        ans = float(str(value).strip())
        return ans if math.isfinite(ans) else None
    except (TypeError, ValueError):
        return None


def _rows(workspace: Path) -> tuple[dict[str, dict[str, str]], list[str]]:
    path = workspace / "output" / "guide_mapping.csv"
    if not path.is_file():
        return {}, []
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            keys = {_norm(name): name for name in fields}
            required = ["targetgene", "queryguideid", "score", "runnerupscore", "confidence"]
            if any(key not in keys for key in required):
                return {}, fields
            grouped: dict[str, list[dict[str, str]]] = {}
            for raw in reader:
                target = str(raw.get(keys["targetgene"], "")).strip().upper()
                if target in TRUTH:
                    grouped.setdefault(target, []).append(
                        {key: str(raw.get(keys[key], "")).strip() for key in required}
                    )
            # Multiple rows for one requested target are scientifically ambiguous,
            # even when a candidate repeats the same self-reported mapping.
            return {target: values[0] for target, values in grouped.items() if len(values) == 1}, fields
    except Exception:
        return {}, []


def _chain(node: ast.AST) -> list[str]:
    if isinstance(node, ast.Name):
        return [node.id]
    if isinstance(node, ast.Attribute):
        return _chain(node.value) + [node.attr]
    if isinstance(node, ast.Subscript):
        out = _chain(node.value)
        key = node.slice
        if isinstance(key, ast.Constant) and isinstance(key.value, str):
            out.append(key.value)
        return out
    return []


def _query_target_leakage(source: str) -> bool:
    """Detect executable reads of target_gene from query observations; never execute source."""
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return False  # Syntax/executability belongs to the fixed script checker.

    query_vars = {"q", "query", "query_adata", "adata_q", "q_adata", "query_obs", "q_obs"}
    changed = True
    while changed:
        changed = False
        for node in ast.walk(tree):
            if not isinstance(node, (ast.Assign, ast.AnnAssign)):
                continue
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            value = node.value
            value_text = ast.unparse(value).lower() if value is not None else ""
            rhs_chain = _chain(value) if value is not None else []
            from_query = "query" in value_text or (rhs_chain and rhs_chain[0] in query_vars)
            if from_query:
                for target in targets:
                    if isinstance(target, ast.Name) and target.id not in query_vars:
                        query_vars.add(target.id)
                        changed = True

    query_names = {name.lower() for name in query_vars}
    for node in ast.walk(tree):
        chain = _chain(node)
        lowered = [part.lower() for part in chain]
        if lowered and lowered[0] in query_names and "obs" in lowered and "target_gene" in lowered:
            return True
        if isinstance(node, ast.Call):
            owner_low = [part.lower() for part in _chain(node.func)]
            args = [arg.value.lower() for arg in node.args if isinstance(arg, ast.Constant) and isinstance(arg.value, str)]
            if owner_low and owner_low[0] in query_names and "obs" in owner_low and "target_gene" in args:
                return True
    return False


def _segments(text: str) -> list[str]:
    return [part.strip() for part in re.split(r"(?:[\r\n]+|(?<=[.!?;])\s+)", text) if part.strip()]


def _affirmed_pair(text: str, left: str, right: str) -> bool:
    """Require positive evidence and reject any explicit contradiction for the pair."""
    relevant = [s for s in _segments(text) if re.search(rf"\b{re.escape(left)}\b", s, re.I) and re.search(rf"\b{re.escape(right)}\b", s, re.I)]
    return bool(relevant) and any(not NEGATION.search(s) for s in relevant) and not any(NEGATION.search(s) for s in relevant)


def _affirmed_leo1_ambiguity(text: str) -> bool:
    terms = re.compile(r"\b(?:ambig(?:uous|uity)?|weakest|smallest|lowest|least|runner[- ]?up|margin|confidence)\b", re.I)
    relevant = [s for s in _segments(text) if re.search(r"\bLEO1\b", s, re.I) and terms.search(s)]
    return bool(relevant) and any(not NEGATION.search(s) for s in relevant) and not any(NEGATION.search(s) for s in relevant)


def _confidence_rank(value: str) -> float | None:
    label = value.strip().lower()
    named = {"low": 0.0, "lowest": 0.0, "weak": 0.0, "medium": 0.5, "moderate": 0.5,
             "high": 1.0, "strong": 1.0}
    if label in named:
        return named[label]
    return _number(label)


def check(workspace: Path):
    workspace = Path(workspace)
    rows, _fields = _rows(workspace)
    mapping_ok: dict[str, bool] = {}
    numeric_ok: dict[str, bool] = {}
    margins: dict[str, float] = {}

    for target, expected in TRUTH.items():
        row = rows.get(target)
        mapping_ok[target] = bool(row and row["queryguideid"].lower() == expected["guide"].lower())
        score = _number(row.get("score")) if row else None
        runner = _number(row.get("runnerupscore")) if row else None
        numeric_ok[target] = bool(
            mapping_ok[target] and score is not None and runner is not None and score > runner
            and abs(score - expected["score"]) <= 0.15
            and abs(runner - expected["runner"]) <= 0.15
        )
        if score is not None and runner is not None and score > runner:
            margins[target] = (score - runner) / max(abs(score), 1e-12)

    mapping_gate = all(mapping_ok.values()) and len(rows) == len(TRUTH)
    analysis_path = workspace / "output" / "analysis.py"
    try:
        analysis_source = analysis_path.read_text(encoding="utf-8", errors="replace") if analysis_path.is_file() else ""
    except Exception:
        analysis_source = ""
    leakage = _query_target_leakage(analysis_source)
    leakage_gate = not leakage

    core = sum(14 if target in {"PABPC1", "NUDT21"} else 12 for target, ok in mapping_ok.items() if ok)
    numeric_points = 2 * sum(numeric_ok.values())
    ambiguity_ok = (
        all(numeric_ok.values()) and len(margins) == 3
        and margins["LEO1"] < margins["NUDT21"] < margins["PABPC1"]
    )
    ranks = {target: _confidence_rank(rows[target]["confidence"]) if target in rows else None for target in TRUTH}
    confidence_ok = (
        mapping_gate and all(value is not None for value in ranks.values())
        and ranks["PABPC1"] >= ranks["NUDT21"] >= ranks["LEO1"]
        and ranks["PABPC1"] > ranks["LEO1"]
    )
    leakage_points = 2 if mapping_gate and bool(analysis_source.strip()) and leakage_gate else 0
    direction = numeric_points + (5 if ambiguity_ok else 0) + (2 if confidence_ok else 0) + leakage_points

    report_path = workspace / "output" / "report.md"
    try:
        report = report_path.read_text(encoding="utf-8", errors="replace") if report_path.is_file() else ""
    except Exception:
        report = ""
    report_mapping = {
        target: mapping_ok[target] and _affirmed_pair(report, target, expected["guide"])
        for target, expected in TRUTH.items()
    }
    report_ambiguity = ambiguity_ok and _affirmed_leo1_ambiguity(report)
    summary = sum(report_mapping.values()) + (2 if report_ambiguity else 0)

    criteria = {
        "fatal_gates": {
            "FATAL_MAPPING_IDENTITY": mapping_gate,
            "FATAL_NO_QUERY_TARGET_LEAKAGE": leakage_gate,
        },
        "mapping_identity": mapping_ok,
        "numeric_scores_grounded": numeric_ok,
        "relative_margins": {key: round(value, 8) for key, value in margins.items()},
        "leo1_is_most_ambiguous": ambiguity_ok,
        "confidence_order_grounded": confidence_ok,
        "query_target_metadata_leakage_detected": leakage,
        "report_mapping_affirmed": report_mapping,
        "report_ambiguity_affirmed": report_ambiguity,
    }
    failures: list[str] = []
    if not mapping_gate:
        failures.append("FATAL_MAPPING_IDENTITY")
    if not leakage_gate:
        failures.append("FATAL_NO_QUERY_TARGET_LEAKAGE")
    if not all(numeric_ok.values()):
        failures.append("NUMERIC_ALIGNMENT_EVIDENCE_UNGROUNDED")
    if not ambiguity_ok:
        failures.append("AMBIGUITY_QUANTIFICATION_MISMATCH")
    if not confidence_ok:
        failures.append("CONFIDENCE_ORDER_MISMATCH")
    if summary < 5:
        failures.append("REPORT_SCIENTIFIC_SUMMARY_UNGROUNDED")

    return {
        "core_science": int(max(0, min(40, core))),
        "direction": int(max(0, min(15, direction))),
        "summary": int(max(0, min(5, summary))),
        "hardgate_pass": bool(mapping_gate and leakage_gate),
        "criteria": criteria,
        "failure_codes": failures,
    }
