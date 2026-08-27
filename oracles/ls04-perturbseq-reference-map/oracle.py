#!/usr/bin/env python3
"""Self-contained deterministic Oracle for ls04-perturbseq-reference-map."""
from __future__ import annotations
# Task-specific scientific scoring implementation.
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
# Shared deterministic artifact validation and command-line runner.
"""Static, fail-closed oracle runner for life-science task artifacts.

This module never imports or executes output/analysis.py. A task becomes scoreable
only after its local scientific_checks.py has passed the acceptance suite and
ACCEPTED is set to True there.
"""
import argparse
import ast
import csv
import importlib.util
import json
import math
from pathlib import Path
from typing import Any


REQUIRED_OUTPUTS = {
    "ls01-grna-offtarget-rank": ["ranked_guides.csv", "analysis.py", "report.md"],
    "ls01-primer-transcript-audit": ["primer_audit.csv", "analysis.py", "report.md"],
    "ls01-vector-orf-audit": ["construct_audit.csv", "analysis.py", "report.md"],
    "ls02-deleterious-mutation": ["variant.tsv", "evidence.json", "analysis.py", "report.md"],
    "ls02-find-deletion": ["deletion.tsv", "qc.json", "analysis.py", "report.md"],
    "ls02-infer-genome-build": ["build_call.json", "analysis.py", "report.md"],
    "ls03-cryptic-exon": ["cryptic_exon.tsv", "junctions.tsv", "analysis.py", "report.md"],
    "ls03-atac-sample-swap": ["swap_call.json", "sample_similarity.csv", "analysis.py", "report.md"],
    "ls03-genome-coordinates": ["cell_metrics.csv", "lag_analysis.csv", "analysis.py", "report.md"],
    "ls04-differential-composition": ["composition.csv", "depleted_call.json", "analysis.py", "report.md"],
    "ls04-perturbseq-reference-map": ["guide_mapping.csv", "analysis.py", "report.md"],
    "ls04-spatial-deconvolution": ["spot_710_composition.csv", "analysis.py", "report.md"],
    "ls05-protein-shape": ["shape_call.json", "shape_view.png"],
    "ls05-structure-model-ranking": ["model_ranking.csv", "analysis.py", "report.md"],
    "ls05-low-confidence-pocket": ["mutation_priorities.csv", "pocket_assessment.json", "analysis.py", "report.md"],
    "ls06-eno1-effect-size": ["eno1_effect.json", "analysis.py", "report.md"],
    "ls06-eno1-significance-audit": ["eno1_significance.json", "analysis.py", "report.md"],
    "ls07-combination-treatment-deg": ["differential_expression.csv", "summary.json", "analysis.py", "report.md"],
    "ls07-combination-treatment-mechanism": ["pathway_enrichment.csv", "mechanism_call.json", "resource_manifest.json", "analysis.py", "report.md"],
    "ls08-multiome-column-match": ["column_mapping.csv", "score_matrix.csv", "analysis.py", "report.md"],
    "ls08-enhancer-promoter-integration": ["pair_evidence.csv", "least_supported.json", "analysis.py", "report.md"],
    "ls09-opentrons-sop": ["protocol.py", "transfer_plan.csv", "simulation.txt", "report.md"],
    "ls09-plate-dilution-recovery": ["root_cause.json", "recovery_plan.csv", "analysis.py", "report.md"],
    "ls10-neun-power-analysis": ["power_result.json", "analysis.py", "report.md"],
    "ls10-treatment-response-model": ["model_coefficients.csv", "model_metadata.json", "analysis.py", "report.md"],
}


def _parse_artifact(path: Path) -> None:
    suffix = path.suffix.lower()
    if path.stat().st_size == 0:
        raise ValueError("empty file")
    if suffix == ".json":
        value = json.loads(path.read_text(encoding="utf-8"))
        _reject_nonfinite(value)
    elif suffix in {".csv", ".tsv"}:
        delimiter = "\t" if suffix == ".tsv" else ","
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            if not reader.fieldnames:
                raise ValueError("missing header")
            list(reader)
    elif suffix == ".png":
        if path.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
            raise ValueError("invalid PNG signature")
    else:
        path.read_text(encoding="utf-8")


def _reject_nonfinite(value: Any) -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("non-finite JSON number")
    if isinstance(value, dict):
        for child in value.values():
            _reject_nonfinite(child)
    elif isinstance(value, list):
        for child in value:
            _reject_nonfinite(child)


def _coverage(output: Path, required: list[str]) -> tuple[int, list[str]]:
    failures: list[str] = []
    for name in required:
        path = output / name
        if not path.is_file():
            failures.append(f"MISSING:{name}")
            continue
        try:
            _parse_artifact(path)
        except Exception as exc:  # stable failure code plus audit detail
            failures.append(f"UNPARSEABLE:{name}:{type(exc).__name__}")
    return (10 if not failures else 0), failures


def _script(output: Path, required: list[str]) -> tuple[int, list[str]]:
    script_names = [name for name in required if name in {"analysis.py", "protocol.py"}]
    if not script_names:  # L1 visual health check uses reproducible view metadata in task checker.
        return 0, []
    failures: list[str] = []
    for name in script_names:
        path = output / name
        if not path.is_file():
            failures.append(f"SCRIPT_MISSING:{name}")
            continue
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=name)
        except Exception as exc:
            failures.append(f"SCRIPT_SYNTAX:{name}:{type(exc).__name__}")
            continue
        literals = [node.value for node in ast.walk(tree) if isinstance(node, ast.Constant) and isinstance(node.value, str)]
        if any(value.startswith(("C:\\", "/Users/", "/home/")) for value in literals):
            failures.append(f"SCRIPT_ABSOLUTE_PATH:{name}")
    return (10 if not failures else 0), failures


def _load_scientific_checker(oracle_dir: Path):
    path = oracle_dir / "scientific_checks.py"
    if not path.is_file():
        return None
    spec = importlib.util.spec_from_file_location("task_scientific_checks", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def run(task_id: str) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", required=True)
    parser.add_argument("--out", "--json-out", dest="json_out")
    args = parser.parse_args()
    workspace = Path(args.workspace).resolve()
    output = workspace / "output"
    required = REQUIRED_OUTPUTS[task_id]
    coverage_score, failures = _coverage(output, required)
    script_score, script_failures = _script(output, required)
    failures.extend(script_failures)

    accepted = bool(ACCEPTED)
    if not accepted:
        result = {
            "task_id": task_id,
            "grader_status": "blocked",
            "blocked_reason": "Scientific checker is absent or has not passed 3/3 reference, empty, and wrong-answer acceptance tests.",
            "hardgate_pass": False,
            "scores": {"coverage": coverage_score, "core_science": None, "direction": None, "summary": None, "script": script_score},
            "deterministic_score": None,
            "failure_codes": sorted(failures + ["ORACLE_NOT_ACCEPTED"]),
        }
    else:
        scientific = check(workspace)  # local task checker; never imports submission code
        core = int(scientific["core_science"])
        direction = int(scientific["direction"])
        summary = int(scientific["summary"])
        if not (0 <= core <= 40 and 0 <= direction <= 15 and 0 <= summary <= 5):
            raise ValueError("Scientific checker returned an out-of-range component")
        failures.extend(scientific.get("failure_codes", []))
        score = coverage_score + core + direction + summary + script_score
        hardgate = not failures and bool(scientific.get("hardgate_pass", False))
        result = {
            "task_id": task_id,
            "grader_status": "scored",
            "hardgate_pass": hardgate,
            "scores": {"coverage": coverage_score, "core_science": core, "direction": direction, "summary": summary, "script": script_score},
            "deterministic_score": score,
            "failure_codes": sorted(set(failures)),
            "criteria": scientific.get("criteria", {}),
        }
    payload = json.dumps(result, ensure_ascii=False, indent=2)
    print(payload)
    if args.json_out:
        Path(args.json_out).write_text(payload + "\n", encoding="utf-8")
    return 0 if result["grader_status"] == "scored" else 2
if __name__ == "__main__":
    raise SystemExit(run("ls04-perturbseq-reference-map"))