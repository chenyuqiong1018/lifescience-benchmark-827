#!/usr/bin/env python3
"""Self-contained deterministic Oracle for ls05-low-confidence-pocket."""
from __future__ import annotations
import base64

_EMBEDDED_FILES_B64 = {
    "gold.json": "ew0KICAib3JkZXIiOiBbIlkyMThGIiwgIlIyNDRBIiwgIkwyNjVXIiwgIlMzMDFBIl0sDQogICJyb3dzIjogew0KICAgICJZMjE4RiI6IHsicmFuayI6IDEsICJwb2NrZXRfc3VwcG9ydCI6ICJ1bnN1cHBvcnRlZF9sb3dfY29uZmlkZW5jZSIsICJjb25maWRlbmNlX3BlbmFsdHkiOiAicGxkZHRfbHRfNTA7cGFlX2d0XzEwQSIsICJkZWNpc2lvbiI6ICJkZWZlcl9zdHJ1Y3R1cmVfdmFsaWRhdGlvbiJ9LA0KICAgICJSMjQ0QSI6IHsicmFuayI6IDIsICJwb2NrZXRfc3VwcG9ydCI6ICJ1bnN1cHBvcnRlZF9sb3dfY29uZmlkZW5jZSIsICJjb25maWRlbmNlX3BlbmFsdHkiOiAicGxkZHRfbHRfNTA7cGFlX2d0XzEwQSIsICJkZWNpc2lvbiI6ICJkZWZlcl9zdHJ1Y3R1cmVfdmFsaWRhdGlvbiJ9LA0KICAgICJMMjY1VyI6IHsicmFuayI6IDMsICJwb2NrZXRfc3VwcG9ydCI6ICJ1bnN1cHBvcnRlZF9sb3dfY29uZmlkZW5jZSIsICJjb25maWRlbmNlX3BlbmFsdHkiOiAicGxkZHRfbHRfNTA7cGFlX2d0XzEwQSIsICJkZWNpc2lvbiI6ICJkZWZlcl9zdHJ1Y3R1cmVfdmFsaWRhdGlvbiJ9LA0KICAgICJTMzAxQSI6IHsicmFuayI6IDQsICJwb2NrZXRfc3VwcG9ydCI6ICJjYXV0aW91c19zdXBwb3J0IiwgImNvbmZpZGVuY2VfcGVuYWx0eSI6ICJub25lIiwgImRlY2lzaW9uIjogIm91dF9vZl9zY29wZV9ub25fcG9ja2V0In0NCiAgfSwNCiAgImFzc2Vzc21lbnQiOiB7InBvY2tldF9yZWxpYWJsZSI6IGZhbHNlLCAicHJpb3JpdGl6ZV9tdXRhdGlvbnMiOiBmYWxzZX0NCn0NCg==",
}

def _embedded_bytes(name: str) -> bytes:
    return base64.b64decode(_EMBEDDED_FILES_B64[name])

def _embedded_json(name: str):
    return json.loads(_embedded_bytes(name).decode("utf-8"))

# Task-specific scientific scoring implementation.
"""Deterministic checker for the LS05 pocket-uncertainty local extension."""
import csv
import json
from pathlib import Path

ACCEPTED = True


def _load_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def check(workspace: Path) -> dict:
    output = workspace / "output"
    gold = _embedded_json("gold.json")
    rows = _load_rows(output / "mutation_priorities.csv")
    by_mut = {row.get("mutation", ""): row for row in rows}
    ids_ok = len(rows) == len(by_mut) == 4 and set(by_mut) == set(gold["order"])
    def rank(value):
        try: return int(value)
        except (TypeError, ValueError): return -1
    ranks_ok = ids_ok and sorted(rank(row.get("rank", "")) for row in rows) == [1, 2, 3, 4]
    order_ok = ranks_ok and all(rank(by_mut[mut].get("rank", "")) == index for index, mut in enumerate(gold["order"], 1))

    evidence_checks = []
    decision_checks = []
    for mutation, expected in gold["rows"].items():
        row = by_mut.get(mutation, {})
        evidence_checks.extend([
            row.get("pocket_support") == expected["pocket_support"],
            row.get("confidence_penalty") == expected["confidence_penalty"],
        ])
        decision_checks.append(row.get("decision") == expected["decision"])

    assessment_path = output / "pocket_assessment.json"
    try:
        assessment = json.loads(assessment_path.read_text(encoding="utf-8")) if assessment_path.is_file() else {}
    except (json.JSONDecodeError, OSError):
        assessment = {}
    assessment_ok = (
        assessment.get("pocket_reliable") is False
        and assessment.get("prioritize_mutations") is False
        and isinstance(assessment.get("reason"), str)
        and bool(assessment["reason"].strip())
    )
    core = (10 if order_ok else 0) + round(20 * sum(evidence_checks) / len(evidence_checks)) + (10 if assessment_ok else 0)
    direction = round(15 * sum(decision_checks) / len(decision_checks))
    report_path = output / "report.md"
    report = report_path.read_text(encoding="utf-8", errors="replace").lower() if report_path.is_file() else ""
    prohibited = ("experimentally validated", "measured binding", "measured activity")
    summary_ok = ("low confidence" in report or "低置信" in report) and ("defer" in report or "暂缓" in report) and not any(term in report for term in prohibited)
    failures = []
    if not ids_ok: failures.append("MUTATION_COVERAGE_OR_DUPLICATE")
    if not order_ok: failures.append("ORDER_MISMATCH")
    if not all(evidence_checks): failures.append("CONFIDENCE_EVIDENCE_MISMATCH")
    if not assessment_ok: failures.append("ASSESSMENT_MISMATCH")
    if not all(decision_checks): failures.append("DECISION_MISMATCH")
    if not summary_ok: failures.append("REPORT_OVERCLAIM_OR_INCONSISTENT")
    return {
        "core_science": core,
        "direction": direction,
        "summary": 5 if summary_ok else 0,
        "hardgate_pass": ids_ok and ranks_ok,
        "failure_codes": failures,
        "criteria": {"order": order_ok, "evidence_correct": sum(evidence_checks), "evidence_total": len(evidence_checks), "assessment": assessment_ok, "decisions_correct": sum(decision_checks), "report_consistent": summary_ok},
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
    raise SystemExit(run("ls05-low-confidence-pocket"))