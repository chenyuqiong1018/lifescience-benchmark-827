#!/usr/bin/env python3
"""Self-contained deterministic Oracle for ls10-neun-power-analysis."""
from __future__ import annotations
# Task-specific scientific scoring implementation.
import json
import math
import re
from pathlib import Path

ACCEPTED = True  # BixBench range/exact gold; acceptance suite in scripts/acceptance-ls06-ls10.py.
READY_CANDIDATE = True


EXPECTED = {"KD_mean": 214.5, "CTRL_mean": 210.625, "KD_sd": 10.9414023651,
            "CTRL_sd": 22.8531804601, "pooled_sd": 17.9162236933, "d": 0.2162844172}


def close(v, e, tol=1e-5):
    try: return math.isclose(float(v), e, rel_tol=tol, abs_tol=tol)
    except (TypeError, ValueError): return False


def pair(obj, labels, key):
    value = obj.get(key, {})
    if isinstance(value, dict): return value.get(labels[0]), value.get(labels[1])
    if isinstance(value, list) and len(value) == 2: return value
    return None, None


def check(workspace: Path):
    failures, criteria = [], {}
    try: data = json.loads((workspace / "output" / "power_result.json").read_text(encoding="utf-8"))
    except Exception: data = {}
    labels = data.get("group_labels", [])
    valid_labels = isinstance(labels, list) and len(labels) == 2 and set(map(str.upper, labels)) == {"KD", "CTRL"}
    labels = labels if valid_labels else ["KD", "CTRL"]
    m1,m2=pair(data,labels,"means"); s1,s2=pair(data,labels,"sds")
    mapped_m={str(labels[0]).upper():m1,str(labels[1]).upper():m2}; mapped_s={str(labels[0]).upper():s1,str(labels[1]).upper():s2}
    tests={"labels":valid_labels,
           "means":close(mapped_m.get("KD"),EXPECTED["KD_mean"]) and close(mapped_m.get("CTRL"),EXPECTED["CTRL_mean"]),
           "sds":close(mapped_s.get("KD"),EXPECTED["KD_sd"]) and close(mapped_s.get("CTRL"),EXPECTED["CTRL_sd"]),
           "pooled_sd":close(data.get("pooled_sd"),EXPECTED["pooled_sd"]),
           "d":close(abs(float(data.get("cohens_d",math.nan))),EXPECTED["d"],tol=5e-3),
           "alpha":close(data.get("alpha"),.05), "power":close(data.get("power"),.8),
           "alternative":str(data.get("alternative","")).lower() in {"two-sided","two sided","two_sided","two.sided"},
           "required_n":data.get("required_n_per_group")==337}
    core=(8 if tests["means"] else 0)+(8 if tests["sds"] else 0)+(6 if tests["pooled_sd"] else 0)+(8 if tests["d"] else 0)+(10 if tests["required_n"] else 0)
    decision_ok=tests["alpha"] and tests["power"] and tests["alternative"] and tests["required_n"]
    report=workspace/"output"/"report.md"; text=report.read_text(encoding="utf-8",errors="replace").lower() if report.is_file() else ""
    summary_ok=bool(re.search(r"0\.21[56]",text) and re.search(r"\b337\b",text) and "0.05" in text and ("0.8" in text or "80%" in text))
    criteria.update(tests,report_consistent=summary_ok)
    if not tests["labels"]:failures.append("GROUP_LABELS_INVALID")
    if core<40:failures.append("POWER_CALCULATION_MISMATCH")
    if not decision_ok:failures.append("POWER_SPECIFICATION_MISMATCH")
    if not summary_ok:failures.append("REPORT_MISMATCH")
    return {"core_science":core,"direction":15 if decision_ok else 0,"summary":5 if summary_ok else 0,
            "hardgate_pass":tests["labels"] and core==40 and decision_ok,"failure_codes":failures,"criteria":criteria}
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
    raise SystemExit(run("ls10-neun-power-analysis"))