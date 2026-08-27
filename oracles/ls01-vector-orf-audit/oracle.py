#!/usr/bin/env python3
"""Self-contained deterministic Oracle for ls01-vector-orf-audit."""
from __future__ import annotations
# Task-specific scientific scoring implementation.
import csv
import re
from pathlib import Path

ACCEPTED = True

EXPECTED = {
    "c01": (True, True, True, True),
    "c02": (False, True, False, True),
    "c03": (True, True, True, True),
}
FIELDS = ("frame_ok", "start_ok", "stop_ok", "tag_ok")


def _truth(value):
    value = str(value or "").strip().lower()
    if value in {"true", "t", "yes", "y", "1", "ok", "pass", "compatible", "present", "intact"}:
        return True
    if value in {"false", "f", "no", "n", "0", "fail", "incompatible", "absent", "missing"}:
        return False
    return None


def _rows(path: Path):
    try:
        with path.open(encoding="utf-8-sig", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle)
            rows = {}
            for raw in reader:
                row = {str(k or "").strip().lower(): v for k, v in raw.items()}
                cid = str(row.get("construct_id", "")).strip().lower()
                if cid and cid not in rows:
                    rows[cid] = row
            return rows
    except Exception:
        return {}


def _has(text, pattern):
    return bool(re.search(pattern, str(text or "").lower()))


def _context(text, cid):
    match = re.search(rf"\b{cid}\b(.{{0,360}}?)(?=\bc0[123]\b|$)", text, re.S)
    return (match.group(0) if match else "").lower()


def check(workspace: Path):
    rows = _rows(Path(workspace) / "output" / "construct_audit.csv")
    criteria = {}
    correct = 0
    for cid, expected in EXPECTED.items():
        row = rows.get(cid, {})
        for field, wanted in zip(FIELDS, expected):
            ok = _truth(row.get(field)) is wanted
            criteria[f"{cid}.{field}"] = ok
            correct += int(ok)
    core = round(40 * correct / 12)

    positive = r"\b(pass(?:ed|es)?|ok|valid|acceptable|compliant|clear|good)\b"
    flagged = r"\b(fail(?:ed)?|error|invalid|reject(?:ed)?|non.?compliant|review|warn(?:ing)?|flagged|attention|inconsistent|issue)\b"
    status_checks = {
        "c01.status_direction": _has(rows.get("c01", {}).get("overall_status"), positive),
        "c02.status_direction": _has(rows.get("c02", {}).get("overall_status"), flagged),
        "c03.status_direction": _has(rows.get("c03", {}).get("overall_status"), flagged),
    }
    c01_issues = str(rows.get("c01", {}).get("issues", "")).strip().lower()
    issue_checks = {
        "c01.issues_clean": "c01" in rows and (not c01_issues or bool(re.fullmatch(r"(?:none|n/?a|no issues?|clear|ok|pass)", c01_issues))),
        "c02.issues_frame_and_stop": _has(rows.get("c02", {}).get("issues"), r"frame|triplet|divisib|modulo|multiple\s+of\s+3|length")
        and _has(rows.get("c02", {}).get("issues"), r"stop|terminat|end\s+codon"),
        "c03.issues_claim_conflict": _has(rows.get("c03", {}).get("issues"), r"claim|declar|metadata|annotat|report|out.?of.?frame|disagree|mismatch|inconsisten|conflict|contrar"),
    }
    criteria.update(status_checks, **issue_checks)
    direction = 3 * sum(status_checks.values()) + 2 * sum(issue_checks.values())

    report_path = Path(workspace) / "output" / "report.md"
    try:
        report = report_path.read_text(encoding="utf-8", errors="replace").lower()
    except Exception:
        report = ""
    c01, c02, c03 = (_context(report, cid) for cid in EXPECTED)
    report_checks = {
        "report_all_constructs": all(re.search(rf"\b{cid}\b", report) for cid in EXPECTED),
        "report_c01_pass": _has(c01, positive),
        "report_c02_frame": _has(c02, r"frame|triplet|divisib|multiple\s+of\s+3|length"),
        "report_c02_stop": _has(c02, r"no\s+(?:terminal\s+)?stop|missing\s+(?:a\s+)?(?:terminal\s+)?stop|lacks?\s+(?:a\s+)?(?:terminal\s+)?stop|stop.+(?:absent|missing|false|fail)"),
        "report_c03_claim_conflict": _has(c03, r"claim|declar|metadata|annotat|out.?of.?frame|disagree|mismatch|inconsisten|conflict|contrar"),
    }
    criteria.update(report_checks)
    summary = sum(report_checks.values())

    failures = []
    if core < 40:
        failures.append("SCIENTIFIC_AUDIT_MISMATCH")
    if direction < 15:
        failures.append("AUDIT_DIRECTION_MISMATCH")
    if summary < 5:
        failures.append("REPORT_MISMATCH")
    return {
        "core_science": core,
        "direction": direction,
        "summary": summary,
        "hardgate_pass": core == 40 and direction == 15,
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
    raise SystemExit(run("ls01-vector-orf-audit"))