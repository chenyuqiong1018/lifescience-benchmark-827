#!/usr/bin/env python3
"""Self-contained deterministic Oracle for ls01-grna-offtarget-rank."""
from __future__ import annotations
# Task-specific scientific scoring implementation.
import csv
import math
import re
from pathlib import Path


ACCEPTED = True


GUIDES = {
    "g01": {"activity": 0.82, "risk": "high", "mismatches": 1, "region": "coding", "bystanders": 1},
    "g02": {"activity": 0.67, "risk": "low", "mismatches": 3, "region": "intron", "bystanders": 0},
    "g03": {"activity": 0.74, "risk": "medium", "mismatches": 2, "region": "coding", "bystanders": 0},
    "g04": {"activity": 0.59, "risk": "low", "mismatches": 4, "region": "intergenic", "bystanders": 0},
    "g05": {"activity": 0.78, "risk": "high", "mismatches": 1, "region": "coding", "bystanders": 2},
    "g06": {"activity": 0.64, "risk": "low", "mismatches": 3, "region": "intergenic", "bystanders": 0},
}

RISK_WORDS = {
    "low": {"low", "minimal", "lower"},
    "medium": {"medium", "moderate", "intermediate"},
    "high": {"high", "severe", "critical"},
}


def _risk(value: object) -> str:
    text = re.sub(r"[^a-z]+", " ", str(value).lower()).strip()
    words = set(text.split())
    for level, aliases in RISK_WORDS.items():
        if words & aliases:
            return level
    return ""


def _numbered_evidence(text: str, number: int, noun: str) -> bool:
    words = {0: "zero", 1: "one", 2: "two", 3: "three", 4: "four"}
    number_token = rf"(?:{number}|{words.get(number, str(number))})"
    if noun == "mismatch":
        return bool(re.search(rf"\b{number_token}\s*(?:[- ]?mm\b|[- ]?mismatch(?:es)?\b)", text))
    return bool(re.search(rf"\b{number_token}\s+(?:coding\s+)?bystander(?:s)?\b", text))


def _decision_kind(value: object) -> str:
    text = re.sub(r"\s+", " ", str(value).lower()).strip()
    if re.search(r"do not advance|not recommend|no[- ]go|exclude|reject|avoid|deprioriti|hold", text):
        return "negative"
    if re.search(r"caution|conditional|review|secondary|backup", text):
        return "conditional"
    if re.search(r"advance|recommend|select|preferred|primary|lead|shortlist", text):
        return "positive"
    return "unclear"


def _has_tradeoff(text: str) -> bool:
    if re.search(r"trade[ -]?off|balanc(?:e|ing)", text):
        return True
    return bool(
        re.search(r"activit.{0,100}(?:risk|safe)", text, re.DOTALL)
        or re.search(r"(?:risk|safe).{0,100}activit", text, re.DOTALL)
    )


def _explicit_report_top(text: str) -> str:
    patterns = (
        r"(?:top[- ]ranked(?:\s+guide)?|rank(?:ed)?\s*#?\s*1|primary recommendation|top choice|best guide)\s*(?:is|:|-)?\s*`?\b(g0[1-6])\b",
        r"\b(g0[1-6])\b\s+(?:is|as)\s+(?:the\s+)?(?:top[- ]ranked(?:\s+guide)?|rank(?:ed)?\s*#?\s*1|primary recommendation|top choice|best guide)",
    )
    for pattern in patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            return match.group(1).lower()
    return ""


def _read_rows(path: Path) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))
    except Exception:
        return []


def check(workspace: Path):
    output = Path(workspace) / "output"
    rows = _read_rows(output / "ranked_guides.csv")
    failures: list[str] = []

    by_id: dict[str, dict[str, str]] = {}
    ranks: dict[str, float] = {}
    for row in rows:
        guide_id = str(row.get("guide_id", "")).strip().lower()
        if guide_id not in GUIDES:
            continue
        by_id[guide_id] = row
        try:
            ranks[guide_id] = float(row.get("rank", ""))
        except (TypeError, ValueError):
            pass

    activity_ok: dict[str, bool] = {}
    risk_ok: dict[str, bool] = {}
    mismatch_ok: dict[str, bool] = {}
    context_ok: dict[str, bool] = {}
    decision_kinds: dict[str, str] = {}

    for guide_id, expected in GUIDES.items():
        row = by_id.get(guide_id, {})
        try:
            observed_activity = float(row.get("on_target_score", "nan"))
        except (TypeError, ValueError):
            observed_activity = math.nan
        activity_ok[guide_id] = math.isclose(
            observed_activity, expected["activity"], rel_tol=0.0, abs_tol=1e-6
        )
        risk_ok[guide_id] = _risk(row.get("risk_class", "")) == expected["risk"]

        rationale = str(row.get("rationale", "")).lower()
        mismatch_ok[guide_id] = _numbered_evidence(
            rationale, int(expected["mismatches"]), "mismatch"
        )
        region_ok = expected["region"] in rationale or (
            expected["region"] == "coding" and "exon" in rationale
        )
        bystander_ok = expected["bystanders"] == 0 or _numbered_evidence(
            rationale, int(expected["bystanders"]), "bystander"
        )
        context_ok[guide_id] = bool(region_ok and bystander_ok)
        decision_kinds[guide_id] = _decision_kind(row.get("decision", ""))

    # g01 strictly dominates g05 on the supplied evidence: higher activity,
    # the same one-mismatch high-risk coding hit, and fewer coding bystanders.
    pareto_ok = (
        "g01" in ranks and "g05" in ranks and ranks["g01"] < ranks["g05"]
    )

    activity_points = 2 * sum(activity_ok.values())
    risk_points = 2 * sum(risk_ok.values())
    evidence_points = sum(mismatch_ok.values()) + sum(context_ok.values())
    core_science = activity_points + risk_points + evidence_points + (4 if pareto_ok else 0)

    high_safe = {
        guide_id: decision_kinds[guide_id] == "negative" for guide_id in ("g01", "g05")
    }
    medium_safe = decision_kinds["g03"] in {"conditional", "negative"}
    top_guide = min(ranks, key=ranks.get) if ranks else ""
    top_safe = bool(
        top_guide
        and GUIDES[top_guide]["risk"] != "high"
        and decision_kinds[top_guide] in {"positive", "conditional"}
    )

    try:
        report_text = (output / "report.md").read_text(encoding="utf-8", errors="replace").lower()
    except Exception:
        report_text = ""
    combined_text = report_text + "\n" + "\n".join(
        str(row.get("rationale", "")).lower() for row in by_id.values()
    )
    tradeoff_ok = _has_tradeoff(combined_text)

    direction = (
        (6 if all(high_safe.values()) else 0)
        + (3 if medium_safe else 0)
        + (3 if top_safe else 0)
        + (3 if tradeoff_ok else 0)
    )

    report_top = _explicit_report_top(report_text)
    report_top_consistent = bool(report_top and report_top == top_guide)
    report_mismatch_context = bool(
        re.search(r"\b(?:1|2|3|4|one|two|three|four)\s*(?:[- ]?mm\b|[- ]?mismatch(?:es)?\b)", report_text)
        and re.search(r"coding|exon|intron|intergenic", report_text)
    )
    report_tradeoff = _has_tradeoff(report_text)
    summary = (3 if report_top_consistent else 0) + (1 if report_mismatch_context else 0) + (1 if report_tradeoff else 0)

    if not rows:
        failures.append("MISSING_OR_UNREADABLE_RANKING")
    if not all(activity_ok.values()):
        failures.append("ACTIVITY_SCORE_MISMATCH")
    if not all(risk_ok.values()):
        failures.append("RISK_CLASS_MISMATCH")
    if not all(mismatch_ok.values()):
        failures.append("RATIONALE_MISMATCH_EVIDENCE_MISSING")
    if not all(context_ok.values()):
        failures.append("RATIONALE_CONTEXT_EVIDENCE_MISSING")
    if not pareto_ok:
        failures.append("DOMINATED_GUIDE_ORDER")
    if not all(high_safe.values()):
        failures.append("UNSAFE_HIGH_RISK_DECISION")
    if not medium_safe:
        failures.append("MEDIUM_RISK_DECISION_UNQUALIFIED")
    if not top_safe:
        failures.append("UNSAFE_OR_UNCLEAR_TOP_GUIDE")
    if not tradeoff_ok:
        failures.append("TRADEOFF_NOT_STATED")
    if not report_top:
        failures.append("REPORT_TOP_NOT_EXPLICIT")
    elif not report_top_consistent:
        failures.append("REPORT_TOP_CONTRADICTS_RANKING")
    if not report_mismatch_context:
        failures.append("REPORT_SAFETY_EVIDENCE_INCOMPLETE")
    if not report_tradeoff:
        failures.append("REPORT_TRADEOFF_MISSING")

    criteria = {
        "activity_score_matches": activity_ok,
        "risk_class_matches": risk_ok,
        "rationale_mismatch_evidence": mismatch_ok,
        "rationale_context_evidence": context_ok,
        "pareto_g01_ahead_of_g05": pareto_ok,
        "high_risk_decisions_safe": high_safe,
        "medium_risk_decision_qualified": medium_safe,
        "decision_kinds": decision_kinds,
        "top_guide": top_guide,
        "top_is_safety_acceptable": top_safe,
        "explicit_activity_risk_tradeoff": tradeoff_ok,
        "report_top_guide": report_top,
        "report_top_consistent": report_top_consistent,
        "report_has_mismatch_and_region_evidence": report_mismatch_context,
        "report_has_tradeoff": report_tradeoff,
    }
    hardgate_pass = bool(core_science == 40 and direction == 15)
    return {
        "core_science": core_science,
        "direction": direction,
        "summary": summary,
        "hardgate_pass": hardgate_pass,
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
    raise SystemExit(run("ls01-grna-offtarget-rank"))