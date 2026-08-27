#!/usr/bin/env python3
"""Self-contained deterministic Oracle for ls09-plate-dilution-recovery."""
from __future__ import annotations
# Task-specific scientific scoring implementation.
import csv
import json
from pathlib import Path

ACCEPTED = True


def _failure_mode_is_tip_pickup_before_aspirate(value):
    text = str(value or "").strip().lower().replace("-", " ").replace("_", " ")
    tokens = set(text.split())
    tip_pickup = "tip" in tokens and ("pickup" in tokens or {"pick", "up"}.issubset(tokens))
    failed = bool(tokens & {"failed", "failure", "error"})
    before_aspirate = "before" in tokens and bool(tokens & {"aspirate", "aspiration"})
    return tip_pickup and failed and before_aspirate


def check(workspace: Path):
    output = workspace / "output"
    failures = []
    criteria = {}
    try:
        root = json.loads((output / "root_cause.json").read_text(encoding="utf-8"))
    except Exception:
        root = {}
    root_ok = (
        root.get("failed_well") == "B2"
        and _failure_mode_is_tip_pickup_before_aspirate(root.get("failure_mode"))
        and root.get("liquid_moved") is False
        and set(root.get("completed_wells", [])) == {"A1", "A2", "A3", "B1"}
        and set(root.get("recovery_wells", [])) == {"B2", "B3"}
    )

    try:
        with (output / "recovery_plan.csv").open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except Exception:
        rows = []
    required = {"step", "source", "destination", "transfer_uL", "transfer_pipette", "diluent_source", "diluent_uL", "diluent_pipette", "final_concentration", "final_volume_uL"}
    schema_ok = len(rows) == 2 and required.issubset(rows[0] if rows else {})
    plan_ok = schema_ok
    concentration_ok = schema_ok
    pipette_ok = schema_ok
    no_overdraw = schema_ok
    seen = set()
    used = 0.0
    if schema_ok:
        for row in rows:
            try:
                transfer = float(row["transfer_uL"])
                diluent = float(row["diluent_uL"])
                final_c = float(row["final_concentration"])
                final_v = float(row["final_volume_uL"])
            except (TypeError, ValueError):
                plan_ok = concentration_ok = pipette_ok = no_overdraw = False
                continue
            well = row["destination"]
            seen.add(well)
            used += transfer
            if row["source"] != "source:A2" or row["diluent_source"] != "diluent:R1" or well not in {"B2", "B3"} or abs(transfer - 2) > 1e-6 or abs(diluent - 98) > 1e-6:
                plan_ok = False
            if abs(final_v - 100) > 1e-6 or abs(final_c - 0.5) > 1e-6 or abs(25 * transfer / final_v - final_c) > 1e-6:
                concentration_ok = False
            if row["transfer_pipette"] != "P20" or row["diluent_pipette"] != "P300" or not (2 <= transfer <= 20) or not (20 <= diluent <= 300):
                pipette_ok = False
        plan_ok = plan_ok and seen == {"B2", "B3"}
        no_overdraw = no_overdraw and used <= 1998

    report_text = ""
    try:
        report_text = (output / "report.md").read_text(encoding="utf-8").lower()
    except Exception:
        pass
    summary_ok = all(term in report_text for term in ["b2", "b3", "tip pickup", "0.5", "100"])

    criteria.update({"root_cause_from_log": root_ok, "two_well_recovery_contract": plan_ok,
                     "dilution_mass_balance": concentration_ok, "pipette_feasibility": pipette_ok,
                     "source_not_overdrawn": no_overdraw, "report_consistency": summary_ok})
    for name, ok in criteria.items():
        if not ok:
            failures.append("LS09_RECOVERY_" + name.upper())
    core = (14 if root_ok else 0) + (10 if plan_ok else 0) + (10 if concentration_ok else 0) + (6 if pipette_ok and no_overdraw else 0)
    direction = 15 if root_ok and plan_ok and concentration_ok and pipette_ok and no_overdraw else 0
    summary = 5 if summary_ok else 0
    return {"core_science": core, "direction": direction, "summary": summary,
            "hardgate_pass": root_ok and plan_ok and concentration_ok and pipette_ok and no_overdraw,
            "failure_codes": failures, "criteria": criteria}
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
    raise SystemExit(run("ls09-plate-dilution-recovery"))