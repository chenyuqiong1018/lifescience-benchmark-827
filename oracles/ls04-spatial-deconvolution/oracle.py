#!/usr/bin/env python3
"""Self-contained deterministic Oracle for ls04-spatial-deconvolution."""
from __future__ import annotations
# Task-specific scientific scoring implementation.
import csv
import math
import re
from pathlib import Path

ACCEPTED = True

# Independent truth source: active-set non-negative least squares of the raw
# Spot_710-1 expression vector against the six cell-type mean profiles in the
# authorized spatial.sim.tar.gz (archive SHA-256 f290570ad3d3230c0672e3050b54b37f9d89e5673f2a10216ccd359233aa417a).
# The six coefficients below are the NNLS coefficients normalized to sum to one.
EXPECTED = {
    "B_Cell": 0.323040493130,
    "Endothelial": 0.346842122837,
    "Fibroblast_Stroma": 0.0,
    "Macrophage": 0.318705728025,
    "T_Cell": 0.0114116560073,
    "Tumor_Core": 0.0,
}
MAJOR = ("B_Cell", "Endothelial", "Macrophage")
OTHER = ("Fibroblast_Stroma", "T_Cell", "Tumor_Core")

ALIASES = {
    "B_Cell": (r"\bb[-_ ]?cells?\b", r"\bb[-_ ]?lymphocytes?\b"),
    "Endothelial": (r"\bendothelial(?:[-_ ]?cells?)?\b",),
    "Fibroblast_Stroma": (r"\bfibroblasts?\b", r"\bstroma(?:l)?(?:[-_ ]?cells?)?\b"),
    "Macrophage": (r"\bmacrophages?\b",),
    "T_Cell": (r"\bt[-_ ]?cells?\b", r"\bt[-_ ]?lymphocytes?\b"),
    "Tumor_Core": (r"\btumou?r(?:[-_ ]?core|[-_ ]?cells?)?\b", r"\bmalignant(?:[-_ ]?cells?)?\b"),
}


def _canon(value: object) -> str | None:
    s = re.sub(r"[^a-z0-9]+", "", str(value).lower())
    if s.startswith(("bcell", "blymph")):
        return "B_Cell"
    if s.startswith("endothelial"):
        return "Endothelial"
    if s.startswith(("fibro", "stroma")):
        return "Fibroblast_Stroma"
    if s.startswith("macroph"):
        return "Macrophage"
    if s.startswith(("tcell", "tlymph")):
        return "T_Cell"
    if s.startswith(("tumor", "tumour", "malignant", "cancer")):
        return "Tumor_Core"
    return None


def _header(fieldnames: list[str] | None, choices: set[str]) -> str | None:
    for name in fieldnames or []:
        if re.sub(r"[^a-z0-9]+", "", name.lower()) in choices:
            return name
    return None


def _composition(path: Path):
    raw = {k: 0.0 for k in EXPECTED}
    unknown = 0.0
    problems: list[str] = []
    rows = 0
    try:
        text = path.read_text(encoding="utf-8-sig", errors="replace")
        dialect = csv.Sniffer().sniff(text[:4096], delimiters=",;\t")
        reader = csv.DictReader(text.splitlines(), dialect=dialect)
        type_col = _header(reader.fieldnames, {"celltype", "type", "cellidentity", "label"})
        weight_col = _header(reader.fieldnames, {"weight", "proportion", "fraction", "composition", "abundance"})
        if not type_col or not weight_col:
            raise ValueError("required semantic columns not found")
        for row in reader:
            if not any(str(v).strip() for v in row.values() if v is not None):
                continue
            rows += 1
            try:
                value = float(row.get(weight_col, ""))
            except (TypeError, ValueError):
                problems.append("non-numeric weight")
                continue
            if not math.isfinite(value) or value < 0:
                problems.append("non-finite or negative weight")
                continue
            label = _canon(row.get(type_col, ""))
            if label is None:
                unknown += value
            else:
                raw[label] += value
    except Exception as exc:
        problems.append(type(exc).__name__)

    total = sum(raw.values()) + unknown
    usable = rows > 0 and not problems and math.isfinite(total) and total > 0
    vector = {k: (v / total if usable else 0.0) for k, v in raw.items()}
    unknown_fraction = unknown / total if usable else 0.0
    normalized = usable and abs(total - 1.0) <= 0.010000001
    recognized = usable and unknown_fraction <= 0.01
    return vector, total, unknown_fraction, usable, normalized and recognized, problems


def _mention_state(text: str, label: str) -> tuple[bool, bool]:
    positive = negative = False
    for pattern in ALIASES[label]:
        for match in re.finditer(pattern, text, flags=re.I):
            before = text[max(0, match.start() - 40):match.start()].lower()
            after = text[match.end():match.end() + 35].lower()
            pre_neg = bool(re.search(r"(?:\bno\b|\bnot\b|\bwithout\b|\babsent\b|\bmissing\b|\black(?:s|ed|ing)?\b)(?:(?![.!?;]).){0,28}$", before))
            if re.search(r"\bnot\s+only\s*$", before):
                pre_neg = False
            post_neg = bool(re.match(r"[ \t]*(?:(?:is|are|was|were|seems?|appears?)[ \t]*)?(?:not\b|absent\b|missing\b|excluded\b|unsupported\b|negligible\b|zero\b)", after))
            negative |= pre_neg or post_neg
            positive |= not (pre_neg or post_neg)
    return positive, negative


def _reported_number_conflict(text: str, label: str, submitted: float) -> bool:
    number = r"(?P<n>(?:0(?:\.\d+)?|1(?:\.0+)?|\d{1,3}(?:\.\d+)?)\s*%?)"
    for alias in ALIASES[label]:
        patterns = (alias + r"[^\n.!?;]{0,18}?" + number, number + r"[^\n.!?;]{0,18}?" + alias)
        for pattern in patterns:
            for match in re.finditer(pattern, text, flags=re.I):
                token = match.group("n").replace(" ", "")
                value = float(token.rstrip("%")) / (100.0 if token.endswith("%") else 1.0)
                if 0 <= value <= 1 and abs(value - submitted) > 0.15:
                    return True
    return False


def check(workspace: Path):
    vector, total, unknown, usable, valid, problems = _composition(
        Path(workspace) / "output" / "spot_710_composition.csv"
    )
    background = sum(vector[k] for k in OTHER) + unknown
    major_total = sum(vector[k] for k in MAJOR)
    l1 = sum(abs(vector[k] - EXPECTED[k]) for k in EXPECTED) + unknown if usable else 2.0

    # Exactly two explicitly named fatal scientific gates.
    gate_valid = valid
    gate_major = (
        valid
        and all(vector[k] >= 0.18 for k in MAJOR)
        and major_total >= 0.80
        and background <= 0.20
    )
    hardgate = gate_valid and gate_major

    structure_points = 8.0 if valid else 0.0
    major_points = sum(6.0 * max(0.0, 1.0 - abs(vector[k] - EXPECTED[k]) / 0.25) for k in MAJOR) if usable else 0.0
    global_points = 8.0 * max(0.0, 1.0 - l1 / 0.80) if usable else 0.0
    expected_background = sum(EXPECTED[k] for k in OTHER)
    background_points = 6.0 * max(0.0, 1.0 - abs(background - expected_background) / 0.30) if usable else 0.0
    core = max(0, min(40, round(structure_points + major_points + global_points + background_points)))

    recovered = sum(2 for k in MAJOR if vector[k] >= 0.15) if usable else 0
    balanced = 4 if usable and all(0.20 <= vector[k] <= 0.50 for k in MAJOR) else 0
    mixed = 2 if usable and sum(v >= 0.20 for v in vector.values()) >= 2 and max(vector.values()) <= 0.70 else 0
    clean_background = 3 if usable and background <= 0.10 else 0
    direction = recovered + balanced + mixed + clean_background

    report_path = Path(workspace) / "output" / "report.md"
    report = report_path.read_text(encoding="utf-8", errors="replace") if report_path.is_file() else ""
    states = {k: _mention_state(report, k) for k in EXPECTED}
    mixture_stated = bool(re.search(r"\b(?:mix(?:ed|ture)?|admix(?:ed|ture)?|combination|compos(?:ed|ition)|contribut(?:e|es|ion))\b", report, re.I))
    types_supported = all(states[k][0] and not states[k][1] for k in MAJOR)
    unsupported_major_claim = any(
        states[k][0]
        and any(re.search(r"\b(?:major|dominant|substantial|primary|main)\b", report[max(0, m.start()-35):m.end()+35], re.I)
                for alias in ALIASES[k] for m in re.finditer(alias, report, re.I))
        for k in OTHER
    )
    numbers_consistent = not any(_reported_number_conflict(report, k, vector[k]) for k in MAJOR)
    summary = 0
    if hardgate:
        summary += 2 if mixture_stated else 0
        summary += 3 if types_supported and not unsupported_major_claim and numbers_consistent else 0

    criteria = {
        "fatal_gate_valid_normalized_composition": gate_valid,
        "fatal_gate_supported_major_mixture": gate_major,
        "composition_usable": usable,
        "weight_sum": total,
        "unknown_weight_fraction": unknown,
        "submitted_normalized_weights": vector,
        "truth_weights": EXPECTED,
        "l1_distance_to_truth": l1,
        "major_type_points": round(major_points, 4),
        "global_fidelity_points": round(global_points, 4),
        "background_points": round(background_points, 4),
        "direction_major_recovery": recovered,
        "direction_balanced_three_way_mixture": balanced,
        "direction_not_forced_single_type": mixed,
        "direction_unsupported_types_minor": clean_background,
        "report_mixture_stated": mixture_stated,
        "report_expected_types_affirmed": types_supported,
        "report_no_unsupported_major_claim": not unsupported_major_claim,
        "report_numbers_consistent_with_submission": numbers_consistent,
        "parse_problems": problems,
    }
    failures: list[str] = []
    if not gate_valid:
        failures.append("FATAL_INVALID_NORMALIZED_COMPOSITION")
    if not gate_major:
        failures.append("FATAL_SUPPORTED_MAJOR_MIXTURE_NOT_RECOVERED")
    if usable and l1 > 0.30:
        failures.append("QUANTITATIVE_COMPOSITION_MISMATCH")
    if direction < 15:
        failures.append("MIXTURE_DIRECTION_INCOMPLETE")
    if summary < 5:
        failures.append("REPORT_SUMMARY_NOT_SCIENTIFICALLY_SUPPORTED")
    return {
        "core_science": core,
        "direction": direction,
        "summary": summary,
        "hardgate_pass": hardgate,
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
    raise SystemExit(run("ls04-spatial-deconvolution"))