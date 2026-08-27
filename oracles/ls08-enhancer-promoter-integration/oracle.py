#!/usr/bin/env python3
"""Self-contained deterministic Oracle for ls08-enhancer-promoter-integration."""
from __future__ import annotations
# Task-specific scientific scoring implementation.
import csv
import hashlib
import json
import math
import re
from pathlib import Path


ACCEPTED = True

INPUT_ROOT = Path(__file__).resolve().parents[2] / "inputs" / "ls08-enhancer-promoter-integration"
INPUT_HASHES = {
    "ep.interactions.q1.expr.csv": "5707e13db4d8252a25db70d0b62e9c11fa8f5721e5ee5320c5be501f3b955092",
    "ep.interactions.q1.hic.csv": "662ec81e31f7699cc5eaa7e0803377fef942320ac9a610e121f825cf8607768f",
}
EXPECTED_IDS = tuple(f"EP{i}" for i in range(1, 8))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _truth() -> dict[str, dict[str, float]]:
    """Recompute truth from the two immutable, hash-verified input modalities."""
    for name, digest in INPUT_HASHES.items():
        path = INPUT_ROOT / name
        if not path.is_file() or _sha256(path) != digest:
            raise ValueError(f"authorized input integrity failure: {name}")

    with (INPUT_ROOT / "ep.interactions.q1.hic.csv").open(
        encoding="utf-8-sig", newline=""
    ) as handle:
        hic = list(csv.DictReader(handle))

    def contact_mean(row: dict[str, str]) -> float:
        return sum(float(row[f"count_rep{i}"]) for i in (1, 2, 3)) / 3.0

    background = [row for row in hic if row.get("set") == "background"]
    x = [math.log10(float(row["distance_bp"])) for row in background]
    y = [contact_mean(row) for row in background]
    xbar, ybar = sum(x) / len(x), sum(y) / len(y)
    slope = sum((a - xbar) * (b - ybar) for a, b in zip(x, y)) / sum(
        (a - xbar) ** 2 for a in x
    )
    intercept = ybar - slope * xbar

    truth: dict[str, dict[str, float]] = {}
    for row in hic:
        if row.get("set") != "candidate":
            continue
        mean = contact_mean(row)
        fitted = intercept + slope * math.log10(float(row["distance_bp"]))
        truth[row["pair_id"]] = {
            "contact_mean": mean,
            "contact": mean - fitted,
        }

    with (INPUT_ROOT / "ep.interactions.q1.expr.csv").open(
        encoding="utf-8-sig", newline=""
    ) as handle:
        expression = list(csv.DictReader(handle))
    for pair_id in truth:
        control = [
            float(row["rna_count"])
            for row in expression
            if row.get("pair_id") == pair_id and row.get("condition") == "control"
        ]
        perturbed = [
            float(row["rna_count"])
            for row in expression
            if row.get("pair_id") == pair_id and row.get("condition") == "perturbed"
        ]
        control_mean = sum(control) / len(control)
        perturbed_mean = sum(perturbed) / len(perturbed)
        truth[pair_id]["effect"] = abs(
            math.log2((perturbed_mean + 0.5) / (control_mean + 0.5))
        )

    if set(truth) != set(EXPECTED_IDS):
        raise ValueError("authorized candidate set is not EP1 through EP7")

    def ranks(field: str) -> dict[str, int]:
        ordered = sorted(truth, key=lambda pair_id: (truth[pair_id][field], pair_id))
        return {pair_id: index for index, pair_id in enumerate(ordered, 1)}

    contact_ranks, effect_ranks = ranks("contact"), ranks("effect")
    for pair_id in truth:
        truth[pair_id]["contact_rank"] = float(contact_ranks[pair_id])
        truth[pair_id]["effect_rank"] = float(effect_ranks[pair_id])
        truth[pair_id]["combined"] = (
            contact_ranks[pair_id] + effect_ranks[pair_id]
        ) / 2.0
    return truth


def _finite(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _close(value: object, expected: float) -> bool:
    number = _finite(value)
    return number is not None and math.isclose(
        number, expected, rel_tol=1e-5, abs_tol=5e-5
    )


def _field(row: dict[str, str], names: tuple[str, ...]) -> object:
    lowered = {str(key).strip().lower(): value for key, value in row.items()}
    for name in names:
        if name in lowered:
            return lowered[name]
    return None


def _read_rows(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))
    except Exception:
        return []


def _called_pair(path: Path) -> str:
    try:
        data = json.loads(path.read_text(encoding="utf-8-sig"))
    except Exception:
        return ""
    if not isinstance(data, dict):
        return ""
    for key in ("pair_id", "least_supported", "least_supported_pair", "pair"):
        value = data.get(key)
        if isinstance(value, str):
            match = re.fullmatch(r"\s*(EP[1-7])\s*", value, re.IGNORECASE)
            if match:
                return match.group(1).upper()
    return ""


def _positive_ep6_summary(text: str) -> bool:
    for segment in re.split(r"[.!?\n;]+", text):
        lowered = segment.lower()
        if "ep6" not in lowered:
            continue
        if not re.search(r"least\s+supported|weakest|lowest\s+(?:combined\s+)?support", lowered):
            continue
        if re.search(r"\b(?:not|isn't|isnt|never|incorrect|false)\b", lowered):
            continue
        if re.search(r"(?:rather\s+than|but)\s+ep6", lowered):
            continue
        return True
    return False


def _report_numeric_consistent(text: str, truth: dict[str, dict[str, float]]) -> bool:
    """Reject explicit fabricated EP6 modality values; qualitative prose is allowed."""
    ep6 = truth["EP6"]
    lowered = text.lower()
    if "ep6" not in lowered:
        return True
    numbers = [
        float(token)
        for token in re.findall(r"(?<![A-Za-z])[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", text)
        if token not in {"6", "+6", "-6"}
    ]
    if not numbers:
        return True
    if re.search(r"contact|hi-c|hic|residual", lowered):
        allowed = (ep6["contact"], ep6["contact_mean"], ep6["contact_rank"])
        if not any(any(math.isclose(n, a, rel_tol=2e-3, abs_tol=5e-3) for a in allowed) for n in numbers):
            return False
    if re.search(r"perturb|crispr|expression|log2|effect", lowered):
        allowed = (ep6["effect"], ep6["effect_rank"])
        if not any(any(math.isclose(n, a, rel_tol=2e-3, abs_tol=5e-3) for a in allowed) for n in numbers):
            return False
    return True


def check(workspace: Path):
    failures: list[str] = []
    criteria: dict[str, object] = {}
    try:
        truth = _truth()
        truth_ok = True
    except Exception:
        truth, truth_ok = {}, False

    rows = _read_rows(Path(workspace) / "output" / "pair_evidence.csv")
    by_id: dict[str, dict[str, str]] = {}
    duplicate = False
    for row in rows:
        pair_id = str(_field(row, ("pair_id", "pair", "id")) or "").strip().upper()
        if pair_id in by_id:
            duplicate = True
        by_id[pair_id] = row

    aliases = {
        "contact": ("contact_evidence", "contact_residual", "hic_residual"),
        "effect": ("perturbation_effect", "abs_log2fc", "absolute_log2fc"),
        "combined": ("combined_support", "combined_mean_rank", "mean_rank"),
        "rank": ("rank", "combined_rank", "support_rank"),
    }
    evaluable = (
        not duplicate
        and len(rows) == 7
        and set(by_id) == set(EXPECTED_IDS)
        and all(
            _finite(_field(by_id[pair_id], aliases[field])) is not None
            for pair_id in EXPECTED_IDS
            for field in aliases
        )
    )

    contact_ok = effect_ok = combined_ok = rank_ok = False
    if truth_ok and evaluable:
        contact_ok = all(
            _close(_field(by_id[pair_id], aliases["contact"]), truth[pair_id]["contact"])
            for pair_id in EXPECTED_IDS
        )
        effect_ok = all(
            _close(_field(by_id[pair_id], aliases["effect"]), truth[pair_id]["effect"])
            for pair_id in EXPECTED_IDS
        )
        combined_ok = all(
            _close(_field(by_id[pair_id], aliases["combined"]), truth[pair_id]["combined"])
            for pair_id in EXPECTED_IDS
        )
        reported_rank = {
            pair_id: float(_field(by_id[pair_id], aliases["rank"]))
            for pair_id in EXPECTED_IDS
        }
        rank_ok = reported_rank["EP6"] == min(reported_rank.values()) and all(
            reported_rank[a] < reported_rank[b]
            for a in EXPECTED_IDS
            for b in EXPECTED_IDS
            if truth[a]["combined"] < truth[b]["combined"]
        )

    called = _called_pair(Path(workspace) / "output" / "least_supported.json")
    call_ok = truth_ok and called == "EP6"
    table_direction_ok = combined_ok and truth_ok and min(
        EXPECTED_IDS, key=lambda pair_id: (truth[pair_id]["combined"], pair_id)
    ) == "EP6"

    report_path = Path(workspace) / "output" / "report.md"
    try:
        report = report_path.read_text(encoding="utf-8-sig", errors="replace")
    except Exception:
        report = ""
    positive_summary = _positive_ep6_summary(report)
    dual_modality = bool(
        re.search(r"contact|hi-c|hic|physical", report, re.IGNORECASE)
        and re.search(r"perturb|crispr|expression|log2|effect", report, re.IGNORECASE)
    )
    report_numbers_ok = truth_ok and _report_numeric_consistent(report, truth)

    core = (15 if contact_ok else 0) + (15 if effect_ok else 0)
    core += 6 if combined_ok else 0
    core += 4 if rank_ok else 0
    direction = (12 if call_ok else 0) + (3 if table_direction_ok else 0)
    summary = (3 if positive_summary and report_numbers_ok else 0)
    summary += 2 if positive_summary and dual_modality and report_numbers_ok else 0

    fatal_gates = {
        "fatal_gate_authorized_truth_integrity": truth_ok,
        "fatal_gate_all_pairs_scientifically_evaluable": evaluable,
        "fatal_gate_both_modalities_grounded": contact_ok and effect_ok,
        "fatal_gate_integration_and_call_grounded": combined_ok and call_ok,
    }
    criteria.update(
        fatal_gates,
        contact_residuals_match_authorized_recomputation=contact_ok,
        perturbation_effects_match_authorized_recomputation=effect_ok,
        combined_mean_ranks_match_authorized_recomputation=combined_ok,
        reported_rank_respects_grounded_support_order=rank_ok,
        least_supported_json_is_EP6=call_ok,
        grounded_table_direction_is_EP6=table_direction_ok,
        report_has_non_negated_EP6_conclusion=positive_summary,
        report_treats_modalities_distinctly=dual_modality,
        report_has_no_fabricated_EP6_modality_number=report_numbers_ok,
    )

    if not truth_ok:
        failures.append("AUTHORIZED_TRUTH_INTEGRITY_FAILURE")
    if not evaluable:
        failures.append("PAIR_EVIDENCE_NOT_SCIENTIFICALLY_EVALUABLE")
    if not contact_ok:
        failures.append("CONTACT_EVIDENCE_MISMATCH")
    if not effect_ok:
        failures.append("PERTURBATION_EFFECT_MISMATCH")
    if not combined_ok:
        failures.append("COMBINED_SUPPORT_MISMATCH")
    if not rank_ok:
        failures.append("SUPPORT_RANK_MISMATCH")
    if not call_ok:
        failures.append("LEAST_SUPPORTED_CALL_MISMATCH")
    if not positive_summary or not report_numbers_ok:
        failures.append("REPORT_CONCLUSION_UNSUPPORTED")
    if not dual_modality:
        failures.append("REPORT_MODALITIES_NOT_DISTINCT")

    return {
        "core_science": core,
        "direction": direction,
        "summary": summary,
        "hardgate_pass": all(fatal_gates.values()),
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
    raise SystemExit(run("ls08-enhancer-promoter-integration"))