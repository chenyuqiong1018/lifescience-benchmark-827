#!/usr/bin/env python3
"""Self-contained deterministic Oracle for ls02-find-deletion."""
from __future__ import annotations
# Task-specific scientific scoring implementation.
import csv
import json
import math
import re
from pathlib import Path

ACCEPTED = True

# Immutable truth from grounding-manifest.json (input hashes verified by the author).
TRUTH = {
    "chrom": "22",
    "start": 20_000_000,
    "end": 21_000_000,
    "size": 1_000_000,
    "depth_ratio": 0.4922519150926846,
}
TOLERANCE_BP = 100_000
FATAL_GATES = ("FATAL_GROUNDED_LOCUS", "FATAL_GROUNDED_DEPTH_EVIDENCE")


def _norm(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def _number(value: object) -> float:
    text = str(value).strip().lower().replace(",", "")
    match = re.search(r"[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?", text)
    if not match:
        return math.nan
    number = float(match.group())
    if re.search(r"\bmb\b", text):
        number *= 1_000_000
    elif re.search(r"\bkb\b", text):
        number *= 1_000
    return number


def _read_text(path: Path, limit: int = 1_000_000) -> str:
    try:
        if not path.is_file() or path.stat().st_size > limit:
            return ""
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""


ALIASES = {
    "chrom": {"chrom", "chromosome", "chr", "contig"},
    "start": {"start_100kb", "start", "start_bp", "breakpoint_start", "deletion_start"},
    "end": {"end_100kb", "end", "end_bp", "breakpoint_end", "deletion_end"},
    "size": {"size_bp", "size", "length", "length_bp", "deletion_size"},
    "support": {"supporting_signals", "supporting_signal", "support", "signals", "evidence"},
}


def _deletion_rows(path: Path) -> list[dict[str, str]]:
    text = _read_text(path)
    if not text.strip():
        return []
    try:
        dialect = csv.Sniffer().sniff(text[:4096], delimiters="\t,;")
        reader = csv.DictReader(text.splitlines(), dialect=dialect)
        headers = {_norm(h): h for h in (reader.fieldnames or []) if h}
        selected = {key: next((headers[a] for a in names if a in headers), None)
                    for key, names in ALIASES.items()}
        if not all(selected[key] for key in ("chrom", "start", "end", "size")):
            return []
        rows = []
        for raw in reader:
            row = {key: str(raw.get(header, "") or "").strip()
                   for key, header in selected.items() if header}
            if any(row.get(key) for key in ("chrom", "start", "end", "size")):
                rows.append(row)
        return rows
    except (csv.Error, TypeError):
        return []


def _chrom(value: object) -> str:
    return re.sub(r"^chr", "", str(value).strip().lower())


def _close(value: object, expected: float, tolerance: float = TOLERANCE_BP) -> bool:
    number = _number(value)
    return math.isfinite(number) and abs(number - expected) <= tolerance


def _json(path: Path) -> object:
    try:
        if not path.is_file() or path.stat().st_size > 1_000_000:
            return {}
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def _candidate_ratios(obj: object) -> list[float]:
    found: list[float] = []
    if isinstance(obj, dict):
        for key, value in obj.items():
            name = _norm(key)
            if "ratio" in name and any(token in name for token in ("depth", "coverage", "flank")):
                try:
                    number = float(value)
                    if 1 < number <= 100:
                        number /= 100
                    found.append(number)
                except (TypeError, ValueError):
                    found.append(math.nan)
            elif isinstance(value, (dict, list)):
                found.extend(_candidate_ratios(value))
    elif isinstance(obj, list):
        for value in obj:
            found.extend(_candidate_ratios(value))
    return found


def _text_ratios(text: str) -> list[float]:
    values = []
    patterns = (
        r"(?:depth|coverage)(?:[-_ ]to[-_ ]flank)?[-_ ]ratio\s*(?:is|=|:|of)?\s*([0-9]+(?:\.[0-9]+)?%?)",
        r"ratio\s+(?:of\s+)?([0-9]+(?:\.[0-9]+)?%?)\s*(?:for\s+)?(?:depth|coverage)",
    )
    for pattern in patterns:
        for raw in re.findall(pattern, text, flags=re.I):
            number = float(raw.rstrip("%"))
            if raw.endswith("%") or 1 < number <= 100:
                number /= 100
            values.append(number)
    return values


def _evidence_state(text: str, qc: object) -> tuple[bool, bool, bool]:
    lowered = re.sub(r"\s+", " ", text.lower())
    positive = bool(
        re.search(r"(?:depth|coverage).{0,70}(?:deplet|drop|reduc|decreas|deficit|loss|low|lower|half)", lowered)
        or re.search(r"(?:deplet|drop|reduc|decreas|deficit|loss|low|lower|half).{0,70}(?:depth|coverage)", lowered)
    )
    contradicted = bool(
        re.search(r"(?:no|not|without|absent|lacking).{0,30}(?:depth|coverage).{0,35}(?:deplet|drop|reduc|loss|support|evidence|change)", lowered)
        or re.search(r"(?:no|not|without|absent|lacking).{0,55}(?:deplet|drop|reduc|decreas|deficit|loss|low|lower).{0,25}(?:depth|coverage)", lowered)
        or re.search(r"(?:depth|coverage).{0,30}(?:normal|unchanged|not reduced|not depleted|no drop)", lowered)
    )
    ratios = _candidate_ratios(qc) + _text_ratios(lowered)
    numeric_ok = all(math.isfinite(v) and abs(v - TRUTH["depth_ratio"]) <= 0.20 for v in ratios)
    return positive, contradicted, numeric_ok


def check(workspace: Path):
    output = Path(workspace) / "output"
    rows = _deletion_rows(output / "deletion.tsv")
    singleton = len(rows) == 1
    row = rows[0] if singleton else {}

    chrom_ok = singleton and _chrom(row.get("chrom", "")) == TRUTH["chrom"]
    start_ok = singleton and _close(row.get("start", ""), TRUTH["start"])
    end_ok = singleton and _close(row.get("end", ""), TRUTH["end"])
    size_ok = singleton and _close(row.get("size", ""), TRUTH["size"])
    grounded_locus = chrom_ok and start_ok and end_ok and size_ok

    report = _read_text(output / "report.md")
    combined = " ".join((row.get("support", ""), report))
    qc = _json(output / "qc.json")
    positive, contradicted, numeric_ok = _evidence_state(combined, qc)
    manifest_depth_support = TRUTH["depth_ratio"] < 0.65
    grounded_evidence = grounded_locus and manifest_depth_support and positive and not contradicted and numeric_ok

    deletion_claim = bool(re.search(
        r"\bdelet(?:ion|ed)\b|\bheterozygous\s+(?:cnv\s+)?loss\b|\b(?:copy[- ]number|cnv)\s+loss\b",
        combined, flags=re.I,
    ))
    duplication_claim = bool(re.search(r"\bduplication\b|\bcopy[- ]number\s+gain\b", combined, flags=re.I))
    direction_class = grounded_evidence and deletion_claim and not duplication_claim

    precision_ok = bool(re.search(
        r"(?:nearest|rounded|resolution|precision|approx(?:imate|imately)?).{0,45}100\s*[- ]?kb|"
        r"100\s*[- ]?kb.{0,45}(?:nearest|rounded|resolution|precision|approx(?:imate|imately)?)",
        report, flags=re.I,
    ))
    exact_precision_phrase = bool(re.search(
        r"(?:exact|base[- ]pair|single[- ]base).{0,35}(?:breakpoint|precision|resolution)", report, flags=re.I,
    ))
    exact_precision_negated = bool(re.search(
        r"(?:no|not|cannot|can't|do not|does not|without|insufficient).{0,45}(?:exact|base[- ]pair|single[- ]base)",
        report, flags=re.I,
    ))
    false_precision = exact_precision_phrase and not exact_precision_negated
    precision_ok = precision_ok and not false_precision

    core = (8 if chrom_ok else 0) + (8 if start_ok else 0) + (8 if end_ok else 0)
    core += (6 if size_ok else 0) + (10 if grounded_evidence else 0)
    direction = (8 if direction_class else 0) + (7 if grounded_locus and precision_ok else 0)

    compact_report = report.lower().replace(",", "")
    report_locus = "chr22" in compact_report and bool(
        re.search(r"20(?:\.0)?\s*(?:[-–—]|to)\s*21(?:\.0)?\s*mb", compact_report)
        or ("20000000" in compact_report and "21000000" in compact_report)
    )
    report_size = bool(re.search(r"\b1(?:\.0)?\s*mb\b|\b1000000\s*bp\b", compact_report))
    report_consistent = report_locus and report_size and direction_class and precision_ok
    summary = 5 if report_consistent else 0

    criteria = {
        "singleton_deletion_record": singleton,
        "chromosome_matches_ground_truth": chrom_ok,
        "start_within_100kb": start_ok,
        "end_within_100kb": end_ok,
        "size_within_100kb": size_ok,
        "depth_evidence_positive": positive,
        "depth_evidence_not_negated": not contradicted,
        "reported_depth_ratios_concordant": numeric_ok,
        "deletion_direction": direction_class,
        "precision_limit_100kb": precision_ok,
        "report_consistent": report_consistent,
        "fatal_gates": {
            FATAL_GATES[0]: grounded_locus,
            FATAL_GATES[1]: grounded_evidence,
        },
    }
    failures = []
    if not grounded_locus:
        failures.append(FATAL_GATES[0])
    if not grounded_evidence:
        failures.append(FATAL_GATES[1])
    if not direction_class:
        failures.append("DELETION_DIRECTION_UNSUPPORTED")
    if not precision_ok:
        failures.append("PRECISION_LIMIT_NOT_ESTABLISHED")
    if not report_consistent:
        failures.append("REPORT_SCIENCE_MISMATCH")

    return {
        "core_science": core,
        "direction": direction,
        "summary": summary,
        "hardgate_pass": grounded_locus and grounded_evidence,
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
    raise SystemExit(run("ls02-find-deletion"))