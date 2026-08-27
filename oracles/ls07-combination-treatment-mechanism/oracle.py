#!/usr/bin/env python3
"""Self-contained deterministic Oracle for ls07-combination-treatment-mechanism."""
from __future__ import annotations
# Task-specific scientific scoring implementation.
import csv
import json
import math
import re
from pathlib import Path

ACCEPTED = True

CANONICAL_ID = "R-HSA-6791312"
EXPECTED_P = 0.00014600657788625928  # Hypergeometric tail: N=10489, K=336, M=49, x=8.
REQUIRED_COLUMNS = ("pathway_id", "pathway_name", "overlap", "p_value", "padj", "direction")
ALIASES = {
    "pathway_id": ("pathway_id", "reactome_id", "term_id", "id"),
    "pathway_name": ("pathway_name", "pathway", "term", "term_name"),
    "overlap": ("overlap", "overlap_fraction", "gene_ratio", "genes_ratio"),
    "p_value": ("p_value", "pvalue", "p_val", "p"),
    "padj": ("padj", "adjusted_p_value", "adj_p", "q_value", "fdr"),
    "direction": ("direction", "regulation", "effect_direction"),
}


def _key(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def _read_rows(path: Path) -> list[dict[str, str]]:
    try:
        sample = path.read_text(encoding="utf-8-sig", errors="strict")
        dialect = csv.Sniffer().sniff(sample[:4096], delimiters=",\t;")
        raw = list(csv.DictReader(sample.splitlines(), dialect=dialect))
    except Exception:
        return []
    rows: list[dict[str, str]] = []
    for source in raw:
        normalized = {_key(k): ("" if v is None else str(v).strip()) for k, v in source.items() if k}
        row = {}
        for target, aliases in ALIASES.items():
            row[target] = next((normalized[a] for a in aliases if a in normalized), "")
        rows.append(row)
    return rows


def _flatten(value: object) -> str:
    if isinstance(value, dict):
        return " ".join(_flatten(v) for v in value.values())
    if isinstance(value, list):
        return " ".join(_flatten(v) for v in value)
    return str(value)


def _mechanism_match(text: str) -> bool:
    t = text.lower()
    positive = bool(re.search(r"\b(?:tp53|p53)\b", t) and re.search(r"cell[\s_-]*cycle", t))
    negated = bool(
        re.search(r"\b(?:not|no|without)\s+(?:an?\s+)?(?:tp53|p53)\b", t)
        or re.search(r"\b(?:tp53|p53)\b.{0,35}\b(?:is|was|are|were)?\s*not\s+(?:enriched|supported|the\s+mechanism)", t)
    )
    return positive and not negated


def _number(value: object) -> float | None:
    try:
        result = float(str(value).strip())
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def _overlap_is_official(value: str) -> bool:
    match = re.fullmatch(r"\s*(\d+)\s*/\s*(\d+)\s*", value)
    return bool(match and (int(match.group(1)), int(match.group(2))) == (8, 49))


def _direction_supported(text: str) -> bool:
    t = text.lower()
    supported = bool(re.search(r"\b(mixed|bidirectional|dysregulat|cell[\s_-]*cycle\s+arrest|arrest|repress|suppress|inhibit|downregulat|decreas(?:e|ed|ing))\b", t))
    proliferation_claim = bool(
        re.search(r"\b(?:increase[sd]?|promote[sd]?|activate[sd]?|upregulat\w*)\s+(?:the\s+)?(?:cell[\s_-]*cycle|proliferation)\b", t)
        or re.search(r"\b(?:cell[\s_-]*cycle|proliferation)\s+(?:is\s+)?(?:increase[sd]?|promote[sd]?|activate[sd]?|upregulat\w*)\b", t)
    )
    return supported and not proliferation_claim


def _causal_overclaim(text: str) -> bool:
    for sentence in re.split(r"[.!?\n]+", text.lower()):
        if not re.search(r"\b(?:prove[sd]?|demonstrate[sd]?|establish(?:es|ed)?|confirm(?:s|ed)?)\b", sentence):
            continue
        if re.search(r"\b(?:does|do|did|can|cannot|can't|is|was)\s+not\s+(?:prove|demonstrate|establish|confirm)|\b(?:cannot|can't)\s+(?:prove|demonstrate|establish|confirm)", sentence):
            continue
        if re.search(r"\b(?:caus\w*|mechanism|mediate[sd]?|drive[sn]?|responsible)\b", sentence):
            return True
    return False


def _causal_caveat(text: str) -> bool:
    t = text.lower()
    return bool(
        re.search(r"\b(?:does|do|did|can)\s+not\s+(?:prove|demonstrate|establish)\b", t)
        or re.search(r"\b(?:cannot|can't)\s+(?:prove|demonstrate|establish)\b", t)
        or re.search(r"\b(?:association|hypothesis|consistent\s+with|supports?\s+but|not\s+causation|not\s+causal)\b", t)
    )


def check(workspace: Path):
    out = Path(workspace) / "output"
    rows = _read_rows(out / "pathway_enrichment.csv")
    try:
        call = json.loads((out / "mechanism_call.json").read_text(encoding="utf-8"))
    except Exception:
        call = {}
    try:
        report = (out / "report.md").read_text(encoding="utf-8", errors="replace")
    except Exception:
        report = ""

    call_text = _flatten(call)
    canonical = next(
        (
            row for row in rows
            if CANONICAL_ID.lower() in (row["pathway_id"] + " " + row["pathway_name"]).lower()
            and "tp53" in row["pathway_name"].lower()
            and bool(re.search(r"cell[\s_-]*cycle", row["pathway_name"].lower()))
        ),
        None,
    )
    overlap_ok = bool(canonical and _overlap_is_official(canonical["overlap"]))
    p = _number(canonical["p_value"]) if canonical else None
    q = _number(canonical["padj"]) if canonical else None
    stats_ok = bool(
        p is not None and q is not None
        and math.isclose(p, EXPECTED_P, rel_tol=0.10, abs_tol=1e-8)
        and p <= q <= 0.05
    )
    mechanism_ok = _mechanism_match(call_text)
    direction_ok = bool(canonical and _direction_supported(canonical["direction"] + " " + call_text))
    overclaim = _causal_overclaim(report + " " + call_text)
    report_ok = bool(
        report
        and _mechanism_match(report)
        and "enrich" in report.lower()
        and _causal_caveat(report)
        and not overclaim
    )

    table_parseable = bool(rows)
    canonical_ok = canonical is not None
    core = (10 if canonical_ok else 0) + (10 if overlap_ok else 0) + (8 if stats_ok else 0) + (12 if mechanism_ok else 0)
    direction = 15 if canonical_ok and overlap_ok and direction_ok else 0
    summary = 5 if report_ok else 0

    gates = {
        "FATAL_GROUNDED_ENRICHMENT_EVIDENCE": bool(canonical_ok and overlap_ok and stats_ok),
        "FATAL_PRIMARY_MECHANISM_TRUTH": mechanism_ok,
        "FATAL_NO_CAUSAL_OVERCLAIM": not overclaim,
    }
    criteria = {
        "table_parseable": table_parseable,
        "canonical_reactome_pathway": canonical_ok,
        "official_overlap_8_of_49": overlap_ok,
        "grounded_hypergeometric_p_and_valid_fdr": stats_ok,
        "primary_tp53_cell_cycle_mechanism": mechanism_ok,
        "mixed_or_cell_cycle_repressive_direction": direction_ok,
        "report_scientifically_consistent": report_ok,
        "fatal_gates": gates,
    }
    failures = []
    if not table_parseable: failures.append("ENRICHMENT_TABLE_UNPARSEABLE")
    if not canonical_ok: failures.append("CANONICAL_PATHWAY_NOT_FOUND")
    if canonical_ok and not overlap_ok: failures.append("OFFICIAL_OVERLAP_MISMATCH")
    if canonical_ok and not stats_ok: failures.append("ENRICHMENT_STATISTICS_INVALID")
    if not mechanism_ok: failures.append("PRIMARY_MECHANISM_MISMATCH")
    if not direction_ok: failures.append("DIRECTION_UNSUPPORTED")
    if overclaim: failures.append("CAUSAL_OVERCLAIM")
    if not report_ok: failures.append("REPORT_SCIENCE_INCOMPLETE")
    return {
        "core_science": core,
        "direction": direction,
        "summary": summary,
        "hardgate_pass": all(gates.values()),
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
    raise SystemExit(run("ls07-combination-treatment-mechanism"))