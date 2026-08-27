#!/usr/bin/env python3
"""Self-contained deterministic Oracle for ls01-primer-transcript-audit."""
from __future__ import annotations
# Task-specific scientific scoring implementation.
import csv
import re
from pathlib import Path


ACCEPTED = True


EXPECTED = {
    "p01": ({"TX_CANONICAL"}, 102),
    "p02": ({"TX_ALT"}, 99),
    "p03": (set(), None),
}


def _key(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


ALIASES = {
    "pair_id": {"pairid", "id", "primerpair", "primerpairid"},
    "transcripts_matched": {
        "transcriptsmatched", "matchedtranscripts", "transcriptmatches",
        "transcriptmatched",
    },
    "amplicon_length": {
        "ampliconlength", "ampliconbp", "productlength", "productbp",
    },
    "cds_compatible": {"cdscompatible", "cdscompatibility"},
    "status": {"status", "assessment", "disposition", "result"},
    "reason": {"reason", "rationale", "explanation", "notes"},
}


def _read_rows(path: Path) -> tuple[dict[str, dict[str, str]], bool]:
    try:
        raw = path.read_text(encoding="utf-8-sig", errors="replace")
    except Exception:
        return {}, False
    if not raw.strip():
        return {}, False
    try:
        dialect = csv.Sniffer().sniff(raw[:4096], delimiters=",;\t")
        parsed = list(csv.DictReader(raw.splitlines(), dialect=dialect))
    except Exception:
        return {}, False
    if not parsed or not parsed[0]:
        return {}, False
    header_map: dict[str, str] = {}
    for original in parsed[0]:
        normalized = _key(original)
        for canonical, names in ALIASES.items():
            if normalized in names and canonical not in header_map:
                header_map[canonical] = original
    if not set(ALIASES).issubset(header_map):
        return {}, False
    rows: dict[str, dict[str, str]] = {}
    for raw_row in parsed:
        row = {
            canonical: str(raw_row.get(original, "") or "").strip()
            for canonical, original in header_map.items()
        }
        pair_id = row["pair_id"].lower()
        if pair_id in EXPECTED and pair_id not in rows:
            rows[pair_id] = row
    return rows, True


def _transcripts(value: str) -> set[str]:
    text = value.strip()
    if not text or _key(text) in {
        "none", "na", "n/a", "nomatch", "nomatches", "notfound",
        "noamplicon", "null", "zero", "0",
    }:
        return set()
    named = set(re.findall(r"\bTX_[A-Z0-9_]+\b", text.upper()))
    if named:
        return named
    return {
        token.upper()
        for token in re.split(r"[|,;/\s]+", text)
        if token and _key(token) not in {"none", "na", "n/a"}
    }


def _integer(value: str) -> int | None:
    match = re.search(r"(?<!\d)(\d+)(?!\d)", value)
    return int(match.group(1)) if match else None


def _no_product_length(value: str) -> bool:
    if not value.strip():
        return True
    normalized = _key(value)
    return normalized in {
        "0", "none", "na", "n/a", "noamplicon", "noproduct", "notapplicable",
        "notdetected", "null",
    }


def _nonaffirmative_cds(value: str) -> bool:
    normalized = _key(value)
    if not normalized:
        return False
    if normalized in {
        "false", "no", "unknown", "undetermined", "unavailable", "invalid",
        "na", "n/a", "notassessable", "notapplicable", "cannotassess",
    }:
        return True
    return bool(re.search(
        r"\b(no|not|cannot|can't|unknown|undetermined|invalid|unavailable)\b",
        value.lower(),
    ))


def _metadata_warning(text: str) -> bool:
    lower = re.sub(r"\s+", " ", text.lower())
    has_subject = bool(re.search(r"\b(cds|coding sequence|metadata|annotation|coordinates?)\b", lower))
    has_problem = bool(re.search(
        r"\b(invalid|malformed|inconsistent|impossible|out[- ]of[- ]range|outside|"
        r"exceed(?:s|ed)?|beyond|longer than|shorter than|past the end|does not fit|"
        r"cannot be assessed)\b",
        lower,
    ))
    has_length_context = bool(re.search(
        r"\b(sequence|transcript|length|end|range|coordinate|bounds?)\b", lower
    ))
    has_scope = (
        ("tx_canonical" in lower and "tx_alt" in lower)
        or bool(re.search(r"\b(both|all|each|supplied)\b.{0,80}\btranscripts?\b", lower))
        or bool(re.search(r"\btranscripts?\b.{0,80}\b(both|all|each|supplied)\b", lower))
        or ("700" in lower and "640" in lower)
    )
    return has_subject and has_problem and has_length_context and has_scope


def _mismatch_reason(row: dict[str, str], actual: int, claimed: int) -> bool:
    text = f"{row.get('status', '')} {row.get('reason', '')}".lower()
    explicit_numbers = re.search(rf"\b{actual}\b", text) and re.search(rf"\b{claimed}\b", text)
    semantic = (
        re.search(r"\b(mismatch|discrepan|differ|inconsisten|incorrect|wrong|rather than|"
                  r"not equal|off by|expected|claimed)\w*\b", text)
        and re.search(r"\b(amplicon|product|length|size|bp|expected|claimed)\b", text)
    )
    return bool(explicit_numbers or semantic)


def _no_amplicon_reason(row: dict[str, str]) -> bool:
    text = f"{row.get('status', '')} {row.get('reason', '')}".lower()
    return bool(re.search(
        r"\b(no|none|not|failed to|without)\b.{0,45}\b(amplicon|amplif\w*|product|match\w*)\b"
        r"|\b(amplicon|amplif\w*|product|match\w*)\b.{0,45}\b(no|none|not|absent|failed)\b",
        text,
    ))


def _problem_disposition(row: dict[str, str], reason_ok: bool) -> bool:
    status = row.get("status", "").strip()
    positive = _key(status) in {
        "pass", "passed", "ok", "okay", "valid", "compatible", "success",
        "accepted", "match", "matched", "correct",
    }
    return bool(status) and not positive and reason_ok


def _near(text: str, left: str, right_pattern: str, distance: int = 140) -> bool:
    return bool(re.search(
        rf"(?:{re.escape(left)}.{{0,{distance}}}{right_pattern}|"
        rf"{right_pattern}.{{0,{distance}}}{re.escape(left)})",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    ))


def check(workspace: Path):
    failures: list[str] = []
    rows, audit_readable = _read_rows(Path(workspace) / "output" / "primer_audit.csv")
    report_path = Path(workspace) / "output" / "report.md"
    try:
        report = report_path.read_text(encoding="utf-8", errors="replace") if report_path.is_file() else ""
    except Exception:
        report = ""

    p01 = rows.get("p01", {})
    p02 = rows.get("p02", {})
    p03 = rows.get("p03", {})

    tests = {
        "audit_readable": audit_readable,
        "all_primer_rows_present": set(rows) == set(EXPECTED),
        "p01_transcripts": _transcripts(p01.get("transcripts_matched", "")) == {"TX_CANONICAL"},
        "p01_amplicon_102": _integer(p01.get("amplicon_length", "")) == 102,
        "p02_transcripts": _transcripts(p02.get("transcripts_matched", "")) == {"TX_ALT"},
        "p02_amplicon_99": _integer(p02.get("amplicon_length", "")) == 99,
        "p03_no_amplicon": (
            bool(p03)
            and _transcripts(p03.get("transcripts_matched", "")) == set()
            and _no_product_length(p03.get("amplicon_length", ""))
        ),
    }
    evidence = report + "\n" + "\n".join(
        f"{row.get('reason', '')} {row.get('cds_compatible', '')}" for row in rows.values()
    )
    tests["invalid_cds_metadata_identified"] = (
        all(_nonaffirmative_cds(rows.get(pid, {}).get("cds_compatible", "")) for pid in EXPECTED)
        and _metadata_warning(evidence)
    )

    core = (
        (7 if tests["p01_transcripts"] else 0)
        + (6 if tests["p01_amplicon_102"] else 0)
        + (7 if tests["p02_transcripts"] else 0)
        + (6 if tests["p02_amplicon_99"] else 0)
        + (8 if tests["p03_no_amplicon"] else 0)
        + (6 if tests["invalid_cds_metadata_identified"] else 0)
    )

    tests["p01_disposition"] = _problem_disposition(p01, _mismatch_reason(p01, 102, 108))
    tests["p02_disposition"] = _problem_disposition(p02, _mismatch_reason(p02, 99, 104))
    tests["p03_disposition"] = _problem_disposition(p03, _no_amplicon_reason(p03))
    direction = 5 * sum(bool(tests[name]) for name in (
        "p01_disposition", "p02_disposition", "p03_disposition"
    ))

    report_lower = report.lower()
    tests["report_consistent"] = bool(
        _near(report_lower, "p01", r"\b102\s*(?:bp|base pairs?)?\b")
        and _near(report_lower, "p02", r"\b99\s*(?:bp|base pairs?)?\b")
        and _near(report_lower, "p03", r"\b(?:no|not|none|failed|without)\b.{0,40}\b(?:amplicon|amplif\w*|product|match\w*)\b")
        and _metadata_warning(report)
    )
    summary = 5 if tests["report_consistent"] else 0

    if not audit_readable or not tests["all_primer_rows_present"]:
        failures.append("MISSING_OR_UNREADABLE_AUDIT")
    if not tests["p01_transcripts"] or not tests["p02_transcripts"] or not tests["p03_no_amplicon"]:
        failures.append("PRIMER_MATCH_MISMATCH")
    if not tests["p01_amplicon_102"] or not tests["p02_amplicon_99"]:
        failures.append("AMPLICON_LENGTH_MISMATCH")
    if not tests["invalid_cds_metadata_identified"]:
        failures.append("CDS_METADATA_AUDIT_MISSING")
    if direction < 15:
        failures.append("DISPOSITION_MISMATCH")
    if not tests["report_consistent"]:
        failures.append("REPORT_MISMATCH")

    return {
        "core_science": core,
        "direction": direction,
        "summary": summary,
        "hardgate_pass": core == 40 and direction == 15,
        "criteria": tests,
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
    raise SystemExit(run("ls01-primer-transcript-audit"))