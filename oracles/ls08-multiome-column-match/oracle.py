#!/usr/bin/env python3
"""Self-contained deterministic Oracle for ls08-multiome-column-match."""
from __future__ import annotations
import base64

_EMBEDDED_FILES_B64 = {
    "grounding-manifest.json": "ew0KICAic2NoZW1hX3ZlcnNpb24iOiAxLA0KICAidGFza19pZCI6ICJsczA4LW11bHRpb21lLWNvbHVtbi1tYXRjaCIsDQogICJzb3VyY2UiOiB7DQogICAgImRhdGFzZXQiOiAiR2VuZW50ZWNoL2NvbXBiaW9iZW5jaC1kYXRhLXYxIiwNCiAgICAiZGF0YXNldF9yZXZpc2lvbiI6ICJjNjczZjA4NTVmY2UwOWQzMjBmMTY3N2YxNjhmNzg2NGVlYzUyYzFhIiwNCiAgICAicXVlc3Rpb25faWQiOiAibXVsdGlvbWUtbWF0Y2gtYXRhYy1ybmEtcTEiLA0KICAgICJyZXRyaWV2ZWRfYXQiOiAiMjAyNi0wOC0yNlQxMTozODowNSswODowMCINCiAgfSwNCiAgImF1dGhvcml6ZWRfaW5wdXRfaW50ZWdyaXR5IjogWw0KICAgIHsibmFtZSI6ICJtdWx0aW9tZS5tYXRjaC5hdGFjLnJuYS5xMS5hdGFjLnRzdi5neiIsICJieXRlcyI6IDE1MjU5MzU1LCAic2hhMjU2IjogIjJiNzJmODRkNGE5Yzc1ZGM3ODYxMWIzYjZmZjcyN2EwZDIxOTBlM2Q5OGY0OWYyZmQwOGU4ZDFlYmZlMDBkMTYifSwNCiAgICB7Im5hbWUiOiAibXVsdGlvbWUubWF0Y2guYXRhYy5ybmEucTEucm5hLnRzdi5neiIsICJieXRlcyI6IDE0MzIzNTIsICJzaGEyNTYiOiAiY2FjZWVkNTRhOTgzZmI3YzBiNThjYjVkMzRhNjFhOTk5YjBlN2RlYzllZmEyZGQwMTc1NDUwZmEwZjA4MDc4NCJ9DQogIF0sDQogICJhbm5vdGF0aW9uIjogew0KICAgICJuYW1lIjogIk1BTkUuR1JDaDM4LnYxLjMucmVmc2VxX2dlbm9taWMuZ3RmLmd6IiwNCiAgICAic291cmNlX3VybCI6ICJodHRwczovL2h1Z2dpbmdmYWNlLmNvL2RhdGFzZXRzL0dlbmVudGVjaC9jb21wYmlvYmVuY2gtZGF0YS12MS9yZXNvbHZlL2M2NzNmMDg1NWZjZTA5ZDMyMGYxNjc3ZjE2OGY3ODY0ZWVjNTJjMWEvZGF0YS9NQU5FLkdSQ2gzOC52MS4zLnJlZnNlcV9nZW5vbWljLmd0Zi5neiIsDQogICAgImxvY2FsX3BhdGgiOiAiZ3JvdW5kaW5nLXJlc291cmNlcy9NQU5FLkdSQ2gzOC52MS4zLnJlZnNlcV9nZW5vbWljLmd0Zi5neiIsDQogICAgImJ5dGVzIjogODQxNjQ3MywNCiAgICAic2hhMjU2IjogImU2OGM1ZTQ5Mjg5MWFkOGU4MjQ2NGZiNzRiOTkyNmY3NTEzNDA3NzlhZGVkNzI3ZDMzNzdiNzVmYmU2Nzg4ZDgiDQogIH0sDQogICJmcm96ZW5fbWF0Y2hpbmdfcnVsZSI6IHsNCiAgICAiZ2VuZV9hY3Rpdml0eSI6ICJVc2luZyBNQU5FIEdSQ2gzOCB2MS4zIGdlbmUgcmVjb3Jkcywgc3VtIHRoZSBzdXBwbGllZCAxMC1rYiBBVEFDIGJpbnMgb3ZlcmxhcHBpbmcgZWFjaCBnZW5lIGJvZHkgZm9yIGV2ZXJ5IEFUQUMgY29sdW1uLiIsDQogICAgImZlYXR1cmVfZmlsdGVyIjogIktlZXAgc2hhcmVkIGFubm90YXRlZCBnZW5lcyB3aXRoIG5vbnplcm8gdmFyaWFuY2UgaW4gYm90aCBtb2RhbGl0aWVzIGFuZCB3aXRoIFJOQSBtZWFuID4wLjEgb3IgQVRBQyBtZWFuID4xOyAxNywyODEgZ2VuZXMgcGFzcy4iLA0KICAgICJ0cmFuc2Zvcm0iOiAiQXBwbHkgbG9nMXAgc2VwYXJhdGVseSB0byBSTkEgVFBNIGFuZCBBVEFDIGdlbmUgYWN0aXZpdHksIHRoZW4gei1zY29yZSBlYWNoIGdlbmUgYWNyb3NzIHRoZSBlaWdodCBwb3B1bGF0aW9ucy4iLA0KICAgICJzaW1pbGFyaXR5IjogIkZvciBlYWNoIFJOQS9BVEFDIGNvbHVtbiBwYWlyLCB1c2UgdGhlIG1lYW4gcHJvZHVjdCBvZiBnZW5lLXdpc2Ugei1zY29yZXMgKFBlYXJzb24tc3R5bGUgc2ltaWxhcml0eSkuIiwNCiAgICAiYXNzaWdubWVudCI6ICJDaG9vc2UgdGhlIG9uZS10by1vbmUgcGVybXV0YXRpb24gbWF4aW1pemluZyB0aGUgc3VtIG9mIHRoZSBlaWdodCBtYXRjaGVkIHNpbWlsYXJpdGllcy4gUmVzb2x2ZSBleGFjdCB0aWVzIGxleGljb2dyYXBoaWNhbGx5IGZvciBkZXRlcm1pbmlzbS4iDQogIH0sDQogICJyZWNvbXB1dGVkX3RydXRoIjogew0KICAgICJybmFfdG9fYXRhY19tYXBwaW5nIjogWzUsIDEsIDQsIDAsIDYsIDMsIDcsIDJdLA0KICAgICJhbnN3ZXJfc3RyaW5nIjogIjU7MTs0OzA7NjszOzc7MiIsDQogICAgImFzc2lnbm1lbnRfc2NvcmUiOiAyLjk5ODE1MDc1NjUsDQogICAgInJvYnVzdG5lc3MiOiBbDQogICAgICB7ImdlbmVfYWN0aXZpdHlfd2luZG93IjoiZ2VuZV9ib2R5IiwiZmxhbmtfYnAiOjAsIm5fZ2VuZXMiOjE3MjgxLCJtYXBwaW5nIjoiNTsxOzQ7MDs2OzM7NzsyIn0sDQogICAgICB7ImdlbmVfYWN0aXZpdHlfd2luZG93IjoiZ2VuZV9ib2R5IiwiZmxhbmtfYnAiOjEwMDAwLCJuX2dlbmVzIjoxNzM0OSwibWFwcGluZyI6IjU7MTs0OzA7NjszOzc7MiJ9LA0KICAgICAgeyJnZW5lX2FjdGl2aXR5X3dpbmRvdyI6ImdlbmVfYm9keSIsImZsYW5rX2JwIjo1MDAwMCwibl9nZW5lcyI6MTczNTAsIm1hcHBpbmciOiI1OzE7NDswOzY7Mzs3OzIifSwNCiAgICAgIHsiZ2VuZV9hY3Rpdml0eV93aW5kb3ciOiJ0c3MiLCJmbGFua19icCI6NTAwMDAsIm5fZ2VuZXMiOjE3MzUwLCJtYXBwaW5nIjoiNTsxOzQ7MDs2OzM7NzsyIn0sDQogICAgICB7ImdlbmVfYWN0aXZpdHlfd2luZG93IjoidHNzIiwiZmxhbmtfYnAiOjEwMDAwMCwibl9nZW5lcyI6MTczNTAsIm1hcHBpbmciOiI1OzE7NDswOzY7Mzs3OzIifQ0KICAgIF0NCiAgfSwNCiAgImNoZWNrZXJfZ3JvdW5kaW5nX3J1bGUiOiAiUmVjb21wdXRlIGEgb25lLXRvLW9uZSBtYXBwaW5nIGZyb20gYXV0aG9yaXplZCBpbnB1dHMgYW5kIHRoZSBwaW5uZWQgYW5ub3RhdGlvbjsgY2FuZGlkYXRlLXJlcG9ydGVkIGNvcnJlbGF0aW9uIG1hdHJpY2VzIG9yIHBlcm11dGF0aW9ucyBhcmUgbm90IHRydXRoLiBSZXF1aXJlIGV4YWN0bHkgZWlnaHQgdW5pcXVlIDAtYmFzZWQgQVRBQyBpbmRpY2VzIGFuZCB0aGUgZnJvemVuIG1hcHBpbmcuIEFjY2VwdCBlcXVpdmFsZW50IGRldGVybWluaXN0aWMgZ2VuZS1hY3Rpdml0eSB3aW5kb3dzIHdoZW4gdGhlaXIgYXNzaWdubWVudCBpcyBpZGVudGljYWwuIg0KfQ0K",
}

def _embedded_bytes(name: str) -> bytes:
    return base64.b64decode(_EMBEDDED_FILES_B64[name])

def _embedded_json(name: str):
    return json.loads(_embedded_bytes(name).decode("utf-8"))

# Task-specific scientific scoring implementation.
import csv
import hashlib
import json
import re
from pathlib import Path

ACCEPTED = True

# Frozen, provenance-bearing truth in grounding-manifest.json.  The digest makes
# accidental or candidate-side replacement of the truth source fail closed.
GROUNDING_SHA256 = "710b7db1226243d789166e1844e26a126a67793905ae04d7f627abf31e5178fe"
INPUT_SHA256 = {
    "multiome.match.atac.rna.q1.atac.tsv.gz": "2b72f84d4a9c75dc78611b3b6ff727a0d2190e3d98f49f2fd08e8d1ebfe00d16",
    "multiome.match.atac.rna.q1.rna.tsv.gz": "caceed54a983fb7c0b58cb5d34a61a999b0e7dec9efa2dd0175450fa0f080784",
}
EXPECTED = {0: 5, 1: 1, 2: 4, 3: 0, 4: 6, 5: 3, 6: 7, 7: 2}
FATAL_GATES = ("grounding_integrity", "exact_bijective_mapping")


def _grounded() -> bool:
    try:
        raw = _embedded_bytes("grounding-manifest.json")
        if hashlib.sha256(raw).hexdigest() != GROUNDING_SHA256:
            return False
        manifest = json.loads(raw)
        integrity = {x["name"]: x["sha256"] for x in manifest["authorized_input_integrity"]}
        truth = manifest["recomputed_truth"]["rna_to_atac_mapping"]
        rule = manifest["frozen_matching_rule"]
        return (
            manifest.get("task_id") == "ls08-multiome-column-match"
            and manifest.get("source", {}).get("dataset_revision")
            == "c673f0855fce09d320f1677f168f7864eec52c1a"
            and integrity == INPUT_SHA256
            and truth == [EXPECTED[i] for i in range(8)]
            and all(k in rule for k in ("gene_activity", "transform", "similarity", "assignment"))
        )
    except Exception:
        return False


def _index(value: object, modality: str) -> int | None:
    text = str(value).strip()
    prefixes = (
        r"(?:(?:rna|population|pop)[\s:_-]*)?"
        if modality == "rna"
        else r"(?:(?:atac[\s:_-]*(?:column|col)?|column|col)[\s:_-]*)?"
    )
    match = re.fullmatch(prefixes + r"([0-7])", text, flags=re.I)
    return int(match.group(1)) if match else None


def _mapping(workspace: Path) -> tuple[dict[int, int], bool]:
    path = workspace / "output" / "column_mapping.csv"
    try:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except Exception:
        return {}, False
    pairs: list[tuple[int | None, int | None]] = [
        (_index(row.get("rna_population", ""), "rna"), _index(row.get("atac_column", ""), "atac"))
        for row in rows
    ]
    valid = (
        len(pairs) == 8
        and all(a is not None and b is not None for a, b in pairs)
        and {a for a, _ in pairs} == set(range(8))
        and {b for _, b in pairs} == set(range(8))
    )
    return ({int(a): int(b) for a, b in pairs} if valid else {}), valid


def _positive(text: str, patterns: tuple[str, ...]) -> bool:
    negation = re.compile(r"\b(?:no|not|never|without|didn['’]?t|fabricat(?:e|ed|ion)|invent(?:ed)?)\b", re.I)
    for clause in re.split(r"[.!?;\n]+", text):
        if not negation.search(clause) and all(re.search(pattern, clause, re.I) for pattern in patterns):
            return True
    return False


def _explanation(text: str) -> dict[str, bool]:
    return {
        "gene_activity_link": _positive(
            text,
            (r"\b(?:atac|peak|bin|accessib)", r"\b(?:gene|tss|promoter)", r"\b(?:overlap|sum|aggregat|link|window)"),
        ),
        "cross_modal_similarity": _positive(
            text,
            (r"\b(?:rna|expression|tpm)", r"\b(?:atac|activity|accessib)", r"\b(?:correl|similar|covari|spearman|pearson|z[ -]?score)"),
        ),
        "global_bijection": _positive(
            text,
            (r"\b(?:one[- ]to[- ]one|biject|hungarian|assignment)", r"\b(?:maximi[sz]|optimal|total|global)"),
        ),
    }


def _report_mapping(text: str) -> bool:
    pairs = {}
    for a, b in re.findall(r"rna[\s:_-]*([0-7])\s*(?:->|→|=|to)\s*atac[\s:_-]*(?:column[\s:_-]*)?([0-7])", text, re.I):
        pairs[int(a)] = int(b)
    if pairs == EXPECTED:
        return True
    for marker in re.finditer(r"\b(?:permutation|mapping|assignment)\b", text, re.I):
        nums = [int(x) for x in re.findall(r"(?<!\d)([0-7])(?!\d)", text[marker.end() : marker.end() + 220])]
        if nums[:8] == [EXPECTED[i] for i in range(8)]:
            return True
    return False


def check(workspace: Path):
    failures: list[str] = []
    grounded = _grounded()
    mapping, bijection = _mapping(Path(workspace))
    correct_pairs = sum(mapping.get(i) == EXPECTED[i] for i in range(8)) if bijection else 0
    exact = bijection and correct_pairs == 8

    report_path = Path(workspace) / "output" / "report.md"
    try:
        report = report_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        report = ""
    explanation = _explanation(report)
    explained = sum(explanation.values())
    report_summary = exact and explained >= 2 and _report_mapping(report)

    core = 5 * correct_pairs if grounded else 0
    direction = (9 + 2 * explained) if grounded and exact else 0
    summary = 5 if grounded and report_summary else 0
    gates = {"grounding_integrity": grounded, "exact_bijective_mapping": exact}

    if not grounded:
        failures.append("GROUNDING_INTEGRITY_FAILED")
    if not bijection:
        failures.append("MAPPING_NOT_BIJECTION")
    elif not exact:
        failures.append("MAPPING_SCIENTIFICALLY_WRONG")
    if exact and explained < 3:
        failures.append("SHARED_SIGNAL_EXPLANATION_INCOMPLETE")
    if exact and not report_summary:
        failures.append("REPORT_TRUTH_SUMMARY_MISSING")

    return {
        "core_science": core,
        "direction": direction,
        "summary": summary,
        "hardgate_pass": all(gates.values()),
        "criteria": {
            "fatal_gates": gates,
            "bijection": bijection,
            "correct_pairs": correct_pairs,
            "expected_pairs": 8,
            "shared_signal_explanation": explanation,
            "report_truth_summary": report_summary,
        },
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
    raise SystemExit(run("ls08-multiome-column-match"))