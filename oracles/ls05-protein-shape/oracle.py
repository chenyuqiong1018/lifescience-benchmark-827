#!/usr/bin/env python3
"""Self-contained deterministic Oracle for ls05-protein-shape."""
from __future__ import annotations
# Task-specific scientific scoring implementation.
import json
import math
import re
from pathlib import Path


ACCEPTED = True

# Independently derived from protein.shape.q1.pdb (SHA-256 recorded in the
# author result): an XY projection, viewed along Z and rotated a quarter turn,
# has a long vertical spine with terminal and medial arms on the same side.
EXPECTED_LETTER = "F"
ALLOWED_LETTERS = frozenset("BDFHJL NPRTVXZ".replace(" ", ""))

_NEGATION = re.compile(
    r"\b(?:not|no|never|without|neither|isn['’]?t|doesn['’]?t|cannot|can['’]?t|fails?\s+to)\b",
    re.IGNORECASE,
)


def _affirmed(text: str, patterns: tuple[str, ...]) -> bool:
    """Find a positive assertion, excluding locally negated evidence."""
    for pattern in patterns:
        for match in re.finditer(pattern, text, re.IGNORECASE):
            prefix = text[max(0, match.start() - 48) : match.start()]
            # A sentence boundary makes an earlier negation irrelevant.
            prefix = re.split(r"[.;!?\n]", prefix)[-1]
            if not _NEGATION.search(prefix):
                return True
    return False


def _orientation_evidence(notes: str) -> tuple[bool, bool, bool]:
    projection = _affirmed(
        notes,
        (
            r"\bxy[-\s]?plane\b",
            r"\bx[-\s]?y\s+projection\b",
            r"\bproject(?:ed|ion)?\b.{0,28}\b(?:along|down|onto)\b.{0,12}\b[+\-±]?z(?:[-\s]?axis)?\b",
            r"\bview(?:ing)?\s+(?:vector|direction)\b.{0,20}\b0\s*[, ]\s*0\s*[, ]\s*[+\-]?1\b",
            r"\b[+\-±]?z(?:[-\s]?axis)?\b.{0,24}\b(?:line of sight|view direction|orthographic)\b",
        ),
    )

    explicit_turn = _affirmed(
        notes,
        (
            r"\b(?:rotat(?:e|ed|ion)|turn(?:ed)?)\b.{0,20}\b(?:90\s*(?:°|degrees?)|quarter[-\s]?turn)\b",
            r"\b(?:90\s*(?:°|degrees?)|quarter[-\s]?turn)\b.{0,20}\b(?:rotat(?:e|ed|ion)|turn(?:ed)?)\b",
        ),
    )
    axis_mapping = (
        _affirmed(notes, (r"\b(?:pdb\s+)?x(?:[-\s]?axis|\s+extent)?\b.{0,28}\bvertical\b",))
        and _affirmed(notes, (r"\bnegative[-\s]?y|\blow[-\s]?y\b",))
        and _affirmed(notes, (r"\bright(?:ward)?\b|\blateral\b",))
    )
    camera_mapping = (
        projection
        and _affirmed(notes, (r"\bup\s+vector\b.{0,20}\b[+]?1\s*[, ]\s*0\s*[, ]\s*0\b",))
    )
    quarter_turn = explicit_turn or axis_mapping or camera_mapping

    spine = _affirmed(
        notes,
        (r"\b(?:long|main|vertical)\b.{0,22}\b(?:vertical\s+)?(?:spine|stem|backbone|shaft|axis)\b",),
    )
    two_named_arms = (
        _affirmed(notes, (r"\b(?:top|upper|terminal)\b.{0,20}\b(?:arm|branch|lobe|protrusion)\b",))
        and _affirmed(notes, (r"\b(?:mid(?:dle)?|medial|central)\b.{0,20}\b(?:arm|branch|lobe|protrusion)\b",))
    )
    two_counted_arms = _affirmed(
        notes,
        (r"\b(?:two|2)\b.{0,28}\b(?:arms|branches|lobes|protrusions)\b",),
    )
    same_side = _affirmed(
        notes,
        (r"\bsame[-\s]?side\b", r"\bright(?:ward)?\b.{0,20}\b(?:arms?|branches?|lobes?|protrusions?)\b",
         r"\b(?:arms?|branches?|lobes?|protrusions?)\b.{0,20}\bright(?:ward)?\b"),
    )
    topology = (spine and (two_named_arms or two_counted_arms)) or ((two_named_arms or two_counted_arms) and same_side)
    return projection, quarter_turn, topology


def check(workspace: Path):
    criteria: dict[str, bool] = {}
    failures: list[str] = []
    try:
        data = json.loads((Path(workspace) / "output" / "shape_call.json").read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            data = {}
    except Exception:
        data = {}

    raw_letter = data.get("letter")
    letter = str(raw_letter).strip().upper() if isinstance(raw_letter, str) else ""
    notes = data.get("orientation_notes")
    notes = notes.strip() if isinstance(notes, str) else ""

    projection, quarter_turn, topology = _orientation_evidence(notes)
    evidence_count = sum((projection, quarter_turn, topology))

    identity = letter == EXPECTED_LETTER
    allowed_call = len(letter) == 1 and letter in ALLOWED_LETTERS
    grounded_orientation = evidence_count >= 2

    # Exactly two explicitly named fatal scientific gates.
    gate_shape_identity = identity
    gate_grounded_orientation = grounded_orientation

    try:
        confidence = float(data.get("confidence"))
        confidence_valid = math.isfinite(confidence) and 0.0 <= confidence <= 1.0
    except (TypeError, ValueError):
        confidence_valid = False

    criteria.update(
        allowed_shape_call=allowed_call,
        expected_shape_F=identity,
        xy_projection_grounded=projection,
        quarter_turn_grounded=quarter_turn,
        f_topology_grounded=topology,
        confidence_bounded=confidence_valid,
        orientation_notes_present=bool(notes),
        GATE_SHAPE_IDENTITY=gate_shape_identity,
        GATE_GROUNDED_ORIENTATION=gate_grounded_orientation,
    )

    core_science = (32 if identity else 0) + (8 if identity and topology else 0)
    direction = (6 if projection else 0) + (5 if quarter_turn else 0) + (4 if topology else 0)
    summary = 5 if identity and grounded_orientation and confidence_valid and bool(notes) else 0
    hardgate_pass = gate_shape_identity and gate_grounded_orientation

    if not allowed_call:
        failures.append("SHAPE_CALL_INVALID")
    if not gate_shape_identity:
        failures.append("FATAL_SHAPE_IDENTITY_MISMATCH")
    if not gate_grounded_orientation:
        failures.append("FATAL_ORIENTATION_UNGROUNDED")
    if not confidence_valid:
        failures.append("CONFIDENCE_INVALID")
    if direction < 15:
        failures.append("ORIENTATION_INCOMPLETE")

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
    raise SystemExit(run("ls05-protein-shape"))