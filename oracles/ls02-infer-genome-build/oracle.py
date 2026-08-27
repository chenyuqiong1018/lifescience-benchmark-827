#!/usr/bin/env python3
"""Self-contained deterministic Oracle for ls02-infer-genome-build."""
from __future__ import annotations
# Task-specific scientific scoring implementation.
import json
import math
import re
from pathlib import Path
from typing import Any

ACCEPTED = True

# Immutable UCSC-sequence grounding from grounding-manifest.json.  The VCF side
# was independently verified against input SHA-256
# 0f035fb4f35b47be7a52b8d632daeb005eff6016fa1653f3008229d2e46f016e.
ANCHORS = {
    61795: ("G", "C", "G", "A", "T"),
    3752643: ("T", "C", "T", "C", "G"),
    6519054: ("C", "T", "C", "T", "A"),
    10363652: ("G", "A", "G", "A", "T"),
    13750251: ("A", "T", "A", "T", "T"),
    17236905: ("C", "A", "C", "G", "T"),
    20276671: ("T", "C", "T", "A", "G"),
    24373039: ("G", "T", "G", "G", "G"),
    30899532: ("A", "A", "A", "G", "G"),
    36242820: ("T", "T", "T", "T", "T"),
    40667442: ("T", "C", "T", "G", "T"),
    44379292: ("T", "A", "T", "G", "C"),
    47817855: ("G", "A", "G", "G", "A"),
    51196192: ("G", "G", "G", "T", "T"),
    54390600: ("C", "C", "C", "C", "C"),
    57160580: ("C", "T", "C", "A", "T"),
    60495147: ("C", "A", "C", "G", "A"),
    62965185: ("C", None, "C", "T", "T"),
}
ASSEMBLIES = ("hg18", "hg19", "hg38", "hs1")
ALIASES = {
    "hg18": "hg18", "grch36": "hg18", "ncbi36": "hg18",
    "hg19": "hg19", "grch37": "hg19",
    "hg38": "hg38", "grch38": "hg38",
    "hs1": "hs1", "t2t": "hs1", "t2tchm13": "hs1", "chm13": "hs1",
}
FATAL_GATE_NAMES = (
    "FATAL_TRUTH_CONCLUSION",
    "FATAL_GROUNDED_POSITIVE_EVIDENCE",
    "FATAL_DISCRIMINATING_REFERENCE_EVIDENCE",
)


def _norm(value: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def _assembly(value: Any) -> str | None:
    return ALIASES.get(_norm(value))


def _base(value: Any) -> str | None:
    value = str(value).strip().upper() if value is not None else ""
    return value if re.fullmatch(r"[ACGTN]", value) else None


def _position(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)) and float(value).is_integer():
        return int(value)
    match = re.search(r"(?i)(?:chr\s*20\s*[:._-]\s*)?(\d{4,})", str(value))
    return int(match.group(1)) if match else None


def _mapping_base(mapping: Any, assembly: str) -> str | None:
    if not isinstance(mapping, dict):
        return None
    for key, value in mapping.items():
        if _assembly(key) == assembly:
            return _base(value)
    return None


def _collect_claims(root: Any) -> dict[int, dict[str, Any]]:
    claims: dict[int, dict[str, Any]] = {}
    position_keys = {"pos", "position", "coordinate", "coordinate1based", "chrompos", "locus"}
    ref_keys = {"vcfref", "observedref", "variantref", "inputref", "allele"}
    map_keys = {"referencebases", "references", "assemblybases", "refbases", "comparison", "comparisons", "assemblyrefs"}

    def walk(obj: Any, inherited_pos: int | None = None, inherited_ref: str | None = None) -> None:
        if isinstance(obj, list):
            for item in obj:
                walk(item, inherited_pos, inherited_ref)
            return
        if not isinstance(obj, dict):
            return

        normalized = {_norm(k): v for k, v in obj.items()}
        pos = inherited_pos
        for key in position_keys:
            if key in normalized:
                pos = _position(normalized[key]) or pos
                break
        observed = inherited_ref
        for key in ref_keys:
            if key in normalized:
                observed = _base(normalized[key]) or observed
                break
        if "ref" in normalized and not any(k in normalized for k in ("assembly", "referencebuild", "targetbuild")):
            observed = _base(normalized["ref"]) or observed

        if pos in ANCHORS:
            claim = claims.setdefault(pos, {"vcf_refs": set(), "bases": {a: set() for a in ASSEMBLIES}, "negative_hg19": False})
            if observed:
                claim["vcf_refs"].add(observed)
            for key, value in obj.items():
                assembly = _assembly(key)
                if assembly and _base(value):
                    claim["bases"][assembly].add(_base(value))
                if _norm(key) in map_keys and isinstance(value, dict):
                    for assembly_name in ASSEMBLIES:
                        found = _mapping_base(value, assembly_name)
                        if found:
                            claim["bases"][assembly_name].add(found)
                    for akey, aval in value.items():
                        if _assembly(akey) == "hg19" and aval is False:
                            claim["negative_hg19"] = True

            assembly_value = next((normalized[k] for k in ("assembly", "referencebuild", "targetbuild", "build") if k in normalized), None)
            assembly_name = _assembly(assembly_value)
            if assembly_name:
                reported_base = next((_base(normalized[k]) for k in ("referencebase", "refbase", "base") if k in normalized and _base(normalized[k])), None)
                if reported_base:
                    claim["bases"][assembly_name].add(reported_base)
            for key in ("match", "matches", "ismatch", "refmatch"):
                if key in normalized and normalized[key] is False and assembly_name == "hg19":
                    claim["negative_hg19"] = True

        for value in obj.values():
            if isinstance(value, (dict, list)):
                walk(value, pos, observed)

    walk(root)
    return claims


def _evidence_metrics(data: dict[str, Any]) -> dict[str, Any]:
    claims = _collect_claims(data.get("evidence", []))
    verified, discriminating, panels, contradictions = [], [], [], []
    for pos, claim in claims.items():
        expected = dict(zip(("vcf",) + ASSEMBLIES, ANCHORS[pos]))
        refs = claim["vcf_refs"]
        hg19 = claim["bases"]["hg19"]
        ref_ok = expected["vcf"] in refs and not any(v != expected["vcf"] for v in refs)
        hg19_ok = expected["hg19"] in hg19 and not any(v != expected["hg19"] for v in hg19)
        contradicted = claim["negative_hg19"] or any(v != expected["vcf"] for v in refs) or any(v != expected["hg19"] for v in hg19)
        if contradicted:
            contradictions.append(pos)
        if ref_ok and hg19_ok and not contradicted:
            verified.append(pos)
            exact_alternatives = [a for a in ("hg18", "hg38", "hs1") if expected[a] is not None and expected[a] in claim["bases"][a]]
            if any(expected[a] != expected["vcf"] for a in exact_alternatives):
                discriminating.append(pos)
            if all(expected[a] is not None and expected[a] in claim["bases"][a] for a in ASSEMBLIES):
                panels.append(pos)
    return {
        "verified_positions": sorted(verified),
        "discriminating_positions": sorted(discriminating),
        "full_panel_positions": sorted(panels),
        "contradictory_positions": sorted(contradictions),
    }


def _confidence_ok(value: Any, evidence_count: int, contradicted: bool) -> bool:
    if evidence_count < 4 or contradicted:
        return False
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        number = float(value)
        number = number / 100 if number > 1 else number
        return math.isfinite(number) and 0.8 <= number <= 1
    return _norm(value) in {"high", "veryhigh", "strong"}


def check(workspace: Path):
    failures: list[str] = []
    try:
        data = json.loads((workspace / "output" / "build_call.json").read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            data = {}
    except Exception:
        data = {}

    metrics = _evidence_metrics(data)
    verified = len(metrics["verified_positions"])
    discriminating = len(metrics["discriminating_positions"])
    panels = len(metrics["full_panel_positions"])
    contradicted = bool(metrics["contradictory_positions"])
    build_ok = _assembly(data.get("build")) == "hg19"
    span_ok = verified >= 3 and max(metrics["verified_positions"], default=0) - min(metrics["verified_positions"], default=0) >= 20_000_000

    fatal_gates = {
        "FATAL_TRUTH_CONCLUSION": build_ok,
        "FATAL_GROUNDED_POSITIVE_EVIDENCE": verified >= 4 and not contradicted,
        "FATAL_DISCRIMINATING_REFERENCE_EVIDENCE": discriminating >= 2,
    }
    hardgate_pass = all(fatal_gates.values())

    core = (6 if build_ok else 0) + min(verified, 10) * 2 + min(discriminating, 4) * 2 + (4 if panels >= 4 else 0) + (2 if span_ok else 0)
    core = min(40, core)
    evidence_direction = 4 if verified >= 4 and discriminating >= 2 else (2 if verified >= 2 else 0)
    confidence_ok = _confidence_ok(data.get("confidence"), verified, contradicted)
    direction = (9 if build_ok else 0) + evidence_direction + (2 if confidence_ok else 0)

    report_path = workspace / "output" / "report.md"
    report = report_path.read_text(encoding="utf-8", errors="replace") if report_path.is_file() else ""
    report_lower = report.lower()
    correct_conclusion = bool(re.search(
        r"\b(?:supported|inferred|selected|called|concluded)\s+(?:genome\s+)?(?:build|assembly)\s*(?:is|:|=)\s*(?:hg19|grch\s*37)\b",
        report_lower,
    ))
    correct_conclusion |= bool(re.search(
        r"\b(?:build|assembly|call|conclusion)\s*(?:is|:|=)\s*(?:hg19|grch\s*37)\b",
        report_lower,
    ))
    wrong_conclusion = bool(re.search(
        r"\b(?:supported|inferred|selected|called|concluded)\s+(?:genome\s+)?(?:build|assembly)\s*(?:is|:|=)\s*(?:hg18|hg38|hs1|t2t|grch\s*(?:36|38))\b",
        report_lower,
    ))
    wrong_conclusion |= bool(re.search(
        r"\b(?:build|assembly|call|conclusion)\s*(?:is|:|=)\s*(?:hg18|hg38|hs1|t2t|grch\s*(?:36|38))\b",
        report_lower,
    ))
    report_build = correct_conclusion and not wrong_conclusion
    mentioned = sum(bool(re.search(rf"\b{pos}\b[^\n]{{0,100}}\b{ANCHORS[pos][0]}\b", report, re.I)) for pos in metrics["verified_positions"])
    report_method = bool(re.search(r"\b(?:ref|reference allele)\b", report_lower) and re.search(r"\b(?:coordinate|position|locus)\b", report_lower))
    summary = (2 if report_build else 0) + (2 if mentioned >= 2 else 0) + (1 if report_method else 0)

    for name, passed in fatal_gates.items():
        if not passed:
            failures.append(name)
    if core < 40:
        failures.append("PARTIAL_GROUNDED_ANCHOR_COVERAGE")
    if direction < 15:
        failures.append("PARTIAL_DIRECTION_OR_CONFIDENCE")
    if summary < 5:
        failures.append("PARTIAL_REPORT_SUMMARY")

    criteria = {
        "fatal_gates": fatal_gates,
        "build_resolves_to_hg19": build_ok,
        **metrics,
        "verified_anchor_count": verified,
        "discriminating_anchor_count": discriminating,
        "full_reference_panel_count": panels,
        "broad_coordinate_span": span_ok,
        "confidence_calibrated_to_verified_evidence": confidence_ok,
        "report_concludes_hg19": report_build,
        "report_verified_anchor_mentions": mentioned,
        "report_describes_coordinate_ref_method": report_method,
        "candidate_self_reported_counts_and_scores_ignored": True,
    }
    return {
        "core_science": core,
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
    raise SystemExit(run("ls02-infer-genome-build"))