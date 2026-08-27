#!/usr/bin/env python3
"""Self-contained deterministic Oracle for ls03-cryptic-exon."""
from __future__ import annotations
# Task-specific scientific scoring implementation.
import csv
import re
from pathlib import Path

ACCEPTED = True

# Immutable truth from grounding-manifest.json (Genentech/compbiobench-data-v1,
# revision c673f0855fce09d320f1677f168f7864eec52c1a; verified artifact hashes).
GENE, CHROM = "GNG10", "chr9"
EXON_0 = (111664536, 111664589)
EXON_1 = (111664537, 111664589)
LEFT = (111661715, 111664536, 40)
RIGHT = (111664589, 111666814, 33)


def _read(path: Path, limit: int = 2_000_000) -> str:
    try:
        if not path.is_file() or path.stat().st_size > limit:
            return ""
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""


def _rows(path: Path) -> list[dict[str, str]]:
    text = _read(path)
    if not text.strip():
        return []
    try:
        reader = csv.DictReader(text.splitlines(), delimiter="\t")
        return [
            {str(k).strip().lower(): str(v or "").strip() for k, v in row.items() if k is not None}
            for row in reader
        ]
    except (csv.Error, TypeError):
        return []


def _value(row: dict[str, str], *names: str) -> str:
    for name in names:
        if name in row and row[name] != "":
            return row[name]
    return ""


def _integer(value: object) -> int | None:
    match = re.search(r"(?<![\d.])-?\d+(?![\d.])", str(value).replace(",", ""))
    return int(match.group()) if match else None


def _chrom(value: str) -> str:
    value = value.strip().lower()
    return "chr" + value[3:] if value.startswith("chr") else "chr" + value


def _explicit_convention(row: dict[str, str], report: str) -> str:
    keys = " ".join(row)
    values = " ".join(row.values()).lower()
    text = f"{keys} {values} {report.lower()}"
    if re.search(r"(?:0|zero)[- _]?based", text) and re.search(r"half[- _]?open|exclusive", text):
        return "zero_half_open"
    if re.search(r"(?:1|one)[- _]?based", text) and re.search(r"inclusive", text):
        return "one_inclusive"
    return ""


def _exon_row(rows: list[dict[str, str]], report: str) -> tuple[dict[str, str], bool, bool, bool]:
    for row in rows:
        gene = _value(row, "gene", "hgnc_gene", "hgnc_symbol", "gene_symbol", "symbol").upper()
        chrom = _chrom(_value(row, "chrom", "chromosome", "chr"))
        start = _integer(_value(row, "start", "exon_start", "exon_start_0based", "start_0based", "start_1based"))
        end = _integer(_value(row, "end", "exon_end", "exon_end_0based", "end_0based", "end_1based"))
        convention = _explicit_convention(row, report)
        interval_ok = (start, end) == EXON_0 and convention == "zero_half_open"
        interval_ok |= (start, end) == EXON_1 and convention == "one_inclusive"
        if gene == GENE or (chrom == CHROM and interval_ok):
            return row, gene == GENE, chrom == CHROM, interval_ok
    return {}, False, False, False


def _junction_tuple(row: dict[str, str]) -> tuple[str, int | None, int | None, int | None]:
    chrom = _value(row, "chrom", "chromosome", "chr")
    start = _integer(_value(row, "intron_start", "junction_start", "donor", "start", "left"))
    end = _integer(_value(row, "intron_end", "junction_end", "acceptor", "end", "right"))
    packed = _value(row, "junction", "junction_id", "coordinates", "coord", "locus")
    if packed:
        match = re.search(r"(?:(chr)?([0-9xy]+)[:_])?(\d{6,})\s*[-:]\s*(\d{6,})", packed, re.I)
        if match:
            chrom = chrom or ((match.group(1) or "") + (match.group(2) or ""))
            start, end = int(match.group(3)), int(match.group(4))
    reads = _integer(_value(row, "junction_reads", "split_read_count", "split_reads", "read_count", "reads", "count", "support"))
    return _chrom(chrom), start, end, reads


def _find_junction(rows: list[dict[str, str]], expected: tuple[int, int, int]) -> tuple[bool, bool, dict[str, str]]:
    for row in rows:
        chrom, start, end, reads = _junction_tuple(row)
        gene = _value(row, "gene", "hgnc_gene", "hgnc_symbol", "gene_symbol", "symbol").upper()
        gene_ok = not gene or gene == GENE
        if gene_ok and chrom == CHROM and (start, end) == expected[:2]:
            return True, reads == expected[2], row
    return False, False, {}


def _novel(row: dict[str, str]) -> bool:
    value = _value(row, "novelty", "status", "is_novel", "novel").strip().lower()
    return value in {"novel", "true", "yes", "1", "unannotated", "not_annotated"}


def check(workspace: Path):
    output = Path(workspace) / "output"
    report = _read(output / "report.md")
    report_l = report.lower()
    exon_rows = _rows(output / "cryptic_exon.tsv")
    junction_rows = _rows(output / "junctions.tsv")

    exon, gene_ok, chrom_ok, interval_ok = _exon_row(exon_rows, report)
    left_geom, left_count, left_row = _find_junction(junction_rows, LEFT)
    right_geom, right_count, right_row = _find_junction(junction_rows, RIGHT)
    start = _integer(_value(exon, "start", "exon_start", "exon_start_0based", "start_0based", "start_1based"))
    end = _integer(_value(exon, "end", "exon_end", "exon_end_0based", "end_0based", "end_1based"))
    convention = _explicit_convention(exon, report)
    length_ok = interval_ok and ((end - start == 53) if convention == "zero_half_open" else (end - start + 1 == 53))
    exon_counts = (
        _integer(_value(exon, "left_junction_reads", "left_reads", "left_support")) == LEFT[2]
        and _integer(_value(exon, "right_junction_reads", "right_reads", "right_support")) == RIGHT[2]
    )
    expression = _value(exon, "expression_evidence", "evidence", "expression", "supporting_evidence").lower()
    expression_ok = interval_ok and exon_counts and bool(re.search(r"\b(?:510|high(?:ly)?|express|junction|split)\b", expression))

    annotation_text = " ".join([report_l] + [" ".join(r.values()).lower() for r in junction_rows])
    annotation_ok = "mane" in annotation_text and bool(re.search(r"v?1\.3\b", annotation_text)) and "grch38" in annotation_text
    negated_novelty = bool(re.search(r"\b(?:not|non)[- ]+novel\b|\bpreviously[- ]+annotated\b|\bis[- ]+annotated\b", annotation_text))
    novelty_ok = left_geom and right_geom and (_novel(left_row) and _novel(right_row) or "both" in report_l and "novel" in report_l) and not negated_novelty
    protein_coding_ok = gene_ok and bool(re.search(r"protein[- ]coding", report_l))

    criteria = {
        "truth_gene_gng10": gene_ok,
        "truth_chromosome_chr9": chrom_ok,
        "truth_exon_interval_with_explicit_convention": interval_ok,
        "truth_exon_length_53bp": length_ok,
        "truth_left_junction_geometry": left_geom,
        "truth_left_junction_40_reads": left_count,
        "truth_right_junction_geometry": right_geom,
        "truth_right_junction_33_reads": right_count,
        "truth_expression_evidence": expression_ok,
        "mane_grch38_v1_3_provenance": annotation_ok,
        "both_junctions_correctly_called_novel": novelty_ok,
        "protein_coding_target_conclusion": protein_coding_ok,
    }
    core = sum((8 * gene_ok, 8 * (chrom_ok and interval_ok), 4 * length_ok,
                5 * left_geom, 4 * left_count, 5 * right_geom, 4 * right_count, 2 * expression_ok))
    direction = 5 * annotation_ok + 6 * novelty_ok + 4 * protein_coding_ok

    summary_facts = {
        "report_gene": bool(re.search(r"\bgng10\b", report_l)),
        "report_53bp_interval": "53" in report_l and ("111664536" in report_l or "111664537" in report_l) and "111664589" in report_l,
        "report_support_counts": bool(re.search(r"\b40\b", report_l) and re.search(r"\b33\b", report_l)),
        "report_novelty_and_annotation": "novel" in report_l and "mane" in report_l and "1.3" in report_l and not negated_novelty,
    }
    criteria.update(summary_facts)
    summary = sum((2 * summary_facts["report_gene"], summary_facts["report_53bp_interval"],
                   summary_facts["report_support_counts"], summary_facts["report_novelty_and_annotation"]))

    # Exactly three fatal scientific gates; all other criteria retain partial credit.
    fatal_gates = {
        "FATAL_TRUTH_GENE": gene_ok,
        "FATAL_TRUTH_EVENT_GEOMETRY": chrom_ok and interval_ok and length_ok,
        "FATAL_TWO_JUNCTION_SUPPORT": left_geom and left_count and right_geom and right_count,
    }
    failures = [name for name, passed in fatal_gates.items() if not passed]
    if direction < 15:
        failures.append("DIRECTION_INCOMPLETE")
    if summary < 5:
        failures.append("SUMMARY_INCOMPLETE")
    criteria["fatal_gates"] = fatal_gates
    return {
        "core_science": int(core),
        "direction": int(direction),
        "summary": int(summary),
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
    raise SystemExit(run("ls03-cryptic-exon"))