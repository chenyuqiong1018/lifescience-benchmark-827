"""Deterministic checker for the LS05 pocket-uncertainty local extension."""
from __future__ import annotations

import csv
import json
from pathlib import Path

ACCEPTED = True


def _load_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def check(workspace: Path) -> dict:
    output = workspace / "output"
    gold = json.loads((Path(__file__).parent / "gold.json").read_text(encoding="utf-8"))
    rows = _load_rows(output / "mutation_priorities.csv")
    by_mut = {row.get("mutation", ""): row for row in rows}
    ids_ok = len(rows) == len(by_mut) == 4 and set(by_mut) == set(gold["order"])
    def rank(value):
        try: return int(value)
        except (TypeError, ValueError): return -1
    ranks_ok = ids_ok and sorted(rank(row.get("rank", "")) for row in rows) == [1, 2, 3, 4]
    order_ok = ranks_ok and all(rank(by_mut[mut].get("rank", "")) == index for index, mut in enumerate(gold["order"], 1))

    evidence_checks = []
    decision_checks = []
    for mutation, expected in gold["rows"].items():
        row = by_mut.get(mutation, {})
        evidence_checks.extend([
            row.get("pocket_support") == expected["pocket_support"],
            row.get("confidence_penalty") == expected["confidence_penalty"],
        ])
        decision_checks.append(row.get("decision") == expected["decision"])

    assessment_path = output / "pocket_assessment.json"
    try:
        assessment = json.loads(assessment_path.read_text(encoding="utf-8")) if assessment_path.is_file() else {}
    except (json.JSONDecodeError, OSError):
        assessment = {}
    assessment_ok = (
        assessment.get("pocket_reliable") is False
        and assessment.get("prioritize_mutations") is False
        and isinstance(assessment.get("reason"), str)
        and bool(assessment["reason"].strip())
    )
    core = (10 if order_ok else 0) + round(20 * sum(evidence_checks) / len(evidence_checks)) + (10 if assessment_ok else 0)
    direction = round(15 * sum(decision_checks) / len(decision_checks))
    report_path = output / "report.md"
    report = report_path.read_text(encoding="utf-8", errors="replace").lower() if report_path.is_file() else ""
    prohibited = ("experimentally validated", "measured binding", "measured activity")
    summary_ok = ("low confidence" in report or "低置信" in report) and ("defer" in report or "暂缓" in report) and not any(term in report for term in prohibited)
    failures = []
    if not ids_ok: failures.append("MUTATION_COVERAGE_OR_DUPLICATE")
    if not order_ok: failures.append("ORDER_MISMATCH")
    if not all(evidence_checks): failures.append("CONFIDENCE_EVIDENCE_MISMATCH")
    if not assessment_ok: failures.append("ASSESSMENT_MISMATCH")
    if not all(decision_checks): failures.append("DECISION_MISMATCH")
    if not summary_ok: failures.append("REPORT_OVERCLAIM_OR_INCONSISTENT")
    return {
        "core_science": core,
        "direction": direction,
        "summary": 5 if summary_ok else 0,
        "hardgate_pass": ids_ok and ranks_ok,
        "failure_codes": failures,
        "criteria": {"order": order_ok, "evidence_correct": sum(evidence_checks), "evidence_total": len(evidence_checks), "assessment": assessment_ok, "decisions_correct": sum(decision_checks), "report_consistent": summary_ok},
    }
