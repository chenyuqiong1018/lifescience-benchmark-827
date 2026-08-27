"""Deterministic checker for the LS05 model-ranking local extension."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

ACCEPTED = True


def _load_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _close(value: str, target: float) -> bool:
    try:
        return math.isclose(float(value), target, rel_tol=0, abs_tol=1e-6)
    except (TypeError, ValueError):
        return False


def _rank(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return -1


def check(workspace: Path) -> dict:
    output = workspace / "output"
    gold = json.loads((Path(__file__).parent / "gold.json").read_text(encoding="utf-8"))
    rows = _load_rows(output / "model_ranking.csv")
    by_id = {row.get("model_id", ""): row for row in rows}
    ids_ok = len(rows) == len(by_id) == 3 and set(by_id) == set(gold["order"])
    ranks_ok = ids_ok and sorted(_rank(row.get("rank", "")) for row in rows) == [1, 2, 3]
    order_ok = ranks_ok and all(_rank(by_id[mid].get("rank", "")) == index for index, mid in enumerate(gold["order"], 1))

    metric_checks = []
    decision_checks = []
    for model_id, expected in gold["rows"].items():
        row = by_id.get(model_id, {})
        metric_checks.extend([
            _close(row.get("global_score", ""), expected["global_score"]),
            _close(row.get("interface_score", ""), expected["interface_score"]),
            _close(row.get("critical_residue_risk", ""), expected["critical_residue_risk"]),
        ])
        decision_checks.append(row.get("decision") == expected["decision"])

    core = (20 if order_ok else 0) + round(20 * sum(metric_checks) / len(metric_checks))
    direction = round(15 * sum(decision_checks) / len(decision_checks))
    report_path = output / "report.md"
    report = report_path.read_text(encoding="utf-8", errors="replace").lower() if report_path.is_file() else ""
    summary_ok = "model_a" in report and ("incomplete" in report or "mapping" in report) and "model_b" in report
    failures = []
    if not ids_ok: failures.append("MODEL_COVERAGE_OR_DUPLICATE")
    if not order_ok: failures.append("ORDER_MISMATCH")
    if not all(metric_checks): failures.append("METRIC_MISMATCH")
    if not all(decision_checks): failures.append("DECISION_MISMATCH")
    if not summary_ok: failures.append("REPORT_INCONSISTENT")
    return {
        "core_science": core,
        "direction": direction,
        "summary": 5 if summary_ok else 0,
        "hardgate_pass": ids_ok and ranks_ok,
        "failure_codes": failures,
        "criteria": {"order": order_ok, "metrics_correct": sum(metric_checks), "metrics_total": len(metric_checks), "decisions_correct": sum(decision_checks), "report_consistent": summary_ok},
    }
