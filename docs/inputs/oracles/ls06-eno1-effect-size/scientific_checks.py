from __future__ import annotations

import json
import math
import re
from pathlib import Path

ACCEPTED = True  # Benchmark gold independently encoded; acceptance suite in scripts/acceptance-ls06-ls10.py.
READY_CANDIDATE = True


def close(value, expected, *, rel=1e-6, abs_=1e-9):
    try:
        return math.isclose(float(value), expected, rel_tol=rel, abs_tol=abs_)
    except (TypeError, ValueError):
        return False


def check(workspace: Path):
    failures, criteria = [], {}
    try:
        data = json.loads((workspace / "output" / "eno1_effect.json").read_text(encoding="utf-8"))
    except Exception:
        data = {}
    tests = {
        "gene": str(data.get("gene", "")).upper() == "ENO1",
        "normal_value": close(data.get("normal_value"), 72896133.2946858, rel=5e-6),
        "tumor_value": close(data.get("tumor_value"), 350385456.451912, rel=5e-6),
        "fold_change": close(data.get("fold_change"), 4.81, rel=2e-3),
        "log2_fold_change": close(data.get("log2_fold_change"), 2.27, abs_=0.011),
        "source_file": "proteomic" in str(data.get("source_file", "")).lower(),
        "source_sheet": str(data.get("source_sheet", "")).strip().lower() == "tumor vs normal",
    }
    criteria.update(tests)
    core = sum((10, 10, 10, 10)[i] for i, key in enumerate(("normal_value", "tumor_value", "fold_change", "log2_fold_change")) if tests[key])
    direction_ok = tests["fold_change"] and tests["log2_fold_change"] and float(data.get("fold_change", 0)) > 1 and float(data.get("log2_fold_change", 0)) > 0
    report = (workspace / "output" / "report.md")
    text = report.read_text(encoding="utf-8", errors="replace").lower() if report.is_file() else ""
    summary_ok = bool(re.search(r"4\.8(?:0|1)?", text) and ("increase" in text or "higher" in text or "up" in text or "升高" in text or "上调" in text))
    criteria.update(direction=direction_ok, report_consistent=summary_ok)
    if not tests["gene"]: failures.append("WRONG_GENE")
    if not tests["source_file"] or not tests["source_sheet"]: failures.append("UNTRACEABLE_SOURCE")
    if core < 40: failures.append("EFFECT_VALUES_MISMATCH")
    if not direction_ok: failures.append("DIRECTION_MISMATCH")
    if not summary_ok: failures.append("REPORT_MISMATCH")
    return {"core_science": core, "direction": 15 if direction_ok else 0, "summary": 5 if summary_ok else 0,
            "hardgate_pass": tests["gene"] and tests["source_file"] and tests["source_sheet"] and core == 40 and direction_ok,
            "failure_codes": failures, "criteria": criteria}
