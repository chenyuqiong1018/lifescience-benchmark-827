from __future__ import annotations

import json
import math
import re
from pathlib import Path

ACCEPTED = True  # Benchmark gold independently encoded; acceptance suite in scripts/acceptance-ls06-ls10.py.
READY_CANDIDATE = True


def check(workspace: Path):
    failures, criteria = [], {}
    try:
        data = json.loads((workspace / "output" / "eno1_significance.json").read_text(encoding="utf-8"))
    except Exception:
        data = {}
    try: p = float(data.get("adjusted_p_value"))
    except (TypeError, ValueError): p = math.nan
    try: threshold = float(data.get("fdr_threshold"))
    except (TypeError, ValueError): threshold = math.nan
    tests = {
        "gene": str(data.get("gene", "")).upper() == "ENO1",
        "adjusted_p": math.isfinite(p) and math.isclose(p, 0.226, abs_tol=0.0005),
        "threshold": math.isfinite(threshold) and math.isclose(threshold, 0.05, abs_tol=1e-9),
        "decision": data.get("significant") is False,
        "source_file": "proteomic" in str(data.get("source_file", "")).lower(),
        "source_sheet": str(data.get("source_sheet", "")).strip().lower() == "tumor vs normal",
    }
    report = workspace / "output" / "report.md"
    text = report.read_text(encoding="utf-8", errors="replace").lower() if report.is_file() else ""
    summary_ok = bool(re.search(r"0\.226\b", text) and ("not significant" in text or "non-significant" in text or "不显著" in text))
    criteria.update(tests, report_consistent=summary_ok)
    core = 40 if tests["adjusted_p"] else 0
    decision_ok = tests["threshold"] and tests["decision"]
    if not tests["gene"]: failures.append("WRONG_GENE")
    if not tests["source_file"] or not tests["source_sheet"]: failures.append("UNTRACEABLE_SOURCE")
    if not tests["adjusted_p"]: failures.append("ADJUSTED_P_MISMATCH")
    if not decision_ok: failures.append("FDR_DECISION_MISMATCH")
    if not summary_ok: failures.append("REPORT_MISMATCH")
    return {"core_science": core, "direction": 15 if decision_ok else 0, "summary": 5 if summary_ok else 0,
            "hardgate_pass": tests["gene"] and tests["source_file"] and tests["source_sheet"] and core == 40 and decision_ok,
            "failure_codes": failures, "criteria": criteria}
