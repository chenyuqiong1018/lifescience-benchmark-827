from __future__ import annotations

import json
import math
import re
from pathlib import Path

ACCEPTED = True  # BixBench range/exact gold; acceptance suite in scripts/acceptance-ls06-ls10.py.
READY_CANDIDATE = True


EXPECTED = {"KD_mean": 214.5, "CTRL_mean": 210.625, "KD_sd": 10.9414023651,
            "CTRL_sd": 22.8531804601, "pooled_sd": 17.9162236933, "d": 0.2162844172}


def close(v, e, tol=1e-5):
    try: return math.isclose(float(v), e, rel_tol=tol, abs_tol=tol)
    except (TypeError, ValueError): return False


def pair(obj, labels, key):
    value = obj.get(key, {})
    if isinstance(value, dict): return value.get(labels[0]), value.get(labels[1])
    if isinstance(value, list) and len(value) == 2: return value
    return None, None


def check(workspace: Path):
    failures, criteria = [], {}
    try: data = json.loads((workspace / "output" / "power_result.json").read_text(encoding="utf-8"))
    except Exception: data = {}
    labels = data.get("group_labels", [])
    valid_labels = isinstance(labels, list) and len(labels) == 2 and set(map(str.upper, labels)) == {"KD", "CTRL"}
    labels = labels if valid_labels else ["KD", "CTRL"]
    m1,m2=pair(data,labels,"means"); s1,s2=pair(data,labels,"sds")
    mapped_m={str(labels[0]).upper():m1,str(labels[1]).upper():m2}; mapped_s={str(labels[0]).upper():s1,str(labels[1]).upper():s2}
    tests={"labels":valid_labels,
           "means":close(mapped_m.get("KD"),EXPECTED["KD_mean"]) and close(mapped_m.get("CTRL"),EXPECTED["CTRL_mean"]),
           "sds":close(mapped_s.get("KD"),EXPECTED["KD_sd"]) and close(mapped_s.get("CTRL"),EXPECTED["CTRL_sd"]),
           "pooled_sd":close(data.get("pooled_sd"),EXPECTED["pooled_sd"]),
           "d":close(abs(float(data.get("cohens_d",math.nan))),EXPECTED["d"],tol=5e-3),
           "alpha":close(data.get("alpha"),.05), "power":close(data.get("power"),.8),
           "alternative":str(data.get("alternative","")).lower() in {"two-sided","two sided","two_sided","two.sided"},
           "required_n":data.get("required_n_per_group")==337}
    core=(8 if tests["means"] else 0)+(8 if tests["sds"] else 0)+(6 if tests["pooled_sd"] else 0)+(8 if tests["d"] else 0)+(10 if tests["required_n"] else 0)
    decision_ok=tests["alpha"] and tests["power"] and tests["alternative"] and tests["required_n"]
    report=workspace/"output"/"report.md"; text=report.read_text(encoding="utf-8",errors="replace").lower() if report.is_file() else ""
    summary_ok=bool(re.search(r"0\.21[56]",text) and re.search(r"\b337\b",text) and "0.05" in text and ("0.8" in text or "80%" in text))
    criteria.update(tests,report_consistent=summary_ok)
    if not tests["labels"]:failures.append("GROUP_LABELS_INVALID")
    if core<40:failures.append("POWER_CALCULATION_MISMATCH")
    if not decision_ok:failures.append("POWER_SPECIFICATION_MISMATCH")
    if not summary_ok:failures.append("REPORT_MISMATCH")
    return {"core_science":core,"direction":15 if decision_ok else 0,"summary":5 if summary_ok else 0,
            "hardgate_pass":tests["labels"] and core==40 and decision_ok,"failure_codes":failures,"criteria":criteria}
