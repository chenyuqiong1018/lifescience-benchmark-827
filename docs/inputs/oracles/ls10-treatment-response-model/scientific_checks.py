from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

ACCEPTED = True  # BixBench published ranges plus frozen coding; acceptance suite in scripts/acceptance-ls06-ls10.py.
READY_CANDIDATE = True


EXPECTED = {"estimate": -0.0795084690, "std_error": 0.0262977303,
            "p_value": 0.0024995441, "odds_ratio": 0.9235701982}


def close(v, e, *, rel=2e-3, abs_=1e-6):
    try:return math.isclose(float(v),e,rel_tol=rel,abs_tol=abs_)
    except (TypeError,ValueError):return False


def check(workspace: Path):
    failures,criteria=[],{}
    try:
        with (workspace/"output"/"model_coefficients.csv").open(encoding="utf-8-sig",newline="") as f: rows=list(csv.DictReader(f))
    except Exception:rows=[]
    ages=[r for r in rows if str(r.get("term","")).strip().lower()=="age"]
    age=ages[0] if len(ages)==1 else {}
    tests={k:close(age.get(k),v,rel=3e-3,abs_=5e-5) for k,v in EXPECTED.items()}
    try:meta=json.loads((workspace/"output"/"model_metadata.json").read_text(encoding="utf-8"))
    except Exception:meta={}
    meta_text=json.dumps(meta,ensure_ascii=False).lower()
    tests.update({"unique_age":len(ages)==1,"terms":all(x in meta_text for x in ("bmi","age","gender")),
                  "outcome":("pr" in meta_text and ("1" in meta_text or "positive" in meta_text)),
                  "gender_reference":("female" in meta_text),"complete_cases":("80" in meta_text or meta.get("n_complete_cases")==80)})
    core=sum(10 for k in ("estimate","std_error","p_value","odds_ratio") if tests[k])
    decision_ok=tests["estimate"] and tests["odds_ratio"] and tests["p_value"] and float(age.get("estimate",0))<0 and float(age.get("odds_ratio",2))<1 and float(age.get("p_value",1))<.05
    report=workspace/"output"/"report.md";text=report.read_text(encoding="utf-8",errors="replace").lower() if report.is_file() else ""
    summary_ok=bool(re.search(r"-0\.0(?:79|80)",text) and re.search(r"0\.002[45]",text) and ("decrease" in text or "lower" in text or "降低" in text or "下降" in text))
    criteria.update(tests,report_consistent=summary_ok)
    if not tests["unique_age"]:failures.append("AGE_TERM_NOT_UNIQUE")
    if not all(tests[k] for k in ("terms","outcome","gender_reference","complete_cases")):failures.append("MODEL_METADATA_INCOMPLETE")
    if core<40:failures.append("AGE_MODEL_VALUES_MISMATCH")
    if not decision_ok:failures.append("AGE_DIRECTION_MISMATCH")
    if not summary_ok:failures.append("REPORT_MISMATCH")
    return {"core_science":core,"direction":15 if decision_ok else 0,"summary":5 if summary_ok else 0,
            "hardgate_pass":tests["unique_age"] and all(tests[k] for k in ("terms","outcome","gender_reference","complete_cases")) and core==40 and decision_ok,
            "failure_codes":failures,"criteria":criteria}
