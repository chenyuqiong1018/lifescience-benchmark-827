from __future__ import annotations

import csv
import re
from pathlib import Path

ACCEPTED = True

EXPECTED = {
    "c01": (True, True, True, True),
    "c02": (False, True, False, True),
    "c03": (True, True, True, True),
}
FIELDS = ("frame_ok", "start_ok", "stop_ok", "tag_ok")


def _truth(value):
    value = str(value or "").strip().lower()
    if value in {"true", "t", "yes", "y", "1", "ok", "pass", "compatible", "present", "intact"}:
        return True
    if value in {"false", "f", "no", "n", "0", "fail", "incompatible", "absent", "missing"}:
        return False
    return None


def _rows(path: Path):
    try:
        with path.open(encoding="utf-8-sig", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle)
            rows = {}
            for raw in reader:
                row = {str(k or "").strip().lower(): v for k, v in raw.items()}
                cid = str(row.get("construct_id", "")).strip().lower()
                if cid and cid not in rows:
                    rows[cid] = row
            return rows
    except Exception:
        return {}


def _has(text, pattern):
    return bool(re.search(pattern, str(text or "").lower()))


def _context(text, cid):
    match = re.search(rf"\b{cid}\b(.{{0,360}}?)(?=\bc0[123]\b|$)", text, re.S)
    return (match.group(0) if match else "").lower()


def check(workspace: Path):
    rows = _rows(Path(workspace) / "output" / "construct_audit.csv")
    criteria = {}
    correct = 0
    for cid, expected in EXPECTED.items():
        row = rows.get(cid, {})
        for field, wanted in zip(FIELDS, expected):
            ok = _truth(row.get(field)) is wanted
            criteria[f"{cid}.{field}"] = ok
            correct += int(ok)
    core = round(40 * correct / 12)

    positive = r"\b(pass(?:ed|es)?|ok|valid|acceptable|compliant|clear|good)\b"
    flagged = r"\b(fail(?:ed)?|error|invalid|reject(?:ed)?|non.?compliant|review|warn(?:ing)?|flagged|attention|inconsistent|issue)\b"
    status_checks = {
        "c01.status_direction": _has(rows.get("c01", {}).get("overall_status"), positive),
        "c02.status_direction": _has(rows.get("c02", {}).get("overall_status"), flagged),
        "c03.status_direction": _has(rows.get("c03", {}).get("overall_status"), flagged),
    }
    c01_issues = str(rows.get("c01", {}).get("issues", "")).strip().lower()
    issue_checks = {
        "c01.issues_clean": "c01" in rows and (not c01_issues or bool(re.fullmatch(r"(?:none|n/?a|no issues?|clear|ok|pass)", c01_issues))),
        "c02.issues_frame_and_stop": _has(rows.get("c02", {}).get("issues"), r"frame|triplet|divisib|modulo|multiple\s+of\s+3|length")
        and _has(rows.get("c02", {}).get("issues"), r"stop|terminat|end\s+codon"),
        "c03.issues_claim_conflict": _has(rows.get("c03", {}).get("issues"), r"claim|declar|metadata|annotat|report|out.?of.?frame|disagree|mismatch|inconsisten|conflict|contrar"),
    }
    criteria.update(status_checks, **issue_checks)
    direction = 3 * sum(status_checks.values()) + 2 * sum(issue_checks.values())

    report_path = Path(workspace) / "output" / "report.md"
    try:
        report = report_path.read_text(encoding="utf-8", errors="replace").lower()
    except Exception:
        report = ""
    c01, c02, c03 = (_context(report, cid) for cid in EXPECTED)
    report_checks = {
        "report_all_constructs": all(re.search(rf"\b{cid}\b", report) for cid in EXPECTED),
        "report_c01_pass": _has(c01, positive),
        "report_c02_frame": _has(c02, r"frame|triplet|divisib|multiple\s+of\s+3|length"),
        "report_c02_stop": _has(c02, r"no\s+(?:terminal\s+)?stop|missing\s+(?:a\s+)?(?:terminal\s+)?stop|lacks?\s+(?:a\s+)?(?:terminal\s+)?stop|stop.+(?:absent|missing|false|fail)"),
        "report_c03_claim_conflict": _has(c03, r"claim|declar|metadata|annotat|out.?of.?frame|disagree|mismatch|inconsisten|conflict|contrar"),
    }
    criteria.update(report_checks)
    summary = sum(report_checks.values())

    failures = []
    if core < 40:
        failures.append("SCIENTIFIC_AUDIT_MISMATCH")
    if direction < 15:
        failures.append("AUDIT_DIRECTION_MISMATCH")
    if summary < 5:
        failures.append("REPORT_MISMATCH")
    return {
        "core_science": core,
        "direction": direction,
        "summary": summary,
        "hardgate_pass": core == 40 and direction == 15,
        "criteria": criteria,
        "failure_codes": failures,
    }
