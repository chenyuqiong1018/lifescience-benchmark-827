from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("ls01_scientific_checks", HERE / "scientific_checks.py")
CHECKER = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(CHECKER)


def fixture(rows=None, report=None, headers=None):
    root = Path(tempfile.mkdtemp(prefix="ls01-acceptance-"))
    output = root / "output"
    output.mkdir()
    sentinel = root / "candidate-code-ran"
    (output / "analysis.py").write_text(
        "from pathlib import Path\n"
        f"Path({str(sentinel)!r}).write_text('ran')\n"
        "raise RuntimeError('candidate code must not execute')\n",
        encoding="utf-8",
    )
    if rows is not None:
        headers = headers or ["construct_id", "frame_ok", "start_ok", "stop_ok", "tag_ok", "overall_status", "issues"]
        with (output / "construct_audit.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=headers)
            writer.writeheader()
            writer.writerows(rows)
    if report is not None:
        (output / "report.md").write_text(report, encoding="utf-8")
    result = CHECKER.check(root)
    assert not sentinel.exists(), "checker executed candidate code"
    return result


def correct_case():
    rows = [
        dict(construct_id="c01", frame_ok="true", start_ok="true", stop_ok="true", tag_ok="true", overall_status="PASS", issues="none"),
        dict(construct_id="c02", frame_ok="false", start_ok="true", stop_ok="false", tag_ok="true", overall_status="FAIL", issues="length is not divisible by 3; missing terminal stop codon; claimed in_frame conflicts"),
        dict(construct_id="c03", frame_ok="true", start_ok="true", stop_ok="true", tag_ok="true", overall_status="REVIEW", issues="sequence is in frame, contrary to the claimed out_of_frame annotation"),
    ]
    report = """# ORF audit
c01 passes all represented checks.
c02 has a length not divisible by 3 and is missing a terminal stop codon.
c03 is in frame; this conflicts with the claimed out_of_frame annotation.
"""
    result = fixture(rows, report)
    assert result["core_science"] == 40 and result["direction"] == 15 and result["summary"] == 5
    assert result["hardgate_pass"] and not result["failure_codes"]


def empty_or_missing_case():
    missing = fixture()
    empty = fixture([], "")
    for result in (missing, empty):
        assert result["core_science"] == 0 and result["direction"] == 0 and result["summary"] == 0
        assert not result["hardgate_pass"]


def wrong_case():
    rows = [
        dict(construct_id=cid, frame_ok="false", start_ok="false", stop_ok="false", tag_ok="false", overall_status="PASS", issues="none")
        for cid in ("c01", "c02", "c03")
    ]
    result = fixture(rows, "c01 pass. c02 pass. c03 pass.")
    assert result["core_science"] < 20 and result["direction"] < 15 and result["summary"] < 5
    assert not result["hardgate_pass"] and "SCIENTIFIC_AUDIT_MISMATCH" in result["failure_codes"]


def alternative_case():
    headers = ["ISSUES", "TAG_OK", "OVERALL_STATUS", "STOP_OK", "CONSTRUCT_ID", "START_OK", "FRAME_OK"]
    rows = [
        {"CONSTRUCT_ID": "C03", "FRAME_OK": "1", "START_OK": "yes", "STOP_OK": "present", "TAG_OK": "compatible", "OVERALL_STATUS": "flagged", "ISSUES": "metadata says out-of-frame although the triplets are complete"},
        {"CONSTRUCT_ID": "C01", "FRAME_OK": "ok", "START_OK": "present", "STOP_OK": "yes", "TAG_OK": "intact", "OVERALL_STATUS": "valid", "ISSUES": "N/A"},
        {"CONSTRUCT_ID": "C02", "FRAME_OK": "0", "START_OK": "1", "STOP_OK": "missing", "TAG_OK": "intact", "OVERALL_STATUS": "warning", "ISSUES": "ORF length has incomplete triplet geometry and termination codon is absent"},
    ]
    report = """## Independent sequence review
C03 — Triplet-complete and actually in frame, inconsistent with its out-of-frame metadata.
C01 — Valid across the represented sequence features.
C02 — Length leaves an incomplete triplet; the termination stop codon is absent.
"""
    result = fixture(rows, report, headers)
    assert result["core_science"] == 40 and result["direction"] == 15 and result["summary"] == 5
    assert result["hardgate_pass"]


def main():
    cases = {
        "reference_like_correct": correct_case,
        "empty_or_missing": empty_or_missing_case,
        "scientifically_wrong": wrong_case,
        "valid_alternative_implementation": alternative_case,
    }
    results = {}
    for name, test in cases.items():
        try:
            test()
            results[name] = {"passed": True}
        except Exception as exc:
            results[name] = {"passed": False, "error": f"{type(exc).__name__}: {exc}"}
    payload = {"all_passed": all(item["passed"] for item in results.values()), "cases": results}
    (HERE / "acceptance-result.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if not payload["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
