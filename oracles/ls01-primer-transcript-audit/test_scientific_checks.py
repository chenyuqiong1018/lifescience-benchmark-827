from __future__ import annotations

import importlib.util
import json
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
CHECKER_PATH = HERE / "scientific_checks.py"
RESULT_PATH = HERE / "acceptance-result.json"
EXPECTED_RESULT_KEYS = {
    "core_science", "direction", "summary", "hardgate_pass", "criteria", "failure_codes"
}


def _load_checker():
    spec = importlib.util.spec_from_file_location("ls01_scientific_checks", CHECKER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("unable to load scientific checker")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _workspace(csv_text: str | None, report_text: str | None):
    temporary = tempfile.TemporaryDirectory(prefix="ls01-acceptance-")
    root = Path(temporary.name)
    output = root / "output"
    output.mkdir()
    if csv_text is not None:
        (output / "primer_audit.csv").write_text(csv_text, encoding="utf-8")
    if report_text is not None:
        (output / "report.md").write_text(report_text, encoding="utf-8")
    # A deliberately explosive candidate script proves the suite never imports or executes it.
    (output / "analysis.py").write_text(
        "raise RuntimeError('candidate code must never execute')\n", encoding="utf-8"
    )
    return temporary, root


CORRECT_CSV = """pair_id,transcripts_matched,amplicon_length,cds_compatible,status,reason
p01,TX_CANONICAL,102,unknown,FAIL,Observed amplicon is 102 bp rather than the claimed 108 bp; TX_CANONICAL CDS metadata are invalid because coordinate 700 exceeds sequence length.
p02,TX_ALT,99,unknown,FAIL,Observed amplicon is 99 bp rather than the claimed 104 bp; TX_ALT CDS metadata are invalid because coordinate 640 exceeds sequence length.
p03,none,NA,not applicable,FAIL,No amplicon matches any supplied transcript; CDS compatibility cannot be assessed.
"""

CORRECT_REPORT = """# Primer audit

- p01 matches only TX_CANONICAL and produces a 102 bp amplicon, not the claimed 108 bp.
- p02 matches only TX_ALT and produces a 99 bp amplicon, not the claimed 104 bp.
- p03 has no amplicon in either supplied transcript.

The CDS metadata are invalid for both transcripts: TX_CANONICAL ends at coordinate 700 and
TX_ALT at 640, but each coordinate exceeds the length of its supplied transcript sequence.
"""

WRONG_CSV = """pair_id,transcripts_matched,amplicon_length,cds_compatible,status,reason
p01,TX_ALT,108,true,PASS,Matches the expected product.
p02,TX_CANONICAL,104,true,PASS,Matches the expected product.
p03,TX_CANONICAL,120,true,PASS,Matches the expected product.
"""

WRONG_REPORT = """p01 is 108 bp, p02 is 104 bp, and p03 is 120 bp. CDS=101-700 and CDS=101-640."""

ALTERNATIVE_CSV = """primer pair;product bp;matched transcripts;CDS compatibility;assessment;explanation
p01;102;TX_CANONICAL;undetermined;requires attention;The product is 102 bp rather than expected 108 bp
p02;99;TX_ALT;cannot assess;discrepant;The observed size differs from the claimed size of 104 bp
p03;n/a;none;not applicable;issue;No product was amplified from any supplied transcript
"""

ALTERNATIVE_REPORT = """An independent exact-match and reverse-complement search found p01 only on
TX_CANONICAL (102 base pairs), while p02 mapped only to TX_ALT with a 99-bp product.
p03 yielded no amplicon. For each supplied transcript, the CDS annotation is inconsistent:
the annotated end coordinate lies beyond the actual transcript sequence length.
"""


def _full_pass(result: dict) -> bool:
    return (
        set(result) == EXPECTED_RESULT_KEYS
        and result.get("core_science") == 40
        and result.get("direction") == 15
        and result.get("summary") == 5
        and result.get("hardgate_pass") is True
        and result.get("failure_codes") == []
    )


def run() -> dict:
    checker = _load_checker()
    cases: dict[str, bool] = {}
    details: dict[str, object] = {}

    temporary, root = _workspace(CORRECT_CSV, CORRECT_REPORT)
    try:
        result = checker.check(root)
        cases["reference_like_correct"] = _full_pass(result)
        details["reference_like_correct"] = result
    finally:
        temporary.cleanup()

    missing_temp, missing_root = _workspace(None, None)
    empty_temp, empty_root = _workspace("", "")
    try:
        missing = checker.check(missing_root)
        empty = checker.check(empty_root)
        cases["empty_or_missing"] = all(
            set(result) == EXPECTED_RESULT_KEYS
            and result.get("hardgate_pass") is False
            and result.get("core_science") == 0
            and result.get("direction") == 0
            and result.get("summary") == 0
            for result in (missing, empty)
        )
        details["empty_or_missing"] = {"missing": missing, "empty": empty}
    finally:
        missing_temp.cleanup()
        empty_temp.cleanup()

    temporary, root = _workspace(WRONG_CSV, WRONG_REPORT)
    try:
        result = checker.check(root)
        cases["scientifically_wrong"] = (
            set(result) == EXPECTED_RESULT_KEYS
            and result.get("hardgate_pass") is False
            and result.get("core_science", 40) < 40
            and result.get("direction", 15) < 15
            and result.get("summary", 5) < 5
            and result.get("criteria", {}).get("invalid_cds_metadata_identified") is False
        )
        details["scientifically_wrong"] = result
    finally:
        temporary.cleanup()

    temporary, root = _workspace(ALTERNATIVE_CSV, ALTERNATIVE_REPORT)
    try:
        result = checker.check(root)
        cases["valid_alternative_implementation"] = _full_pass(result)
        details["valid_alternative_implementation"] = result
    finally:
        temporary.cleanup()

    payload = {
        "accepted": all(cases.values()),
        "cases": cases,
        "candidate_code_executed": False,
        "details": details,
    }
    RESULT_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


if __name__ == "__main__":
    result = run()
    print(json.dumps({"accepted": result["accepted"], "cases": result["cases"]}, indent=2))
    raise SystemExit(0 if result["accepted"] else 1)
