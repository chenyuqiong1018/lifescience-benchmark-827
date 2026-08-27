from __future__ import annotations

import importlib.util
import json
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("ls02_scientific_checks", HERE / "scientific_checks.py")
CHECKER = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(CHECKER)


def _fixture(files: dict[str, str]):
    root = tempfile.TemporaryDirectory(prefix="ls02-accept-")
    output = Path(root.name) / "output"
    output.mkdir()
    for name, content in files.items():
        path = output / name
        path.write_text(content, encoding="utf-8")
    return root, Path(root.name)


def _run(files: dict[str, str]):
    temp, workspace = _fixture(files)
    try:
        result = CHECKER.check(workspace)
        sentinel = workspace / "candidate-code-ran.txt"
        assert not sentinel.exists(), "candidate analysis.py was executed"
        assert set(result) == {"core_science", "direction", "summary", "hardgate_pass", "criteria", "failure_codes"}
        return result
    finally:
        temp.cleanup()


CORRECT = {
    "deletion.tsv": (
        "chrom\tstart_100kb\tend_100kb\tsize_bp\tsupporting_signals\n"
        "chr22\t20000000\t21000000\t1000000\tHeterozygous deletion with read-depth depletion to about half of both flanks\n"
    ),
    "qc.json": json.dumps({"deletion_to_flank_depth_ratio": 0.4923}),
    "report.md": (
        "# Result\nA heterozygous deletion spans chr22:20-21 Mb (1 Mb). Read coverage drops to about half of "
        "the flanks. Breakpoints are approximate and rounded to the nearest 100 kb; shallow data do not justify exact breakpoints.\n"
    ),
    # Deliberately executable-looking candidate code: the suite and checker must never run it.
    "analysis.py": "from pathlib import Path\nPath('../candidate-code-ran.txt').write_text('bad')\n",
}


def reference_like_correct():
    result = _run(CORRECT)
    assert (result["core_science"], result["direction"], result["summary"]) == (40, 15, 5)
    assert result["hardgate_pass"] is True
    return result


def empty_or_missing():
    result = _run({"analysis.py": CORRECT["analysis.py"]})
    assert (result["core_science"], result["direction"], result["summary"]) == (0, 0, 0)
    assert result["hardgate_pass"] is False
    return result


def scientifically_wrong():
    result = _run({
        "deletion.tsv": "chrom\tstart_100kb\tend_100kb\tsize_bp\tsupporting_signals\nchr22\t21000000\t22000000\t1000000\tCopy-number gain with normal depth\n",
        "qc.json": '{"coverage_ratio": 1.0}',
        "report.md": "A chr22:21-22 Mb duplication has exact base-pair breakpoint precision.\n",
    })
    assert result["core_science"] < 30 and result["direction"] == 0 and result["summary"] == 0
    assert result["hardgate_pass"] is False
    return result


def plausible_but_unsupported():
    wrong_locus = _run({
        "deletion.tsv": "chrom\tstart_100kb\tend_100kb\tsize_bp\tsupporting_signals\nchr22\t18000000\t19000000\t1000000\tDeletion; dramatic depth depletion\n",
        "qc.json": '{"depth_ratio": 0.49}',
        "report.md": "A deletion at chr22:18-19 Mb has reduced coverage and 100 kb rounded resolution.\n",
    })
    fabricated_number = _run({
        "deletion.tsv": "chrom\tstart_100kb\tend_100kb\tsize_bp\tsupporting_signals\nchr22\t20000000\t21000000\t1000000\tDeletion supported by a depth ratio of 0.05\n",
        "qc.json": '{"deletion_to_flank_depth_ratio": 0.05}',
        "report.md": "A deletion spans chr22:20-21 Mb (1 Mb); depth ratio is 0.05.\n",
    })
    negated_evidence = _run({
        "deletion.tsv": "chrom\tstart_100kb\tend_100kb\tsize_bp\tsupporting_signals\nchr22\t20000000\t21000000\t1000000\tDeletion, but no evidence of reduced depth\n",
        "qc.json": '{"depth_ratio": 0.49}',
        "report.md": "The chr22:20-21 Mb (1 Mb) interval has normal, unchanged coverage; no depth drop.\n",
    })
    for result in (wrong_locus, fabricated_number, negated_evidence):
        assert result["hardgate_pass"] is False
        assert result["criteria"]["fatal_gates"]["FATAL_GROUNDED_DEPTH_EVIDENCE"] is False
        assert result["criteria"]["deletion_direction"] is False
    assert wrong_locus["core_science"] <= 14
    assert fabricated_number["core_science"] == 30
    return {"wrong_locus": wrong_locus, "fabricated_number": fabricated_number, "negated_evidence": negated_evidence}


def valid_alternative_implementation():
    result = _run({
        "deletion.tsv": (
            "contig;breakpoint_end;support;length_bp;breakpoint_start\n"
            "22;20.95 Mb;CNV loss with local coverage reduced to roughly half of flanks;0.98 Mb;20.05 Mb\n"
        ),
        "qc.json": json.dumps({"metrics": {"coverage_to_flank_ratio": 0.51}}),
        "report.md": (
            "The CNV loss covers chr22:20.0 to 21.0 Mb (1.0 Mb). Local coverage is roughly half the flanks. "
            "Coordinates are approximate at 100-kb resolution because the shallow data limit breakpoint precision.\n"
        ),
    })
    assert (result["core_science"], result["direction"], result["summary"]) == (40, 15, 5)
    assert result["hardgate_pass"] is True
    return result


CASES = {
    "reference_like_correct": reference_like_correct,
    "empty_or_missing": empty_or_missing,
    "scientifically_wrong": scientifically_wrong,
    "plausible_but_unsupported": plausible_but_unsupported,
    "valid_alternative_implementation": valid_alternative_implementation,
}


def main() -> int:
    records = []
    for name, function in CASES.items():
        try:
            function()
            records.append({"case": name, "passed": True})
        except Exception as exc:
            records.append({"case": name, "passed": False, "error": f"{type(exc).__name__}: {exc}"})
    payload = {
        "schema_version": 1,
        "all_passed": all(record["passed"] for record in records),
        "cases": records,
        "required_case_count": 5,
        "adversarial_variants": ["wrong_locus_self_report", "fabricated_depth_ratio", "negated_depth_evidence"],
        "candidate_code_executed": False,
    }
    (HERE / "acceptance-result.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0 if payload["all_passed"] and len(records) == 5 else 1


if __name__ == "__main__":
    raise SystemExit(main())
