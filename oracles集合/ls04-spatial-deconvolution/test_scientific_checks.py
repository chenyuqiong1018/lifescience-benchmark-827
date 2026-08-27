from __future__ import annotations

import importlib.util
import json
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ls04_scientific_checks", HERE / "scientific_checks.py")
checker = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(checker)


def _fixture(composition: str | None, report: str = ""):
    temp = tempfile.TemporaryDirectory()
    root = Path(temp.name)
    output = root / "output"
    output.mkdir()
    sentinel = root / "candidate_code_executed"
    # Deliberately dangerous candidate code: the acceptance suite proves the
    # checker only reads artifacts and never imports/runs analysis.py.
    (output / "analysis.py").write_text(
        f"from pathlib import Path\nPath({str(sentinel)!r}).write_text('EXECUTED')\nraise RuntimeError('must never run')\n",
        encoding="utf-8",
    )
    if composition is not None:
        (output / "spot_710_composition.csv").write_text(composition, encoding="utf-8")
    if report:
        (output / "report.md").write_text(report, encoding="utf-8")
    result = checker.check(root)
    never_executed = not sentinel.exists()
    return temp, result, never_executed


def _run_case(name: str, variants: list[tuple[str, str | None, str, object]]):
    details = []
    passed = True
    for variant, composition, report, predicate in variants:
        temp, result, never_executed = _fixture(composition, report)
        try:
            ok = bool(predicate(result)) and never_executed
            details.append({
                "variant": variant,
                "passed": ok,
                "scores": {k: result[k] for k in ("core_science", "direction", "summary")},
                "hardgate_pass": result["hardgate_pass"],
                "candidate_code_not_executed": never_executed,
            })
            passed &= ok
        finally:
            temp.cleanup()
    return {"case": name, "passed": passed, "variants": details}


CORRECT = """cell_type,weight,evidence
B_Cell,0.323040493130,reference-profile NNLS
Endothelial,0.346842122837,reference-profile NNLS
Fibroblast_Stroma,0.0,reference-profile NNLS
Macrophage,0.318705728025,reference-profile NNLS
T_Cell,0.0114116560073,reference-profile NNLS
Tumor_Core,0.0,reference-profile NNLS
"""

WRONG = """cell_type,weight,evidence
Tumor_Core,1.0,candidate assertion only
B_Cell,0.0,candidate assertion only
Endothelial,0.0,candidate assertion only
Macrophage,0.0,candidate assertion only
"""

ALTERNATIVE = """evidence;fraction;type
robust regression;0.34;B lymphocyte
robust regression;0.33;Endothelial cells
robust regression;0.32;Macrophages
robust regression;0.01;T-cell
"""


def main():
    cases = [
        _run_case("reference_like_correct", [
            ("canonical", CORRECT,
             "Spot_710-1 is a mixed composition with contributions from B cells, endothelial cells, and macrophages; T cells are negligible and tumor cells are absent.",
             lambda r: r["hardgate_pass"] and r["core_science"] == 40 and r["direction"] == 15 and r["summary"] == 5),
        ]),
        _run_case("empty_or_missing", [
            ("missing_composition", None, "No result was available.",
             lambda r: not r["hardgate_pass"] and r["core_science"] == 0 and r["direction"] == 0 and r["summary"] == 0),
            ("empty_composition", "", "B cells, endothelial cells, and macrophages form a mixture.",
             lambda r: not r["hardgate_pass"] and r["summary"] == 0),
        ]),
        _run_case("scientifically_wrong", [
            ("unsupported_tumor_single_type", WRONG, "Tumor cells dominate Spot_710-1 as a single cell type.",
             lambda r: not r["hardgate_pass"] and r["core_science"] <= 10 and r["direction"] == 0 and r["summary"] == 0),
            ("negated_supported_evidence", WRONG,
             "B cells are not present, endothelial cells are absent, and macrophages are unsupported; tumor cells dominate.",
             lambda r: not r["hardgate_pass"] and r["summary"] == 0),
        ]),
        _run_case("plausible_but_unsupported", [
            ("correct_labels_only_wrong_numbers", WRONG,
             "The conclusion is a B-cell, endothelial-cell, and macrophage mixture at Spot_710-1.",
             lambda r: not r["hardgate_pass"] and r["summary"] == 0),
            ("fabricated_reference_like_numbers", WRONG,
             "Spot_710-1 is a mixture: B cells 32.3%, endothelial cells 34.7%, and macrophages 31.9%, despite the submitted tumor-only composition.",
             lambda r: not r["hardgate_pass"] and r["summary"] == 0),
            ("fabricated_report_numbers_despite_valid_composition", CORRECT,
             "Spot_710-1 is a mixture of B cells 90%, endothelial cells 5%, and macrophages 5%.",
             lambda r: r["hardgate_pass"] and r["summary"] < 5 and not r["criteria"]["report_numbers_consistent_with_submission"]),
            ("directly_negated_evidence_despite_valid_composition", CORRECT,
             "This is a mixture, but B cells are not present; endothelial cells and macrophages contribute.",
             lambda r: r["hardgate_pass"] and r["summary"] < 5 and not r["criteria"]["report_expected_types_affirmed"]),
        ]),
        _run_case("valid_alternative_implementation", [
            ("aliases_semicolon_and_unrelated_negation", ALTERNATIVE,
             "B cells, endothelial cells, and macrophages—not tumor cells—form a mixture at Spot_710-1.",
             lambda r: r["hardgate_pass"] and r["core_science"] >= 37 and r["direction"] == 15 and r["summary"] == 5),
        ]),
    ]
    payload = {
        "schema_version": 1,
        "checker_accepted": checker.ACCEPTED is True,
        "all_passed": all(case["passed"] for case in cases),
        "required_case_count": 5,
        "cases": cases,
        "candidate_code_executed": False,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    raise SystemExit(0 if payload["all_passed"] else 1)


if __name__ == "__main__":
    main()
