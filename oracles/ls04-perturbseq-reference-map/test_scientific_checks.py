from __future__ import annotations

import json
import tempfile
from pathlib import Path

from scientific_checks import check

HERE = Path(__file__).resolve().parent
RESULT_PATH = HERE / "acceptance-result.json"

CORRECT_CSV = """target_gene,query_guide_id,score,runner_up_score,confidence
PABPC1,guide18,0.5216191513,0.2819277,high
NUDT21,guide13,0.6313580675,0.4984148,medium
LEO1,guide14,0.3784984638,0.3137100,low
"""
SAFE_ANALYSIS = """from pathlib import Path

def align_profiles(query_adata, reference_adata):
    query_groups = query_adata.obs[\"guide\"]
    reference_labels = reference_adata.obs[\"target_gene\"]
    return query_groups, reference_labels
"""
POSITIVE_REPORT = """# Reference-map result

PABPC1 maps to guide18.
NUDT21 maps to guide13.
LEO1 maps to guide14.
LEO1 has the smallest relative top-versus-runner-up margin and therefore the lowest confidence.
"""


def fixture(files: dict[str, str]):
    root = tempfile.TemporaryDirectory(prefix="ls04-science-")
    workspace = Path(root.name)
    out = workspace / "output"
    out.mkdir()
    for name, content in files.items():
        (out / name).write_text(content, encoding="utf-8")
    return root, workspace


def run_case(files: dict[str, str]):
    temp, workspace = fixture(files)
    try:
        return check(workspace)
    finally:
        temp.cleanup()


def main() -> int:
    details: dict[str, dict] = {}

    correct = run_case({"guide_mapping.csv": CORRECT_CSV, "analysis.py": SAFE_ANALYSIS, "report.md": POSITIVE_REPORT})
    ok_correct = correct == {
        "core_science": 40,
        "direction": 15,
        "summary": 5,
        "hardgate_pass": True,
        "criteria": correct["criteria"],
        "failure_codes": [],
    }
    details["reference_like_correct"] = {"passed": ok_correct, "observed": correct}

    empty = run_case({})
    ok_empty = (empty["core_science"], empty["direction"], empty["summary"], empty["hardgate_pass"]) == (0, 0, 0, False)
    details["empty_or_missing"] = {"passed": ok_empty, "observed": empty}

    wrong_csv = """target_gene,query_guide_id,score,runner_up_score,confidence
PABPC1,guide13,0.52,0.28,high
NUDT21,guide14,0.63,0.50,medium
LEO1,guide18,0.38,0.31,low
"""
    wrong = run_case({"guide_mapping.csv": wrong_csv, "analysis.py": SAFE_ANALYSIS, "report.md": POSITIVE_REPORT})
    leakage = run_case({
        "guide_mapping.csv": CORRECT_CSV,
        "analysis.py": "query = read_h5ad(query_path)\nlabels = query.obs[\"target_gene\"]\n",
        "report.md": POSITIVE_REPORT,
    })
    ok_wrong = (
        wrong["core_science"] == 0 and not wrong["hardgate_pass"]
        and "FATAL_MAPPING_IDENTITY" in wrong["failure_codes"]
        and not leakage["hardgate_pass"]
        and "FATAL_NO_QUERY_TARGET_LEAKAGE" in leakage["failure_codes"]
    )
    details["scientifically_wrong"] = {"passed": ok_wrong, "observed": wrong, "target_leakage_variant": leakage}

    unsupported_csv = """target_gene,query_guide_id,score,runner_up_score,confidence
PABPC1,guide99,0.99,0.98,high
NUDT21,guide98,0.98,0.97,medium
LEO1,guide97,0.97,0.96,low
"""
    unsupported = run_case({"guide_mapping.csv": unsupported_csv, "analysis.py": "# self-report only\n", "report.md": POSITIVE_REPORT})
    fabricated_csv = """target_gene,query_guide_id,score,runner_up_score,confidence
PABPC1,guide18,0.99,0.01,high
NUDT21,guide13,0.98,0.01,medium
LEO1,guide14,0.97,0.01,low
"""
    fabricated = run_case({"guide_mapping.csv": fabricated_csv, "analysis.py": SAFE_ANALYSIS, "report.md": POSITIVE_REPORT})
    negated_report = """PABPC1 does not map to guide18.
NUDT21 does not map to guide13.
LEO1 does not map to guide14.
LEO1 does not have the smallest relative margin and is not the lowest-confidence match.
"""
    negated = run_case({"guide_mapping.csv": CORRECT_CSV, "analysis.py": SAFE_ANALYSIS, "report.md": negated_report})
    ok_unsupported = (
        unsupported["core_science"] == 0 and unsupported["summary"] == 0 and not unsupported["hardgate_pass"]
        and fabricated["criteria"]["numeric_scores_grounded"] == {"PABPC1": False, "NUDT21": False, "LEO1": False}
        and fabricated["direction"] < 15
        and negated["summary"] == 0
        and not any(negated["criteria"]["report_mapping_affirmed"].values())
        and not negated["criteria"]["report_ambiguity_affirmed"]
    )
    details["plausible_but_unsupported"] = {
        "passed": ok_unsupported,
        "observed": unsupported,
        "fabricated_number_variant": fabricated,
        "negated_evidence_variant": negated,
    }

    alternative_csv = """confidence,runner_up_score,extra_note,TARGET GENE,score,QUERY GUIDE ID
0.20,0.34,global assignment,leo1,0.40,GUIDE14
0.90,0.30,global assignment,pabpc1,0.55,Guide18
0.60,0.51,global assignment,nudt21,0.66,guide13
"""
    alternative_report = """| target | matched query group |
|---|---|
| LEO1 | guide14 |
| PABPC1 | guide18 |
| NUDT21 | guide13 |

Across the alternative global assignment, LEO1 has the lowest confidence because its relative runner-up margin is smallest.
"""
    alternative = run_case({"guide_mapping.csv": alternative_csv, "analysis.py": SAFE_ANALYSIS, "report.md": alternative_report})
    ok_alternative = (alternative["core_science"], alternative["direction"], alternative["summary"], alternative["hardgate_pass"]) == (40, 15, 5, True)
    details["valid_alternative_implementation"] = {"passed": ok_alternative, "observed": alternative}

    ordered = [
        "reference_like_correct", "empty_or_missing", "scientifically_wrong",
        "plausible_but_unsupported", "valid_alternative_implementation",
    ]
    all_passed = all(details[name]["passed"] for name in ordered)
    result = {
        "schema_version": 1,
        "task_id": "ls04-perturbseq-reference-map",
        "all_passed": all_passed,
        "cases_passed": sum(bool(details[name]["passed"]) for name in ordered),
        "cases_total": len(ordered),
        "candidate_code_executed": False,
        "cases": details,
    }
    RESULT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"all_passed": all_passed, "cases": {k: details[k]["passed"] for k in ordered}}, sort_keys=True))
    return 0 if all_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
