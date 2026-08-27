from __future__ import annotations

import csv
import json
import tempfile
from pathlib import Path

import scientific_checks as checker

HERE = Path(__file__).resolve().parent

LOG_PAIRS = [
    ("Bladder", "Limb", 0.931977802),
    ("Gill", "Lung", 0.928067276),
    ("Kidney", "Limb", 0.921551695),
    ("Bladder", "Kidney", 0.916855291),
    ("Brain", "Limb", 0.912854981),
    ("Intestine", "Pancreas", 0.905623437),
    ("Intestine", "Stomach", 0.892433037),
    ("Cloaca", "Intestine", 0.890614893),
    ("Liver", "Pancreas", 0.889239530),
    ("GallBladder", "Heart", 0.210824524),
    ("Heart", "Liver", 0.235150218),
    ("Heart", "Stomach", 0.244384278),
]

COSINE_PAIRS = [
    ("Bladder", "Limb", 0.996615028),
    ("Gill", "Lung", 0.984860447),
    ("Kidney", "Limb", 0.962781724),
    ("Brain", "Limb", 0.957717448),
    ("Intestine", "Pancreas", 0.928022497),
    ("Intestine", "Stomach", 0.996091179),
    ("Cloaca", "Intestine", 0.988922599),
    ("Liver", "Pancreas", 0.919487589),
    ("Bladder", "Prostate", 0.738367139),
]

RATIONALE = (
    "There is one profile per organ and no external identity reference or replicate. "
    "Pairwise correlation alone cannot identify a unique transposed label pair."
)
REPORT = (
    "# ATAC label assessment\n\nThe supplied data do not support a unique swap. "
    "There is a single profile per organ and no external reference or replicate, so "
    "correlation alone cannot identify which labels, if any, were transposed.\n"
)


def _base(root: Path, swap=False, organ_a=None, organ_b=None, evidence=RATIONALE, report=REPORT):
    out = root / "output"
    out.mkdir(parents=True)
    (out / "swap_call.json").write_text(json.dumps({
        "swap_detected": swap,
        "organ_a": organ_a,
        "organ_b": organ_b,
        "confidence": "low" if not swap else "high",
        "evidence": evidence,
    }), encoding="utf-8")
    (out / "report.md").write_text(report, encoding="utf-8")
    return out


def _long(out: Path, pairs):
    with (out / "sample_similarity.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["organ_a", "organ_b", "correlation"])
        writer.writerows(pairs)


def _square(out: Path, pairs):
    organs = sorted({x for a, b, _ in pairs for x in (a, b)})
    values = {(a, b): v for a, b, v in pairs}
    values.update({(b, a): v for a, b, v in pairs})
    with (out / "sample_similarity.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([""] + organs)
        for a in organs:
            writer.writerow([a] + [(1.0 if a == b else values.get((a, b), "")) for b in organs])


def _score_ok(result):
    return (
        result["core_science"] == 40
        and result["direction"] == 15
        and result["summary"] == 5
        and result["hardgate_pass"] is True
    )


def main():
    outcomes = {}
    with tempfile.TemporaryDirectory(prefix="ls03-science-acceptance-") as temp:
        root = Path(temp)

        case = root / "reference_like_correct"
        _long(_base(case), LOG_PAIRS)
        result = checker.check(case)
        outcomes[case.name] = {"passed": _score_ok(result), "result": result}

        case = root / "empty_or_missing"
        case.mkdir()
        result = checker.check(case)
        passed = (result["core_science"], result["direction"], result["summary"]) == (0, 0, 0) and not result["hardgate_pass"]
        outcomes[case.name] = {"passed": passed, "result": result}

        case = root / "scientifically_wrong"
        wrong_report = "# Call\n\nBladder and Limb are definitively swapped.\n"
        out = _base(case, True, "Bladder", "Limb", "Their correlation proves the transposition.", wrong_report)
        _long(out, LOG_PAIRS)
        result = checker.check(case)
        passed = not result["hardgate_pass"] and result["direction"] == 0 and "UNSUPPORTED_SWAP_DIRECTION" in result["failure_codes"]
        outcomes[case.name] = {"passed": passed, "result": result}

        case = root / "plausible_but_unsupported_self_report"
        out = _base(case)
        fabricated = [
            ("Bladder", "Limb", .71), ("Gill", "Lung", .67),
            ("Kidney", "Limb", .72), ("Brain", "Limb", .70),
            ("Intestine", "Pancreas", .64), ("Intestine", "Stomach", .65),
            ("Cloaca", "Intestine", .63), ("Liver", "Pancreas", .62),
            ("Bladder", "Prostate", .91),
        ]
        _long(out, fabricated)
        result = checker.check(case)
        passed = (
            not result["hardgate_pass"]
            and not result["criteria"]["similarity_numeric_agreement"]
            and "SIMILARITY_NOT_GROUNDED" in result["failure_codes"]
        )
        outcomes[case.name] = {"passed": passed, "result": result}

        case = root / "valid_alternative_method_layout"
        out = _base(case)
        _square(out, COSINE_PAIRS)
        result = checker.check(case)
        passed = _score_ok(result) and result["criteria"]["matched_method"] == "cosine"
        outcomes[case.name] = {"passed": passed, "result": result}

    all_passed = all(item["passed"] for item in outcomes.values())
    payload = {
        "schema_version": 1,
        "accepted": all_passed,
        "all_cases_passed": all_passed,
        "adversarial_false_positive_passed": outcomes["plausible_but_unsupported_self_report"]["passed"],
        "candidate_code_executed": False,
        "cases": outcomes,
    }
    (HERE / "acceptance-result.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    if not all_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
