from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULT_PATH = HERE / "acceptance-result.json"


def _load_checker():
    spec = importlib.util.spec_from_file_location("ls01_scientific_checks", HERE / "scientific_checks.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load scientific checker")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BASE = {
    "g01": (0.82, "high", "avoid", "0.82 activity; 1 mismatch coding exon hit and one coding bystander."),
    "g02": (0.67, "low", "advance", "0.67 activity; 3 mismatches in an intron give low annotated risk."),
    "g03": (0.74, "medium", "conditional review", "0.74 activity; 2 mismatches in a coding exon create moderate risk."),
    "g04": (0.59, "low", "shortlist", "0.59 activity; 4 mismatches at an intergenic site give low risk."),
    "g05": (0.78, "high", "do not advance", "0.78 activity; 1 mismatch coding exon hit and two coding bystanders."),
    "g06": (0.64, "low", "advance", "0.64 activity; 3 mismatches at an intergenic site give low risk."),
}


def _write_fixture(root: Path, order: list[str], *, risks=None, decisions=None, rationales=None, report_top=None):
    output = root / "output"
    output.mkdir(parents=True)
    # If a checker accidentally executes candidate code, this marker would appear.
    (output / "analysis.py").write_text(
        "from pathlib import Path\nPath(__file__).with_name('EXECUTED').write_text('bad')\n",
        encoding="utf-8",
    )
    fields = ["decision", "guide_id", "rationale", "rank", "risk_class", "on_target_score"]
    with (output / "ranked_guides.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for rank, guide_id in enumerate(order, 1):
            activity, risk, decision, rationale = BASE[guide_id]
            writer.writerow({
                "rank": rank,
                "guide_id": guide_id,
                "on_target_score": activity,
                "risk_class": (risks or {}).get(guide_id, risk),
                "decision": (decisions or {}).get(guide_id, decision),
                "rationale": (rationales or {}).get(guide_id, rationale),
            })
    top = report_top or order[0]
    (output / "report.md").write_text(
        f"# Result\n\nTop-ranked guide: {top}. The ranking balances activity against annotated "
        "off-target risk. Coding/exonic near matches are penalized: the unsafe candidates have "
        "1 mismatch, while intronic and intergenic sites with 3 or 4 mismatches are safer.\n",
        encoding="utf-8",
    )
    return output


def _run_case(checker, name, builder, predicate):
    with tempfile.TemporaryDirectory(prefix=f"ls01-{name}-") as temp_dir:
        root = Path(temp_dir)
        builder(root)
        result = checker.check(root)
        marker_absent = not (root / "output" / "EXECUTED").exists()
        passed = bool(predicate(result) and marker_absent and set(result) == {
            "core_science", "direction", "summary", "hardgate_pass", "criteria", "failure_codes"
        })
        return {
            "passed": passed,
            "scores": {key: result[key] for key in ("core_science", "direction", "summary")},
            "hardgate_pass": result["hardgate_pass"],
            "failure_codes": result["failure_codes"],
            "candidate_code_not_executed": marker_absent,
        }


def main():
    checker = _load_checker()
    cases = {}

    cases["reference_like_correct"] = _run_case(
        checker,
        "correct",
        lambda root: _write_fixture(root, ["g02", "g06", "g04", "g03", "g01", "g05"]),
        lambda r: r["core_science"] == 40 and r["direction"] == 15 and r["summary"] == 5 and r["hardgate_pass"],
    )

    def empty(root):
        (root / "output").mkdir(parents=True)
        (root / "output" / "analysis.py").write_text(
            "from pathlib import Path\nPath(__file__).with_name('EXECUTED').write_text('bad')\n",
            encoding="utf-8",
        )

    cases["empty_or_missing"] = _run_case(
        checker,
        "empty",
        empty,
        lambda r: r["core_science"] == 0 and r["direction"] == 0 and r["summary"] == 0 and not r["hardgate_pass"],
    )

    def wrong(root):
        rationales = {guide_id: "High activity, no safety concern." for guide_id in BASE}
        decisions = {"g01": "advance", "g05": "primary recommendation", "g03": "advance"}
        risks = {"g01": "low", "g02": "high", "g03": "low", "g04": "high", "g05": "low", "g06": "high"}
        output = _write_fixture(
            root,
            ["g05", "g01", "g03", "g02", "g06", "g04"],
            risks=risks,
            decisions=decisions,
            rationales=rationales,
            report_top="g02",
        )
        path = output / "ranked_guides.csv"
        rows = list(csv.DictReader(path.open(encoding="utf-8")))
        for row in rows:
            row["on_target_score"] = "0.99"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

    cases["scientifically_wrong"] = _run_case(
        checker,
        "wrong",
        wrong,
        lambda r: (
            r["core_science"] < 12
            and r["direction"] < 15
            and r["summary"] < 5
            and not r["hardgate_pass"]
            and "REPORT_TOP_CONTRADICTS_RANKING" in r["failure_codes"]
        ),
    )

    def alternative(root):
        risks = {"g01": "severe risk", "g02": "minimal", "g03": "moderate", "g04": "minimal risk", "g05": "critical", "g06": "lower"}
        decisions = {
            "g01": "deprioritize",
            "g02": "shortlist",
            "g03": "conditional review",
            "g04": "primary",
            "g05": "exclude",
            "g06": "recommended",
        }
        _write_fixture(
            root,
            # Safety-first g04 and activity-forward g03 placement intentionally
            # violate the former checker's disputable g02/g06/g04 relations.
            ["g04", "g03", "g02", "g06", "g01", "g05"],
            risks=risks,
            decisions=decisions,
            report_top="g04",
        )

    cases["valid_alternative_implementation"] = _run_case(
        checker,
        "alternative",
        alternative,
        lambda r: r["core_science"] == 40 and r["direction"] == 15 and r["summary"] == 5 and r["hardgate_pass"],
    )

    result = {
        "all_passed": all(case["passed"] for case in cases.values()),
        "candidate_code_executed": any(not case["candidate_code_not_executed"] for case in cases.values()),
        "cases": cases,
    }
    RESULT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
