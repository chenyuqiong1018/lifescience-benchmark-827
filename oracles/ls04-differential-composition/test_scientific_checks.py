from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("ls04_scientific_checks", HERE / "scientific_checks.py")
CHECKER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECKER)


def _write_fixture(root, rows=None, call=None, report=None):
    output = root / "output"
    output.mkdir(parents=True)
    if rows is not None:
        with (output / "composition.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=("sample", "cell_type", "n_cells", "fraction"))
            writer.writeheader()
            writer.writerows(rows)
    if call is not None:
        (output / "depleted_call.json").write_text(json.dumps(call), encoding="utf-8")
    if report is not None:
        (output / "report.md").write_text(report, encoding="utf-8")


def _rows(h1=179, h2=4, total1=6295, total2=5004, labels=("sample 1", "sample 2")):
    return [
        {"sample": labels[0], "cell_type": "Horizontal cells", "n_cells": h1, "fraction": f"{h1/total1:.12f}"},
        {"sample": labels[0], "cell_type": "Photoreceptors", "n_cells": total1-h1, "fraction": f"{(total1-h1)/total1:.12f}"},
        {"sample": labels[1], "cell_type": "Horizontal cells", "n_cells": h2, "fraction": f"{h2/total2:.12f}"},
        {"sample": labels[1], "cell_type": "Photoreceptors", "n_cells": total2-h2, "fraction": f"{(total2-h2)/total2:.12f}"},
    ]


GOOD_CALL = {"depleted_cell_type": "horizontal cells", "depleted_sample": "sample 2"}
GOOD_REPORT = (
    "Horizontal cells are severely depleted in sample 2. Annotation evidence uses PROX1, LHX1, "
    "and ONECUT1 markers. Low-count sampling and QC sensitivity are important uncertainty limitations."
)


def _run(rows=None, call=None, report=None):
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        if rows is not None or call is not None or report is not None:
            _write_fixture(root, rows, call, report)
        return CHECKER.check(root)


def reference_like_correct():
    result = _run(_rows(), GOOD_CALL, GOOD_REPORT)
    assert result["hardgate_pass"] and (result["core_science"], result["direction"], result["summary"]) == (40, 15, 5)


def empty_or_missing():
    result = _run()
    assert not result["hardgate_pass"] and result["core_science"] == 0 and result["direction"] == 0


def scientifically_wrong():
    negated = (
        "Horizontal cells are not depleted in sample 2; the claim is unsupported. "
        "No horizontal-cell marker evidence was established. QC uncertainty remains."
    )
    result = _run(_rows(h1=4, h2=179), GOOD_CALL, negated)
    assert not result["hardgate_pass"] and result["direction"] < 8
    assert not result["criteria"]["report_affirmed_conclusion"]
    assert not result["criteria"]["report_marker_evidence"]


def plausible_but_unsupported():
    # Reviewer regression: a fabricated fraction must not rescue a one-cell count.
    contradictory = _rows(h1=1, h2=0)
    contradictory[0]["fraction"] = "0.028"
    contradictory[1]["fraction"] = "0.972"
    first = _run(contradictory, GOOD_CALL, GOOD_REPORT)
    assert not first["hardgate_pass"] and not first["criteria"]["count_fraction_arithmetic"]

    # Fabricated totals remain unsupported even when the horizontal counts look plausible.
    fabricated_totals = _rows(h1=179, h2=4, total1=7000, total2=6000)
    second = _run(fabricated_totals, GOOD_CALL, GOOD_REPORT)
    assert not second["hardgate_pass"] and not second["criteria"]["grounded_sample_totals"]

    # Reviewer regression: an explicit opposite direction in the call must be fatal.
    opposite = _run(_rows(), {"depleted_cell_type": "horizontal", "depleted_sample": "sample 1"}, GOOD_REPORT)
    assert not opposite["hardgate_pass"] and not opposite["criteria"]["depleted_call_explicitly_sample2"]


def valid_alternative_implementation():
    rows = [
        {"sample": "q1.2", "cell_type": "Rod", "n_cells": 5000, "fraction": f"{5000/5004:.10f}"},
        {"sample": "q1.1", "cell_type": "Horizontal cell H2", "n_cells": 79, "fraction": f"{79/6295:.10f}"},
        {"sample": "q1.1", "cell_type": "Other retinal populations", "n_cells": 6116, "fraction": f"{6116/6295:.10f}"},
        {"sample": "q1.2", "cell_type": "Horizontal cell H1", "n_cells": 4, "fraction": f"{4/5004:.10f}"},
        {"sample": "q1.1", "cell_type": "Horizontal cell H1", "n_cells": 100, "fraction": f"{100/6295:.10f}"},
    ]
    call = {"depleted": {"population": "HC / horizontal-cell family", "depleted_in": "q1.2"},
            "method": "marker-score clustering"}
    report = (
        "The horizontal-cell family shows a severe abundance collapse in the second sample. "
        "We annotated it using marker evidence from LHX1, ONECUT2, and PROX1 expression. "
        "Annotation sensitivity and low-count sampling are limitations."
    )
    result = _run(rows, call, report)
    assert result["hardgate_pass"] and (result["core_science"], result["direction"], result["summary"]) == (40, 15, 5)


CASES = [reference_like_correct, empty_or_missing, scientifically_wrong,
         plausible_but_unsupported, valid_alternative_implementation]


if __name__ == "__main__":
    results = []
    for case in CASES:
        try:
            case()
            results.append({"case": case.__name__, "passed": True})
        except Exception as exc:
            results.append({"case": case.__name__, "passed": False,
                            "error": f"{type(exc).__name__}: {exc}"})
    payload = {"all_passed": all(item["passed"] for item in results), "cases": results,
               "candidate_code_executed": False}
    print(json.dumps(payload, indent=2))
    raise SystemExit(0 if payload["all_passed"] else 1)
