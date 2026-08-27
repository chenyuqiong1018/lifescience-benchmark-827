from __future__ import annotations

import csv
import json
import tempfile
from pathlib import Path

import scientific_checks as checker


CASE_NAMES = (
    "reference_like_correct",
    "empty_or_missing",
    "scientifically_wrong",
    "plausible_but_unsupported",
    "valid_alternative_implementation",
)


def _write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _trap_candidate_code(workspace: Path) -> Path:
    marker = workspace / "candidate_code_was_executed"
    script = workspace / "output" / "analysis.py"
    script.parent.mkdir(parents=True, exist_ok=True)
    script.write_text(
        "from pathlib import Path\n"
        f"Path({str(marker)!r}).write_text('executed', encoding='utf-8')\n"
        "raise RuntimeError('candidate code must never run')\n",
        encoding="utf-8",
    )
    return marker


def _standard_cells(
    *,
    distance_delta: float = 0.0,
    repeated: tuple[float, float, float] | None = None,
) -> list[dict]:
    rows = []
    for cell, expected in sorted(checker._truth()["cell_metrics"].items()):
        distance = expected["mean_distance_nm"] + distance_delta
        contact = expected["contact_fraction"]
        transcription = expected["transcription_fraction"]
        if repeated is not None:
            distance, contact, transcription = repeated
        rows.append(
            {
                "cell_id": cell,
                "mean_3d_distance_nm": distance,
                "contact_fraction": contact,
                "transcription_fraction": transcription,
            }
        )
    return rows


def _standard_lags(
    method: str = "contact_pearson",
    orientation: int = 1,
    transform=None,
) -> list[dict]:
    rows = []
    for lag in range(-5, 6):
        expected = checker._lag_truth(orientation * lag)
        association = expected[method]
        n_observations = expected["n_observations"]
        if transform is not None:
            association, n_observations = transform(
                lag, association, n_observations
            )
        rows.append(
            {
                "lag": lag,
                "association": association,
                "n_observations": n_observations,
            }
        )
    return rows


def _write_standard(workspace: Path, cells: list[dict], lags: list[dict], report: str) -> None:
    _write_csv(
        workspace / "output" / "cell_metrics.csv",
        [
            "cell_id",
            "mean_3d_distance_nm",
            "contact_fraction",
            "transcription_fraction",
        ],
        cells,
    )
    _write_csv(
        workspace / "output" / "lag_analysis.csv",
        ["lag", "association", "n_observations"],
        lags,
    )
    (workspace / "output" / "report.md").write_text(report, encoding="utf-8")


def _correct_report(extra: str = "") -> str:
    return (
        "# Results\n\n"
        "Euclidean enhancer-promoter distance averaged 526.067 nm. "
        "Using the supplied inclusive 260 nm contact threshold, the contact "
        "fraction was 5.336%, and the transcription fraction was 15.391%. "
        "The strongest temporal association was at lag -1 under the stated "
        "contact-at-t versus transcription-at-t+lag convention. "
        "These are observational data and cannot establish causality or causal "
        "direction; causal inference is not possible from this design. "
        + extra
    )


def _case_reference_like_correct(workspace: Path) -> tuple[bool, dict]:
    marker = _trap_candidate_code(workspace)
    _write_standard(
        workspace,
        _standard_cells(),
        _standard_lags(),
        _correct_report(),
    )
    result = checker.check(workspace)
    passed = (
        result["hardgate_pass"]
        and result["core_science"] == 40
        and result["direction"] == 15
        and result["summary"] == 5
        and not marker.exists()
    )
    return passed, {
        "hardgate_pass": result["hardgate_pass"],
        "scores": [
            result["core_science"], result["direction"], result["summary"]
        ],
        "candidate_code_not_executed": not marker.exists(),
    }


def _case_empty_or_missing(workspace: Path) -> tuple[bool, dict]:
    marker = _trap_candidate_code(workspace)
    result = checker.check(workspace)
    passed = (
        not result["hardgate_pass"]
        and result["core_science"] == 0
        and result["direction"] == 0
        and result["summary"] == 0
        and not marker.exists()
    )
    return passed, {
        "hardgate_pass": result["hardgate_pass"],
        "scores": [
            result["core_science"], result["direction"], result["summary"]
        ],
        "candidate_code_not_executed": not marker.exists(),
    }


def _case_scientifically_wrong(workspace: Path) -> tuple[bool, dict]:
    marker = _trap_candidate_code(workspace)
    wrong_cells = []
    for row in _standard_cells(distance_delta=100.0):
        row["contact_fraction"] = 1.0 - row["contact_fraction"]
        row["transcription_fraction"] = 0.75
        wrong_cells.append(row)
    wrong_lags = _standard_lags(
        transform=lambda lag, association, n: (-association, n + 600)
    )
    _write_standard(
        workspace,
        wrong_cells,
        wrong_lags,
        _correct_report(
            "The tables intentionally contain a scientifically wrong calculation."
        ),
    )
    result = checker.check(workspace)
    gate = result["criteria"]["fatal_gate_CORE_RESULTS_GROUNDED"]
    passed = (
        not result["hardgate_pass"]
        and not gate["passed"]
        and result["core_science"] < 20
        and not marker.exists()
    )
    return passed, {
        "hardgate_pass": result["hardgate_pass"],
        "core_science": result["core_science"],
        "grounded_gate_passed": gate["passed"],
        "candidate_code_not_executed": not marker.exists(),
    }


def _case_plausible_but_unsupported(workspace: Path) -> tuple[bool, dict]:
    marker = _trap_candidate_code(workspace)
    repeated = (526.0669983528302, 0.05336, 0.15390666666666666)
    plausible_lags = _standard_lags(
        transform=lambda lag, association, n: (0.05336, n)
    )
    _write_standard(
        workspace,
        _standard_cells(repeated=repeated),
        plausible_lags,
        _correct_report(
            "Despite that limitation, enhancer-promoter contact causes "
            "transcription."
        ),
    )
    plausible = checker.check(workspace)

    with tempfile.TemporaryDirectory(prefix="ls03-fabricated-") as temp:
        fabricated_workspace = Path(temp)
        fabricated_marker = _trap_candidate_code(fabricated_workspace)
        fabricated_lags = _standard_lags(
            transform=lambda lag, association, n: (0.0, 150000)
        )
        _write_standard(
            fabricated_workspace,
            _standard_cells(repeated=(500.0, 0.10, 0.20)),
            fabricated_lags,
            _correct_report(
                "We cannot say contact causes transcription; there is no "
                "evidence here that proves a causal effect."
            ),
        )
        fabricated = checker.check(fabricated_workspace)
        fabricated_code_safe = not fabricated_marker.exists()

    causal_gate = plausible["criteria"]["fatal_gate_CAUSAL_INTERPRETATION_SAFE"]
    grounded_gate = plausible["criteria"]["fatal_gate_CORE_RESULTS_GROUNDED"]
    fabricated_grounded = fabricated["criteria"][
        "fatal_gate_CORE_RESULTS_GROUNDED"
    ]
    passed = (
        not plausible["hardgate_pass"]
        and not causal_gate["passed"]
        and not grounded_gate["passed"]
        and bool(causal_gate["unsupported_causal_claims"])
        and not fabricated["hardgate_pass"]
        and not fabricated_grounded["passed"]
        and not marker.exists()
        and fabricated_code_safe
    )
    return passed, {
        "self_report_rejected": not plausible["hardgate_pass"],
        "bare_causal_claim_rejected": not causal_gate["passed"],
        "fabricated_numbers_rejected": not fabricated["hardgate_pass"],
        "negated_evidence_not_misclassified": fabricated["criteria"][
            "fatal_gate_CAUSAL_INTERPRETATION_SAFE"
        ]["passed"],
        "candidate_code_not_executed": not marker.exists() and fabricated_code_safe,
    }


def _case_valid_alternative_implementation(workspace: Path) -> tuple[bool, dict]:
    marker = _trap_candidate_code(workspace)
    cell_rows = []
    for cell, expected in sorted(checker._truth()["cell_metrics"].items()):
        cell_rows.append(
            {
                "Cell": cell,
                "Average Distance (nm)": expected["mean_distance_nm"],
                "Contact Percent": 100.0 * expected["contact_fraction"],
                "Active Fraction": expected["transcription_fraction"],
            }
        )
    lag_rows = []
    for lag in range(-5, 6):
        expected = checker._lag_truth(-lag)
        lag_rows.append(
            {
                "Time Lag": lag,
                "Correlation": expected["proximity_pearson"],
                "N Pairs": expected["n_observations"],
            }
        )
    _write_csv(
        workspace / "output" / "cell_metrics.csv",
        ["Cell", "Average Distance (nm)", "Contact Percent", "Active Fraction"],
        cell_rows,
    )
    _write_csv(
        workspace / "output" / "lag_analysis.csv",
        ["Time Lag", "Correlation", "N Pairs"],
        lag_rows,
    )
    (workspace / "output" / "report.md").write_text(
        "# Alternative analysis\n\n"
        "The mean 3D distance was 526.067 nm; at the inclusive 260 nm cutoff, "
        "5.336% were contacts and 15.391% were transcription-positive. "
        "With the opposite lag-sign convention, the strongest temporal "
        "association occurs at lag +1. We cannot say contact causes "
        "transcription. Causal inference is not possible from these "
        "observational data, and directionality remains unresolved.",
        encoding="utf-8",
    )
    result = checker.check(workspace)
    causal = result["criteria"]["fatal_gate_CAUSAL_INTERPRETATION_SAFE"]
    lag = result["criteria"]["lag_analysis_recomputed"]
    passed = (
        result["hardgate_pass"]
        and result["core_science"] == 40
        and result["direction"] == 15
        and result["summary"] == 5
        and causal["passed"]
        and not causal["unsupported_causal_claims"]
        and lag["best_method"] == "proximity_pearson"
        and lag["best_orientation"] == "opposite_lag_sign_convention"
        and not marker.exists()
    )
    return passed, {
        "hardgate_pass": result["hardgate_pass"],
        "scores": [
            result["core_science"], result["direction"], result["summary"]
        ],
        "negated_causal_wording_accepted": causal["passed"],
        "alternative_method": lag["best_method"],
        "alternative_orientation": lag["best_orientation"],
        "candidate_code_not_executed": not marker.exists(),
    }


CASES = {
    "reference_like_correct": _case_reference_like_correct,
    "empty_or_missing": _case_empty_or_missing,
    "scientifically_wrong": _case_scientifically_wrong,
    "plausible_but_unsupported": _case_plausible_but_unsupported,
    "valid_alternative_implementation": _case_valid_alternative_implementation,
}


def main() -> int:
    case_results = []
    for name in CASE_NAMES:
        with tempfile.TemporaryDirectory(prefix=f"ls03-{name}-") as temp:
            try:
                passed, details = CASES[name](Path(temp))
                error = None
            except Exception as exc:
                passed, details, error = False, {}, f"{type(exc).__name__}: {exc}"
            case_results.append(
                {
                    "name": name,
                    "passed": bool(passed),
                    "details": details,
                    "error": error,
                }
            )
    result = {
        "schema_version": 1,
        "task_id": "ls03-genome-coordinates",
        "passed": all(case["passed"] for case in case_results),
        "cases_passed": sum(case["passed"] for case in case_results),
        "cases_total": len(case_results),
        "candidate_code_executed": False,
        "cases": case_results,
    }
    output = Path(__file__).with_name("acceptance-result.json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
