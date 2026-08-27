from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("ls08_scientific_checks", HERE / "scientific_checks.py")
CHECKER = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(CHECKER)


def _write_fixture(
    root: Path,
    rows: list[dict[str, object]] | None,
    call: dict[str, object] | None,
    report: str | None,
) -> None:
    output = root / "output"
    output.mkdir(parents=True)
    if rows is not None:
        with (output / "pair_evidence.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    if call is not None:
        (output / "least_supported.json").write_text(json.dumps(call), encoding="utf-8")
    if report is not None:
        (output / "report.md").write_text(report, encoding="utf-8")
    # Deliberately hostile candidate code: the acceptance suite must never execute it.
    (output / "analysis.py").write_text("raise RuntimeError('candidate code executed')\n", encoding="utf-8")


def _truth_rows(alternative: bool = False) -> list[dict[str, object]]:
    truth = CHECKER._truth()
    ordered = sorted(truth, key=lambda pair_id: (truth[pair_id]["combined"], pair_id))
    final_rank = {pair_id: index for index, pair_id in enumerate(ordered, 1)}
    rows = []
    for pair_id in reversed(CHECKER.EXPECTED_IDS):
        item = truth[pair_id]
        if alternative:
            rows.append(
                {
                    "id": pair_id.lower(),
                    "hic_residual": f"{item['contact']:.8f}",
                    "absolute_log2fc": f"{item['effect']:.8f}",
                    "mean_rank": f"{item['combined']:.6f}",
                    "support_rank": final_rank[pair_id],
                }
            )
        else:
            rows.append(
                {
                    "pair_id": pair_id,
                    "contact_evidence": item["contact"],
                    "perturbation_effect": item["effect"],
                    "combined_support": item["combined"],
                    "rank": final_rank[pair_id],
                }
            )
    return rows


def _run(rows, call, report):
    with tempfile.TemporaryDirectory(prefix="ls08-acceptance-") as tmp:
        root = Path(tmp)
        _write_fixture(root, rows, call, report)
        return CHECKER.check(root)


def main() -> int:
    outcomes: dict[str, dict[str, object]] = {}

    correct = _run(
        _truth_rows(),
        {"pair_id": "EP6", "choice": "F"},
        "EP6 is the least supported causal pair. Its Hi-C contact residual is -2.293, "
        "while the CRISPR perturbation effect is 0.249 absolute log2 fold change; these "
        "are distinct evidence modalities.",
    )
    outcomes["reference_like_correct"] = {
        "passed": correct["hardgate_pass"] and (correct["core_science"], correct["direction"], correct["summary"]) == (40, 15, 5),
        "scores": [correct["core_science"], correct["direction"], correct["summary"]],
    }

    empty = _run(None, None, None)
    outcomes["empty_or_missing"] = {
        "passed": not empty["hardgate_pass"] and (empty["core_science"], empty["direction"], empty["summary"]) == (0, 0, 0),
        "scores": [empty["core_science"], empty["direction"], empty["summary"]],
    }

    wrong_rows = _truth_rows()
    for row in wrong_rows:
        row["contact_evidence"] = -float(row["contact_evidence"])
        row["perturbation_effect"] = 9.0 - float(row["perturbation_effect"])
        row["combined_support"] = 8.0 - float(row["combined_support"])
        row["rank"] = 8 - int(row["rank"])
    wrong = _run(
        wrong_rows,
        {"pair_id": "EP1"},
        "EP1 is the least supported pair based on contact and perturbation evidence.",
    )
    outcomes["scientifically_wrong"] = {
        "passed": not wrong["hardgate_pass"] and wrong["core_science"] == 0 and wrong["direction"] == 0,
        "scores": [wrong["core_science"], wrong["direction"], wrong["summary"]],
    }

    fabricated_rows = [
        {
            "pair_id": pair_id,
            "contact_evidence": 1000 + i,
            "perturbation_effect": 2000 + i,
            "combined_support": i,
            "rank": i,
        }
        for i, pair_id in enumerate(CHECKER.EXPECTED_IDS, 1)
    ]
    # Make the candidate's repeated self-report internally agree on EP6.
    fabricated_rows[5]["combined_support"] = 0
    fabricated_rows[5]["rank"] = 0
    unsupported = _run(
        fabricated_rows,
        {"least_supported": "EP6"},
        "EP6 is the least supported pair. Hi-C contact residual 999 and CRISPR effect 999 prove it.",
    )
    negated = _run(
        _truth_rows(),
        {"pair_id": "EP6"},
        "EP6 is not the least supported pair; contact and perturbation instead favor EP5.",
    )
    outcomes["plausible_but_unsupported"] = {
        "passed": (
            not unsupported["hardgate_pass"]
            and unsupported["core_science"] == 0
            and unsupported["summary"] == 0
            and negated["summary"] == 0
        ),
        "variants": {
            "fabricated_numbers_rejected": unsupported["summary"] == 0,
            "negated_evidence_rejected": negated["summary"] == 0,
        },
        "scores": [unsupported["core_science"], unsupported["direction"], unsupported["summary"]],
    }

    alternative = _run(
        _truth_rows(alternative=True),
        {"least_supported_pair": "ep6"},
        "The weakest integrated candidate is EP6. Physical contact and expression-change "
        "evidence were evaluated separately before rank integration.",
    )
    outcomes["valid_alternative_implementation"] = {
        "passed": alternative["hardgate_pass"] and (alternative["core_science"], alternative["direction"], alternative["summary"]) == (40, 15, 5),
        "scores": [alternative["core_science"], alternative["direction"], alternative["summary"]],
    }

    result = {
        "all_passed": all(bool(item["passed"]) for item in outcomes.values()),
        "cases": outcomes,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
