from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("ls07_scientific_checks", HERE / "scientific_checks.py")
CHECKER = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(CHECKER)


def _fixture(table=None, call=None, report=None, delimiter=","):
    root = Path(tempfile.mkdtemp(prefix="ls07-acceptance-"))
    out = root / "output"
    out.mkdir()
    if table is not None:
        headers, rows = table
        with (out / "pathway_enrichment.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=headers, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(rows)
    if call is not None:
        (out / "mechanism_call.json").write_text(json.dumps(call), encoding="utf-8")
    if report is not None:
        (out / "report.md").write_text(report, encoding="utf-8")
    return root


HEADERS = ["pathway_id", "pathway_name", "overlap", "p_value", "padj", "direction"]
GOOD_ROW = {
    "pathway_id": "R-HSA-6791312",
    "pathway_name": "TP53 Regulates Transcription Of Cell Cycle Genes",
    "overlap": "8/49",
    "p_value": "0.00014600658",
    "padj": "0.012",
    "direction": "mixed TP53 response with cell-cycle repression",
}
GOOD_CALL = {
    "primary_mechanism": "TP53-mediated cell cycle regulation",
    "direction": "checkpoint induction with cell-cycle arrest",
}
GOOD_REPORT = (
    "Reactome enrichment supports TP53-mediated cell-cycle regulation, with checkpoint induction "
    "and cell-cycle repression. This enrichment is consistent with the mechanism hypothesis but "
    "does not prove causation."
)


def run():
    results = {}

    correct = CHECKER.check(_fixture((HEADERS, [GOOD_ROW]), GOOD_CALL, GOOD_REPORT))
    results["reference_like_correct"] = (
        correct["core_science"] == 40 and correct["direction"] == 15
        and correct["summary"] == 5 and correct["hardgate_pass"]
    )

    missing = CHECKER.check(_fixture())
    results["empty_or_missing"] = (
        missing["core_science"] == 0 and missing["direction"] == 0
        and missing["summary"] == 0 and not missing["hardgate_pass"]
    )

    wrong_call = CHECKER.check(_fixture(
        (HEADERS, [GOOD_ROW]),
        {"primary_mechanism": "NRF2-mediated oxidative stress", "direction": "increased proliferation"},
        "Reactome enrichment establishes NRF2 as the causal mechanism and activates the cell cycle.",
    ))
    negated = CHECKER.check(_fixture(
        (HEADERS, [GOOD_ROW]),
        {"primary_mechanism": "TP53 cell cycle is not the mechanism; NRF2 is the mechanism"},
        "TP53 cell-cycle enrichment is not supported; NRF2 is instead responsible.",
    ))
    fabricated = dict(GOOD_ROW, overlap="17/49", p_value="1e-40", padj="1e-30")
    fabricated_result = CHECKER.check(_fixture((HEADERS, [fabricated]), GOOD_CALL, GOOD_REPORT))
    causal = CHECKER.check(_fixture(
        (HEADERS, [GOOD_ROW]), GOOD_CALL,
        "Reactome enrichment proves TP53 is the causal mechanism driving cell-cycle arrest.",
    ))
    results["scientifically_wrong"] = (
        not wrong_call["hardgate_pass"] and wrong_call["core_science"] < 40
        and not negated["hardgate_pass"] and negated["core_science"] < 40
        and not fabricated_result["hardgate_pass"] and fabricated_result["core_science"] < 40
        and not causal["hardgate_pass"] and "CAUSAL_OVERCLAIM" in causal["failure_codes"]
    )

    unsupported = CHECKER.check(_fixture(
        (HEADERS, [{
            "pathway_id": "R-HSA-199420", "pathway_name": "Apoptosis", "overlap": "8/49",
            "p_value": "0.0008", "padj": "0.012", "direction": "downregulated",
        }]),
        GOOD_CALL,
        GOOD_REPORT,
    ))
    results["plausible_but_unsupported"] = (
        not unsupported["hardgate_pass"] and unsupported["core_science"] <= 12
        and unsupported["direction"] == 0
    )

    alt_headers = ["Term_ID", "Term", "Overlap_Fraction", "PValue", "FDR", "Regulation"]
    alt_row = {
        "Term_ID": "R-HSA-6791312",
        "Term": "TP53 Regulates Transcription Of Cell Cycle Genes",
        "Overlap_Fraction": "8 / 49",
        "PValue": "1.4600658e-4",
        "FDR": "1.2e-2",
        "Regulation": "bidirectional p53 targets; proliferation suppressed",
    }
    alternative = CHECKER.check(_fixture(
        (alt_headers, [alt_row]),
        {"conclusion": {"call": "p53-dependent cell-cycle arrest"}},
        "The enrichment is consistent with p53-dependent cell-cycle arrest; it cannot demonstrate causation.",
        delimiter=";",
    ))
    results["valid_alternative_implementation"] = (
        alternative["core_science"] == 40 and alternative["direction"] == 15
        and alternative["summary"] == 5 and alternative["hardgate_pass"]
    )

    return results


if __name__ == "__main__":
    outcome = run()
    print(json.dumps({"passed": all(outcome.values()), "cases": outcome}, indent=2))
    raise SystemExit(0 if all(outcome.values()) else 1)
