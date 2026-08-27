from __future__ import annotations

import csv
import json
import math
import shutil
import tempfile
from pathlib import Path

import scientific_checks as checker


def build_rows(alternative=False, wrong_direction=False, fabricated=False):
    _, truth = checker._truth()
    selected = [(g, b, s) for g, (b, s) in truth.items() if b > 11 and s]
    selected.sort(key=lambda x: x[0])
    passing, controls = selected[:677], selected[677:727]
    fields = (["symbol", "adjusted_p_value", "ensembl_id", "is_significant", "log2fc", "p_value", "base_mean"]
              if alternative else ["gene_id", "gene_name", "baseMean", "log2FoldChange", "pvalue", "padj", "pass"])
    rows = []
    for index, (gene, base, sign) in enumerate(passing + controls):
        is_pass = index < 677
        lfc = sign * (1.0 if is_pass else 0.1)
        if wrong_direction:
            lfc *= -1
        if fabricated:
            base += 123.456
        values = {
            "gene_id": gene, "gene_name": "G" + gene[-5:], "baseMean": base,
            "log2FoldChange": lfc, "pvalue": 0.001 if is_pass else 0.8,
            "padj": 0.002 if is_pass else 0.9, "pass": is_pass,
        }
        if alternative:
            values = {"ensembl_id": values["gene_id"], "symbol": values["gene_name"],
                      "base_mean": values["baseMean"], "log2fc": values["log2FoldChange"],
                      "p_value": values["pvalue"], "adjusted_p_value": values["padj"],
                      "is_significant": 1 if is_pass else 0}
        rows.append({k: values[k] for k in fields})
    return fields, rows


def write_fixture(root, mode):
    out = root / "output"
    out.mkdir(parents=True)
    if mode == "empty":
        return
    alternative = mode == "alternative"
    fields, rows = build_rows(alternative=alternative, wrong_direction=mode == "wrong",
                              fabricated=mode == "fabricated")
    with (out / "differential_expression.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    if alternative:
        summary = {"significant_gene_count": 677}
        report = "Combination Cisplatin + CBD relative to DMSO control; positive log2FC denotes higher expression in combination."
    else:
        summary = {"n_passing": 677, "padj_lt": 0.05, "abs_log2fc_gt": 0.5, "baseMean_gt": 10}
        report = "Cisplatin_IC50_CBD_IC50 compared to DMSO control; positive log2 fold change is higher in combination."
    if mode == "unsupported":
        (out / "differential_expression.csv").unlink()
        summary = {"n_passing": 677, "claim": "DESeq2-complete"}
    if mode == "negated":
        report = "This does not compare the combination to DMSO control. The direction is unknown."
    (out / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    (out / "report.md").write_text(report, encoding="utf-8")


def run_case(name, mode, expected_pass):
    root = Path(tempfile.mkdtemp(prefix="ls07_accept_"))
    try:
        write_fixture(root, mode)
        result = checker.check(root)
        actual = result["hardgate_pass"]
        return {"name": name, "passed": actual == expected_pass, "expected_hardgate": expected_pass,
                "actual_hardgate": actual, "failure_codes": result["failure_codes"]}
    finally:
        shutil.rmtree(root)


def main():
    cases = [
        run_case("reference_like_correct", "correct", True),
        run_case("empty_or_missing", "empty", False),
        run_case("scientifically_wrong", "wrong", False),
        run_case("plausible_but_unsupported", "unsupported", False),
        run_case("valid_alternative_implementation", "alternative", True),
        run_case("negated_evidence_variant", "negated", False),
        run_case("fabricated_number_variant", "fabricated", False),
    ]
    result = {"all_required_cases_passed": all(c["passed"] for c in cases[:5]),
              "adversarial_false_positive_passed": all(c["passed"] for c in cases[5:]), "cases": cases,
              "candidate_code_executed": False}
    Path(__file__).with_name("acceptance-result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if all(c["passed"] for c in cases) else 1)


if __name__ == "__main__":
    main()
