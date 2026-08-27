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

ANCHORS = CHECKER.ANCHORS


def _run(build_call=None, report=None):
    with tempfile.TemporaryDirectory(prefix="ls02-accept-") as temp:
        workspace = Path(temp)
        output = workspace / "output"
        output.mkdir()
        if build_call is not None:
            (output / "build_call.json").write_text(json.dumps(build_call), encoding="utf-8")
        if report is not None:
            (output / "report.md").write_text(report, encoding="utf-8")
        return CHECKER.check(workspace)


def _standard_evidence(limit=10):
    rows = []
    for pos, (vcf_ref, hg18, hg19, hg38, hs1) in list(ANCHORS.items())[:limit]:
        rows.append({
            "coordinate": f"chr20:{pos}",
            "vcf_ref": vcf_ref,
            "reference_bases": {"hg18": hg18, "hg19": hg19, "hg38": hg38, "T2T": hs1},
        })
    return rows


def _report(label="hg19"):
    return (
        f"# Genome build call\n\nThe supported build is {label}. We compared the VCF REF allele at each "
        "1-based coordinate to assembly reference sequence.\n\n"
        "- 61795 VCF REF G; hg19 reference G\n"
        "- 3752643 VCF REF T; hg19 reference T\n"
    )


def reference_like_correct():
    result = _run({
        "build": "hg19", "confidence": 0.99,
        "n_variants_checked": 18, "n_ref_matches": 18, "n_ref_mismatches": 0,
        "evidence": _standard_evidence(),
    }, _report())
    assert result["hardgate_pass"] is True
    assert (result["core_science"], result["direction"], result["summary"]) == (40, 15, 5)
    return result


def empty_or_missing():
    result = _run()
    assert result["hardgate_pass"] is False
    assert (result["core_science"], result["direction"], result["summary"]) == (0, 0, 0)
    return result


def scientifically_wrong():
    result = _run({
        "build": "hg38", "confidence": "high", "evidence": _standard_evidence(),
    }, _report("hg38"))
    assert result["hardgate_pass"] is False
    assert "FATAL_TRUTH_CONCLUSION" in result["failure_codes"]
    assert result["direction"] < 15
    assert result["criteria"]["report_concludes_hg19"] is False
    return result


def plausible_but_unsupported():
    copied = [{"position": pos, "vcf_ref": values[0]} for pos, values in list(ANCHORS.items())[:10]]
    unsupported = _run({
        "build": "hg19", "confidence": 0.999,
        "n_variants_checked": 84664, "n_ref_matches": 84664, "n_ref_mismatches": 0,
        "match_scores": {"hg18": 4, "hg19": 18, "hg38": 4, "T2T": 4},
        "evidence": copied,
    }, _report())
    assert unsupported["hardgate_pass"] is False
    assert unsupported["core_science"] <= 6

    fabricated_numbers = _run({
        "build": "hg19", "confidence": "high",
        "n_variants_checked": 18, "n_ref_matches": 18,
        "evidence": [{"position": 123456, "vcf_ref": "A", "hg19": "A"}],
    }, _report())
    assert fabricated_numbers["hardgate_pass"] is False
    assert fabricated_numbers["criteria"]["verified_anchor_count"] == 0

    negated = _standard_evidence()
    for row in negated:
        row["comparisons"] = {"hg19": False}
    negated_result = _run({"build": "hg19", "confidence": "high", "evidence": negated}, _report())
    assert negated_result["hardgate_pass"] is False
    assert "FATAL_GROUNDED_POSITIVE_EVIDENCE" in negated_result["failure_codes"]
    return {
        "unsupported_self_report": unsupported,
        "fabricated_number_variant": fabricated_numbers,
        "negated_evidence_variant": negated_result,
    }


def valid_alternative_implementation():
    loci = []
    for pos, (vcf_ref, hg18, hg19, hg38, hs1) in list(ANCHORS.items())[:10]:
        loci.append({
            "locus": f"chr20:{pos}", "observed_ref": vcf_ref,
            "comparisons": [
                {"assembly": "NCBI36", "base": hg18},
                {"assembly": "GRCh37", "base": hg19},
                {"assembly": "GRCh38", "base": hg38},
                {"assembly": "CHM13", "base": hs1},
            ],
        })
    result = _run({
        "build": "GRCh37", "confidence": "high",
        "evidence": {"method": "indexed FASTA allele lookup", "loci": loci},
    }, _report("GRCh37"))
    assert result["hardgate_pass"] is True
    assert (result["core_science"], result["direction"], result["summary"]) == (40, 15, 5)
    return result


def main():
    cases = [
        ("reference_like_correct", reference_like_correct),
        ("empty_or_missing", empty_or_missing),
        ("scientifically_wrong", scientifically_wrong),
        ("plausible_but_unsupported", plausible_but_unsupported),
        ("valid_alternative_implementation", valid_alternative_implementation),
    ]
    results = []
    for name, case in cases:
        try:
            detail = case()
            results.append({"case": name, "passed": True, "detail": detail})
        except Exception as exc:
            results.append({"case": name, "passed": False, "error": f"{type(exc).__name__}: {exc}"})
    payload = {
        "schema_version": 1,
        "all_passed": all(item["passed"] for item in results) and len(results) == 5,
        "cases_passed": sum(item["passed"] for item in results),
        "cases_required": 5,
        "candidate_code_executed": False,
        "results": results,
    }
    (HERE / "acceptance-result.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    raise SystemExit(0 if payload["all_passed"] else 1)


if __name__ == "__main__":
    main()
