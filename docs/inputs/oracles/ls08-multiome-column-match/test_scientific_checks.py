from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("ls08_scientific_checks", HERE / "scientific_checks.py")
CHECKER = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(CHECKER)
EXPECTED = [5, 1, 4, 0, 6, 3, 7, 2]

GOOD_REPORT = """# Result
The RNA-to-ATAC permutation is 5;1;4;0;6;3;7;2.
We summed ATAC bins overlapping each gene body to construct gene activity.
We compared normalized RNA expression and ATAC activity using gene-wise z-score correlation.
We used a global one-to-one assignment maximizing total similarity.
"""


def write_mapping(root: Path, values: list[int], *, labelled: bool = False, reverse: bool = False) -> None:
    out = root / "output"
    out.mkdir(parents=True, exist_ok=True)
    order = list(reversed(range(8))) if reverse else list(range(8))
    with (out / "column_mapping.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["match_score", "atac_column", "rna_population", "runner_up_score"])
        writer.writeheader()
        for i in order:
            writer.writerow({
                "rna_population": f"RNA_{i}" if labelled else i,
                "atac_column": f"ATAC-column-{values[i]}" if labelled else values[i],
                "match_score": 0.99,
                "runner_up_score": 0.01,
            })


def run_fixture(setup, assertions) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        setup(root)
        assertions(CHECKER.check(root))


def main() -> int:
    results = []

    def case(name, setup, assertions):
        try:
            run_fixture(setup, assertions)
            results.append({"case": name, "passed": True})
        except Exception as exc:
            results.append({"case": name, "passed": False, "error": f"{type(exc).__name__}: {exc}"})

    def correct(root):
        write_mapping(root, EXPECTED)
        (root / "output" / "report.md").write_text(GOOD_REPORT, encoding="utf-8")

    case("reference_like_correct", correct, lambda r: (
        (_ for _ in ()).throw(AssertionError(r))
        if not (r["core_science"] == 40 and r["direction"] == 15 and r["summary"] == 5 and r["hardgate_pass"])
        else None
    ))

    case("empty_or_missing", lambda root: None, lambda r: (
        (_ for _ in ()).throw(AssertionError(r))
        if not (r["core_science"] == r["direction"] == r["summary"] == 0 and not r["hardgate_pass"])
        else None
    ))

    def wrong(root):
        write_mapping(root, [(x + 1) % 8 for x in EXPECTED])
        (root / "output" / "report.md").write_text(
            "The score matrix has perfect 0.999 correlations and proves this mapping. " + GOOD_REPORT,
            encoding="utf-8",
        )

    case("scientifically_wrong", wrong, lambda r: (
        (_ for _ in ()).throw(AssertionError(r))
        if not (r["core_science"] == 0 and r["direction"] == 0 and r["summary"] == 0 and not r["hardgate_pass"]
                and "MAPPING_SCIENTIFICALLY_WRONG" in r["failure_codes"])
        else None
    ))

    def unsupported(root):
        out = root / "output"
        out.mkdir(parents=True)
        (out / "report.md").write_text(
            "I report permutation 5;1;4;0;6;3;7;2 and correlations 0.99, but did not use gene activity, "
            "did not compare RNA with ATAC, and fabricated the numbers.", encoding="utf-8"
        )

    def unsupported_assert(r):
        assert r["core_science"] == r["direction"] == r["summary"] == 0 and not r["hardgate_pass"], r
        # Negated-evidence variant: a real mapping cannot turn denial/fabrication into explanation credit.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_mapping(root, EXPECTED)
            (root / "output" / "report.md").write_text(
                "Permutation 5;1;4;0;6;3;7;2. We did not use gene activity or any RNA-ATAC biological signal; scores were fabricated.",
                encoding="utf-8",
            )
            negated = CHECKER.check(root)
            assert negated["direction"] == 9 and negated["summary"] == 0, negated
        # Fabricated-number variant: impressive self-reported scores do not rescue a false permutation.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            wrong(root)
            fabricated = CHECKER.check(root)
            assert fabricated["core_science"] == fabricated["direction"] == fabricated["summary"] == 0, fabricated

    case("plausible_but_unsupported", unsupported, unsupported_assert)

    def alternative(root):
        write_mapping(root, EXPECTED, labelled=True, reverse=True)
        (root / "output" / "report.md").write_text(
            "Alternative implementation: RNA-to-ATAC assignment 5,1,4,0,6,3,7,2. "
            "ATAC accessibility was aggregated in promoter/TSS windows linked to each gene. "
            "RNA TPM and ATAC gene activity were compared by Spearman similarity after normalization. "
            "A Hungarian one-to-one assignment found the global optimum.", encoding="utf-8"
        )

    case("valid_alternative_implementation", alternative, lambda r: (
        (_ for _ in ()).throw(AssertionError(r))
        if not (r["core_science"] == 40 and r["direction"] == 15 and r["summary"] == 5 and r["hardgate_pass"])
        else None
    ))

    payload = {"schema_version": 1, "all_passed": all(x["passed"] for x in results), "cases": results}
    (HERE / "acceptance-result.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0 if payload["all_passed"] and len(results) == 5 else 1


if __name__ == "__main__":
    raise SystemExit(main())
