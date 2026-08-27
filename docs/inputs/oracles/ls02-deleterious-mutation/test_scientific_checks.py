from __future__ import annotations

import importlib.util
import json
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ls02_science", HERE / "scientific_checks.py")
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


def write_fixture(root, variant=None, evidence=None, report=None):
    out = root / "output"
    out.mkdir(parents=True)
    if variant is not None:
        (out / "variant.tsv").write_text(variant, encoding="utf-8")
    if evidence is not None:
        (out / "evidence.json").write_text(json.dumps(evidence), encoding="utf-8")
    if report is not None:
        (out / "report.md").write_text(report, encoding="utf-8")
    # Deliberately inert; the acceptance suite never imports, runs, or inspects it.
    (out / "analysis.py").write_text("raise RuntimeError('must never execute')\n", encoding="utf-8")


CORRECT_REPORT = """# Result
The high-confidence mosaic call is **chr9:127661125 G>T** in **STXBP1**, a
loss-of-function-intolerant gene. It is a nonsense/stop-gained change,
p.Glu117Ter. There are 18 alternate reads among 93 quality-filtered reads
(allele fraction 0.1935). Coordinates use GRCh38.p14/hg38 and annotation uses
RefSeq transcript NM_003165.6.
"""


def run_case(name, setup, predicate):
    with tempfile.TemporaryDirectory(prefix="ls02_accept_") as tmp:
        root = Path(tmp)
        setup(root)
        result = checker.check(root)
        return {"name": name, "passed": bool(predicate(result)), "result": result}


def main():
    cases = []
    cases.append(run_case(
        "reference_like_correct",
        lambda root: write_fixture(
            root,
            "chrom\tpos\tref\talt\tgene\tconsequence\talt_reads\ttotal_reads\tallele_fraction\n"
            "chr9\t127661125\tG\tT\tSTXBP1\tstop_gained;p.Glu117Ter\t18\t93\t0.193548\n",
            {"variant": {"chrom": "chr9", "pos": 127661125, "ref": "G", "alt": "T"},
             "read_evidence": {"alt_reads": 18, "total_reads": 93, "vaf": 18 / 93},
             "interpretation": "high-confidence mosaic; pLI/LoF-intolerant"},
            CORRECT_REPORT),
        lambda r: r["hardgate_pass"] and (r["core_science"], r["direction"], r["summary"]) == (40, 15, 5),
    ))
    cases.append(run_case(
        "empty_or_missing",
        lambda root: None,
        lambda r: not r["hardgate_pass"] and r["core_science"] == 0 and "VARIANT_TABLE_EMPTY_OR_UNREADABLE" in r["failure_codes"],
    ))
    cases.append(run_case(
        "scientifically_wrong",
        lambda root: write_fixture(
            root,
            "chrom\tpos\tref\talt\tgene\tconsequence\talt_reads\ttotal_reads\tallele_fraction\n"
            "chr9\t127661126\tA\tC\tNOTSTXBP1\tsynonymous_variant\t40\t80\t0.5\n",
            {"interpretation": "germline benign call"},
            "GRCh38 result: chr9:127661126 A>C is synonymous. Ensembl release 113."),
        lambda r: not r["hardgate_pass"] and r["core_science"] <= 2 and "TARGET_VARIANT_MISMATCH" in r["failure_codes"],
    ))
    cases.append(run_case(
        "valid_alternative_implementation",
        lambda root: write_fixture(
            root,
            "contig,coordinate,reference_allele,variant_allele,symbol,effect,variant_reads,depth,vaf\n"
            "9,127661125,G,T,STXBP1,p.E117*,18,94,19.1489%\n",
            {"call": {"chromosome": 9, "position": 127661125, "reference": "G", "alternate": "T", "gene_symbol": "STXBP1", "effect": "nonsense"},
             "pileup": {"supporting_reads": 18, "coverage": 94, "allele_fraction": 0.191489},
             "notes": ["mosaic", "haploinsufficiency/LoF intolerance"]},
            "STXBP1 chr9:127661125 G→T is a mosaic nonsense p.E117* variant; 18/94 reads (19.15%). Reference hg38; GENCODE v47 annotation."),
        lambda r: r["hardgate_pass"] and (r["core_science"], r["direction"], r["summary"]) == (40, 15, 5),
    ))
    result = {"all_passed": all(case["passed"] for case in cases), "cases": cases}
    (HERE / "acceptance-result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["all_passed"] else 1)


if __name__ == "__main__":
    main()
