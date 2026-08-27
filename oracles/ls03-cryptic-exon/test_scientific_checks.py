from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("scientific_checks", HERE / "scientific_checks.py")
SC = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(SC)


def fixture(exon: str = "", junctions: str = "", report: str = ""):
    root = Path(tempfile.mkdtemp(prefix="ls03-acceptance-"))
    out = root / "output"
    out.mkdir()
    if exon:
        (out / "cryptic_exon.tsv").write_text(exon, encoding="utf-8")
    if junctions:
        (out / "junctions.tsv").write_text(junctions, encoding="utf-8")
    if report:
        (out / "report.md").write_text(report, encoding="utf-8")
    # A candidate script is deliberately inert input: the checker must never execute it.
    (out / "analysis.py").write_text("raise RuntimeError('candidate code executed')\n", encoding="utf-8")
    return root


CORRECT_EXON = """gene\tchrom\tstart\tend\tleft_junction_reads\tright_junction_reads\texpression_evidence\tcoordinate_system
GNG10\tchr9\t111664537\t111664589\t40\t33\t510 primary gene alignments; strong split-read expression\tone-based inclusive GRCh38
"""
CORRECT_JUNCTIONS = """gene\tchrom\tintron_start\tintron_end\tjunction_reads\tnovelty\tannotation
GNG10\tchr9\t111661715\t111664536\t40\tnovel\tMANE GRCh38 v1.3
GNG10\tchr9\t111664589\t111666814\t33\tnovel\tMANE GRCh38 v1.3
"""
CORRECT_REPORT = """GNG10 is the protein-coding target. A 53 bp cryptic exon at chr9:111664537-111664589
(one-based inclusive, GRCh38) is bounded by both novel junctions, supported by 40 and 33 reads.
Novelty was assessed against MANE GRCh38 v1.3.
"""


class Acceptance(unittest.TestCase):
    results: dict[str, dict] = {}

    def record(self, name: str, result: dict, passed: bool):
        self.results[name] = {"passed": passed, "result": result}

    def test_01_reference_like_correct(self):
        r = SC.check(fixture(CORRECT_EXON, CORRECT_JUNCTIONS, CORRECT_REPORT))
        ok = r["hardgate_pass"] and (r["core_science"], r["direction"], r["summary"]) == (40, 15, 5)
        self.record("reference_like_correct", r, ok)
        self.assertTrue(ok)

    def test_02_empty_or_missing(self):
        r = SC.check(fixture())
        ok = not r["hardgate_pass"] and r["core_science"] == r["direction"] == r["summary"] == 0
        self.record("empty_or_missing", r, ok)
        self.assertTrue(ok)

    def test_03_scientifically_wrong(self):
        exon = "gene\tchrom\tstart\tend\tleft_junction_reads\tright_junction_reads\texpression_evidence\tcoordinate_system\nCD74\tchr5\t150401000\t150401052\t80\t71\thighly expressed\tone-based inclusive GRCh38\n"
        junctions = "gene\tchrom\tintron_start\tintron_end\tjunction_reads\tnovelty\tannotation\nCD74\tchr5\t150399900\t150401000\t80\tnovel\tMANE GRCh38 v1.3\nCD74\tchr5\t150401052\t150402000\t71\tnovel\tMANE GRCh38 v1.3\n"
        report = "CD74 is a protein-coding gene with a highly expressed 53 bp exon bounded by both novel junctions under MANE GRCh38 v1.3."
        r = SC.check(fixture(exon, junctions, report))
        ok = not r["hardgate_pass"] and r["core_science"] == 0
        self.record("scientifically_wrong", r, ok)
        self.assertTrue(ok)

    def test_04_plausible_but_unsupported(self):
        # Self-consistency without truth geometry; fabricated support; and negated evidence.
        self_report = SC.check(fixture(
            "gene\tchrom\tstart\tend\tleft_junction_reads\tright_junction_reads\texpression_evidence\tcoordinate_system\nABCD1\tchr9\t111700001\t111700053\t120\t95\thighly expressed\tone-based inclusive GRCh38\n",
            "gene\tchrom\tintron_start\tintron_end\tjunction_reads\tnovelty\tannotation\nABCD1\tchr9\t111699000\t111700000\t120\tnovel\tMANE GRCh38 v1.3\nABCD1\tchr9\t111700053\t111702000\t95\tnovel\tMANE GRCh38 v1.3\n",
            "ABCD1 is protein-coding; both junctions form a novel highly expressed exon under MANE GRCh38 v1.3."
        ))
        fabricated = SC.check(fixture(
            CORRECT_EXON.replace("\t40\t33\t", "\t4000\t3300\t"),
            CORRECT_JUNCTIONS.replace("\t40\tnovel", "\t4000\tnovel").replace("\t33\tnovel", "\t3300\tnovel"),
            CORRECT_REPORT.replace("40 and 33", "4000 and 3300")
        ))
        negated = SC.check(fixture(CORRECT_EXON, CORRECT_JUNCTIONS, CORRECT_REPORT.replace("both novel junctions", "both junctions are not novel")))
        ok = (not self_report["hardgate_pass"] and self_report["core_science"] == 0
              and not fabricated["hardgate_pass"] and fabricated["core_science"] < 40
              and negated["direction"] < 15 and negated["summary"] < 5)
        self.record("plausible_but_unsupported", {"self_report": self_report, "fabricated_numbers": fabricated, "negated_evidence": negated}, ok)
        self.assertTrue(ok)

    def test_05_valid_alternative_implementation(self):
        exon = "hgnc_symbol\tchromosome\texon_start_0based\texon_end_0based\tleft_reads\tright_reads\tevidence\tcoordinates\nGNG10\t9\t111664536\t111664589\t40\t33\thighly expressed; 510 primary alignments\tzero-based half-open GRCh38\n"
        junctions = "junction\tsplit_read_count\tstatus\treference_annotation\nchr9:111661715-111664536\t40\tunannotated\tMANE GRCh38 1.3\nchr9:111664589-111666814\t33\tunannotated\tMANE GRCh38 1.3\n"
        report = "Protein-coding GNG10 contains the 53 bp interval chr9:111664536-111664589 (zero-based half-open). Both splice junctions are novel with 40 and 33 split reads relative to MANE GRCh38 1.3."
        r = SC.check(fixture(exon, junctions, report))
        ok = r["hardgate_pass"] and (r["core_science"], r["direction"], r["summary"]) == (40, 15, 5)
        self.record("valid_alternative_implementation", r, ok)
        self.assertTrue(ok)


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Acceptance)
    outcome = unittest.TextTestRunner(verbosity=2).run(suite)
    payload = {
        "schema_version": 1,
        "task_id": "ls03-cryptic-exon",
        "all_passed": outcome.wasSuccessful() and len(Acceptance.results) == 5,
        "tests_run": outcome.testsRun,
        "cases": Acceptance.results,
        "candidate_code_executed": False,
    }
    (HERE / "acceptance-result.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    raise SystemExit(0 if payload["all_passed"] else 1)
