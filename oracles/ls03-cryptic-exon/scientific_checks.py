from __future__ import annotations

import csv
import re
from pathlib import Path

ACCEPTED = True

# Immutable truth from grounding-manifest.json (Genentech/compbiobench-data-v1,
# revision c673f0855fce09d320f1677f168f7864eec52c1a; verified artifact hashes).
GENE, CHROM = "GNG10", "chr9"
EXON_0 = (111664536, 111664589)
EXON_1 = (111664537, 111664589)
LEFT = (111661715, 111664536, 40)
RIGHT = (111664589, 111666814, 33)


def _read(path: Path, limit: int = 2_000_000) -> str:
    try:
        if not path.is_file() or path.stat().st_size > limit:
            return ""
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""


def _rows(path: Path) -> list[dict[str, str]]:
    text = _read(path)
    if not text.strip():
        return []
    try:
        reader = csv.DictReader(text.splitlines(), delimiter="\t")
        return [
            {str(k).strip().lower(): str(v or "").strip() for k, v in row.items() if k is not None}
            for row in reader
        ]
    except (csv.Error, TypeError):
        return []


def _value(row: dict[str, str], *names: str) -> str:
    for name in names:
        if name in row and row[name] != "":
            return row[name]
    return ""


def _integer(value: object) -> int | None:
    match = re.search(r"(?<![\d.])-?\d+(?![\d.])", str(value).replace(",", ""))
    return int(match.group()) if match else None


def _chrom(value: str) -> str:
    value = value.strip().lower()
    return "chr" + value[3:] if value.startswith("chr") else "chr" + value


def _explicit_convention(row: dict[str, str], report: str) -> str:
    keys = " ".join(row)
    values = " ".join(row.values()).lower()
    text = f"{keys} {values} {report.lower()}"
    if re.search(r"(?:0|zero)[- _]?based", text) and re.search(r"half[- _]?open|exclusive", text):
        return "zero_half_open"
    if re.search(r"(?:1|one)[- _]?based", text) and re.search(r"inclusive", text):
        return "one_inclusive"
    return ""


def _exon_row(rows: list[dict[str, str]], report: str) -> tuple[dict[str, str], bool, bool, bool]:
    for row in rows:
        gene = _value(row, "gene", "hgnc_gene", "hgnc_symbol", "gene_symbol", "symbol").upper()
        chrom = _chrom(_value(row, "chrom", "chromosome", "chr"))
        start = _integer(_value(row, "start", "exon_start", "exon_start_0based", "start_0based", "start_1based"))
        end = _integer(_value(row, "end", "exon_end", "exon_end_0based", "end_0based", "end_1based"))
        convention = _explicit_convention(row, report)
        interval_ok = (start, end) == EXON_0 and convention == "zero_half_open"
        interval_ok |= (start, end) == EXON_1 and convention == "one_inclusive"
        if gene == GENE or (chrom == CHROM and interval_ok):
            return row, gene == GENE, chrom == CHROM, interval_ok
    return {}, False, False, False


def _junction_tuple(row: dict[str, str]) -> tuple[str, int | None, int | None, int | None]:
    chrom = _value(row, "chrom", "chromosome", "chr")
    start = _integer(_value(row, "intron_start", "junction_start", "donor", "start", "left"))
    end = _integer(_value(row, "intron_end", "junction_end", "acceptor", "end", "right"))
    packed = _value(row, "junction", "junction_id", "coordinates", "coord", "locus")
    if packed:
        match = re.search(r"(?:(chr)?([0-9xy]+)[:_])?(\d{6,})\s*[-:]\s*(\d{6,})", packed, re.I)
        if match:
            chrom = chrom or ((match.group(1) or "") + (match.group(2) or ""))
            start, end = int(match.group(3)), int(match.group(4))
    reads = _integer(_value(row, "junction_reads", "split_read_count", "split_reads", "read_count", "reads", "count", "support"))
    return _chrom(chrom), start, end, reads


def _find_junction(rows: list[dict[str, str]], expected: tuple[int, int, int]) -> tuple[bool, bool, dict[str, str]]:
    for row in rows:
        chrom, start, end, reads = _junction_tuple(row)
        gene = _value(row, "gene", "hgnc_gene", "hgnc_symbol", "gene_symbol", "symbol").upper()
        gene_ok = not gene or gene == GENE
        if gene_ok and chrom == CHROM and (start, end) == expected[:2]:
            return True, reads == expected[2], row
    return False, False, {}


def _novel(row: dict[str, str]) -> bool:
    value = _value(row, "novelty", "status", "is_novel", "novel").strip().lower()
    return value in {"novel", "true", "yes", "1", "unannotated", "not_annotated"}


def check(workspace: Path):
    output = Path(workspace) / "output"
    report = _read(output / "report.md")
    report_l = report.lower()
    exon_rows = _rows(output / "cryptic_exon.tsv")
    junction_rows = _rows(output / "junctions.tsv")

    exon, gene_ok, chrom_ok, interval_ok = _exon_row(exon_rows, report)
    left_geom, left_count, left_row = _find_junction(junction_rows, LEFT)
    right_geom, right_count, right_row = _find_junction(junction_rows, RIGHT)
    start = _integer(_value(exon, "start", "exon_start", "exon_start_0based", "start_0based", "start_1based"))
    end = _integer(_value(exon, "end", "exon_end", "exon_end_0based", "end_0based", "end_1based"))
    convention = _explicit_convention(exon, report)
    length_ok = interval_ok and ((end - start == 53) if convention == "zero_half_open" else (end - start + 1 == 53))
    exon_counts = (
        _integer(_value(exon, "left_junction_reads", "left_reads", "left_support")) == LEFT[2]
        and _integer(_value(exon, "right_junction_reads", "right_reads", "right_support")) == RIGHT[2]
    )
    expression = _value(exon, "expression_evidence", "evidence", "expression", "supporting_evidence").lower()
    expression_ok = interval_ok and exon_counts and bool(re.search(r"\b(?:510|high(?:ly)?|express|junction|split)\b", expression))

    annotation_text = " ".join([report_l] + [" ".join(r.values()).lower() for r in junction_rows])
    annotation_ok = "mane" in annotation_text and bool(re.search(r"v?1\.3\b", annotation_text)) and "grch38" in annotation_text
    negated_novelty = bool(re.search(r"\b(?:not|non)[- ]+novel\b|\bpreviously[- ]+annotated\b|\bis[- ]+annotated\b", annotation_text))
    novelty_ok = left_geom and right_geom and (_novel(left_row) and _novel(right_row) or "both" in report_l and "novel" in report_l) and not negated_novelty
    protein_coding_ok = gene_ok and bool(re.search(r"protein[- ]coding", report_l))

    criteria = {
        "truth_gene_gng10": gene_ok,
        "truth_chromosome_chr9": chrom_ok,
        "truth_exon_interval_with_explicit_convention": interval_ok,
        "truth_exon_length_53bp": length_ok,
        "truth_left_junction_geometry": left_geom,
        "truth_left_junction_40_reads": left_count,
        "truth_right_junction_geometry": right_geom,
        "truth_right_junction_33_reads": right_count,
        "truth_expression_evidence": expression_ok,
        "mane_grch38_v1_3_provenance": annotation_ok,
        "both_junctions_correctly_called_novel": novelty_ok,
        "protein_coding_target_conclusion": protein_coding_ok,
    }
    core = sum((8 * gene_ok, 8 * (chrom_ok and interval_ok), 4 * length_ok,
                5 * left_geom, 4 * left_count, 5 * right_geom, 4 * right_count, 2 * expression_ok))
    direction = 5 * annotation_ok + 6 * novelty_ok + 4 * protein_coding_ok

    summary_facts = {
        "report_gene": bool(re.search(r"\bgng10\b", report_l)),
        "report_53bp_interval": "53" in report_l and ("111664536" in report_l or "111664537" in report_l) and "111664589" in report_l,
        "report_support_counts": bool(re.search(r"\b40\b", report_l) and re.search(r"\b33\b", report_l)),
        "report_novelty_and_annotation": "novel" in report_l and "mane" in report_l and "1.3" in report_l and not negated_novelty,
    }
    criteria.update(summary_facts)
    summary = sum((2 * summary_facts["report_gene"], summary_facts["report_53bp_interval"],
                   summary_facts["report_support_counts"], summary_facts["report_novelty_and_annotation"]))

    # Exactly three fatal scientific gates; all other criteria retain partial credit.
    fatal_gates = {
        "FATAL_TRUTH_GENE": gene_ok,
        "FATAL_TRUTH_EVENT_GEOMETRY": chrom_ok and interval_ok and length_ok,
        "FATAL_TWO_JUNCTION_SUPPORT": left_geom and left_count and right_geom and right_count,
    }
    failures = [name for name, passed in fatal_gates.items() if not passed]
    if direction < 15:
        failures.append("DIRECTION_INCOMPLETE")
    if summary < 5:
        failures.append("SUMMARY_INCOMPLETE")
    criteria["fatal_gates"] = fatal_gates
    return {
        "core_science": int(core),
        "direction": int(direction),
        "summary": int(summary),
        "hardgate_pass": all(fatal_gates.values()),
        "criteria": criteria,
        "failure_codes": failures,
    }
