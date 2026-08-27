from __future__ import annotations

import csv
import io
import json
import math
import re
from pathlib import Path

ACCEPTED = True

POS = 127661125
ALIASES = {
    "chrom": {"chrom", "chromosome", "chr", "contig"},
    "pos": {"pos", "position", "coordinate", "genomicposition", "start"},
    "ref": {"ref", "reference", "referenceallele", "refallele"},
    "alt": {"alt", "alternate", "alternative", "alternateallele", "altallele", "variantallele"},
    "gene": {"gene", "genesymbol", "symbol", "hgncsymbol"},
    "consequence": {"consequence", "effect", "annotation", "variantconsequence", "molecularconsequence", "hgvsp", "proteinchange"},
    "alt_reads": {"altreads", "alternatereads", "variantreads", "supportingreads", "altcount", "ao"},
    "total_reads": {"totalreads", "depth", "dp", "coverage", "readdepth", "totaldepth"},
    "af": {"allelefraction", "variantallelefraction", "vaf", "af"},
}


def _norm_key(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def _field(record, name):
    normalized = {_norm_key(k): v for k, v in record.items()}
    for key in ALIASES[name]:
        if key in normalized:
            return normalized[key]
    return None


def _number(value, percent=False):
    try:
        text = str(value).strip()
        is_percent = text.endswith("%")
        number = float(text.rstrip("%"))
        if not math.isfinite(number):
            return None
        return number / 100.0 if (percent and is_percent) else number
    except (TypeError, ValueError):
        return None


def _chrom_ok(value):
    return re.sub(r"^chr", "", str(value).strip(), flags=re.I).lstrip("0") == "9"


def _allele(value):
    match = re.search(r"[ACGT]", str(value).upper())
    return match.group(0) if match else ""


def _read(path, limit=4_000_000):
    try:
        return path.read_text(encoding="utf-8", errors="replace")[:limit]
    except Exception:
        return ""


def _tsv_records(path):
    text = _read(path)
    if not text.strip():
        return []
    sample = text[:4096]
    try:
        delimiter = csv.Sniffer().sniff(sample, delimiters="\t,;").delimiter
    except csv.Error:
        delimiter = "\t"
    rows = list(csv.reader(io.StringIO(text), delimiter=delimiter))
    rows = [[cell.strip() for cell in row] for row in rows if any(cell.strip() for cell in row)]
    if not rows:
        return []
    headerish = any(_norm_key(cell) in set().union(*ALIASES.values()) for cell in rows[0])
    if headerish:
        header = rows[0]
        return [dict(zip(header, row)) for row in rows[1:]]
    names = ["chrom", "pos", "ref", "alt", "gene", "consequence", "alt_reads", "total_reads", "allele_fraction"]
    return [dict(zip(names, row)) for row in rows]


def _walk_dicts(value):
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from _walk_dicts(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk_dicts(child)


def _flatten(value, out=None):
    out = {} if out is None else out
    if isinstance(value, dict):
        for key, child in value.items():
            if isinstance(child, (dict, list)):
                _flatten(child, out)
            else:
                out[key] = child
    elif isinstance(value, list):
        for child in value:
            _flatten(child, out)
    return out


def _json_records(path):
    try:
        obj = json.loads(_read(path))
    except Exception:
        return []
    records = list(_walk_dicts(obj))
    flat = _flatten(obj)
    if flat:
        records.append(flat)
    return records


def _identity(record):
    pos = _number(_field(record, "pos"))
    return {
        "chrom": _chrom_ok(_field(record, "chrom")),
        "pos": pos is not None and int(pos) == POS,
        "ref": _allele(_field(record, "ref")) == "G",
        "alt": _allele(_field(record, "alt")) == "T",
    }


def _target_score(record):
    checks = _identity(record)
    return 2 * checks["pos"] + checks["chrom"] + checks["ref"] + checks["alt"]


def _consequence_ok(value):
    text = str(value).lower()
    return bool(re.search(r"stop.?gained|nonsense|glu117(?:ter|\*)|e117\*|p\.e117x", text))


def _quant(record):
    alt = _number(_field(record, "alt_reads"))
    total = _number(_field(record, "total_reads"))
    af = _number(_field(record, "af"), percent=True)
    count_ok = alt is not None and total is not None and 17 <= alt <= 19 and 90 <= total <= 96
    af_ok = af is not None and 0.175 <= af <= 0.205
    consistent = alt is not None and total not in (None, 0) and af is not None and abs(af - alt / total) <= 0.012
    return count_ok, af_ok, consistent, alt, total, af


def check(workspace: Path):
    output = Path(workspace) / "output"
    tsv = _tsv_records(output / "variant.tsv")
    evidence_records = _json_records(output / "evidence.json")
    report = _read(output / "report.md", 1_000_000)
    evidence_text = _read(output / "evidence.json", 1_000_000)

    target = max(tsv, key=_target_score) if tsv else {}
    ident = _identity(target)
    gene_ok = str(_field(target, "gene")).strip().upper() == "STXBP1"
    consequence_ok = _consequence_ok(_field(target, "consequence"))

    tq = _quant(target)
    evidence_target = [r for r in evidence_records if all(_identity(r).values())]
    eq = [_quant(r) for r in evidence_target]
    target_has_quant = any(v is not None for v in tq[3:])
    quantitative_ok = all(tq[:3]) if target_has_quant else any(all(q[:3]) for q in eq)

    core = (2 if ident["chrom"] else 0) + (12 if ident["pos"] else 0)
    core += (5 if ident["ref"] else 0) + (5 if ident["alt"] else 0)
    core += (8 if gene_ok else 0) + (5 if consequence_ok else 0) + (3 if quantitative_ok else 0)

    corpus = (report + "\n" + evidence_text).lower()
    call_ok = all(ident.values()) and gene_ok and consequence_ok
    mosaic_ok = quantitative_ok and "mosaic" in corpus
    intolerance_ok = bool(re.search(r"loss[- ]of[- ]function.{0,35}intoler|lof[- ]?intoler|\bpli\b|\bloeuf\b|haploinsuff", corpus, re.S))
    direction = (8 if call_ok else 0) + (4 if call_ok and mosaic_ok else 0) + (3 if call_ok and intolerance_ok else 0)

    low_report = report.lower()
    report_call = bool("stxbp1" in low_report and str(POS) in low_report and
                       (re.search(r"g\s*>\s*t", low_report) or "g→t" in low_report) and
                       re.search(r"stop.?gained|nonsense|glu117(?:ter|\*)|e117\*", low_report))
    report_quant = bool(re.search(r"18\s*/\s*9[34]", low_report) or
                        re.search(r"\b0\.19\d*\b|\b19(?:\.1|\.2|\.3)?\s*%", low_report))
    reference_version = bool(re.search(r"\b(?:grch38|hg38)(?:\.p\d+)?\b", low_report))
    annotation_version = bool(re.search(r"\b(?:gencode\s*(?:v|release\s*)\d+|ensembl\s*(?:release\s*)?\d+|refseq[^\n]{0,35}(?:release\s*\d+|annotation\s*\d+|nm_\d+\.\d+))", low_report))
    summary = (3 if report_call and report_quant else 0) + (1 if reference_version else 0) + (1 if annotation_version else 0)

    hardgate = call_ok and quantitative_ok
    criteria = {
        "chromosome_9": ident["chrom"], "position_127661125": ident["pos"],
        "ref_G": ident["ref"], "alt_T": ident["alt"], "gene_STXBP1": gene_ok,
        "nonsense_stop_gained": consequence_ok, "read_counts_and_vaf": quantitative_ok,
        "mosaic_interpretation": mosaic_ok, "lof_intolerance_supported": intolerance_ok,
        "report_consistent": report_call and report_quant,
        "reference_version_documented": reference_version,
        "annotation_version_documented": annotation_version,
    }
    failures = []
    if not tsv:
        failures.append("VARIANT_TABLE_EMPTY_OR_UNREADABLE")
    if not all(ident.values()):
        failures.append("TARGET_VARIANT_MISMATCH")
    if not gene_ok:
        failures.append("GENE_MISMATCH")
    if not consequence_ok:
        failures.append("CONSEQUENCE_MISMATCH")
    if not quantitative_ok:
        failures.append("READ_EVIDENCE_MISMATCH")
    if not mosaic_ok:
        failures.append("MOSAIC_INTERPRETATION_MISSING")
    if not intolerance_ok:
        failures.append("LOF_INTOLERANCE_SUPPORT_MISSING")
    if not (report_call and report_quant):
        failures.append("REPORT_MISMATCH")
    if not (reference_version and annotation_version):
        failures.append("VERSION_DOCUMENTATION_MISSING")
    return {"core_science": min(40, core), "direction": min(15, direction),
            "summary": min(5, summary), "hardgate_pass": bool(hardgate),
            "criteria": criteria, "failure_codes": failures}
