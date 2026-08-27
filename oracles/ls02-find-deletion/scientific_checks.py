from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

ACCEPTED = True

# Immutable truth from grounding-manifest.json (input hashes verified by the author).
TRUTH = {
    "chrom": "22",
    "start": 20_000_000,
    "end": 21_000_000,
    "size": 1_000_000,
    "depth_ratio": 0.4922519150926846,
}
TOLERANCE_BP = 100_000
FATAL_GATES = ("FATAL_GROUNDED_LOCUS", "FATAL_GROUNDED_DEPTH_EVIDENCE")


def _norm(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def _number(value: object) -> float:
    text = str(value).strip().lower().replace(",", "")
    match = re.search(r"[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?", text)
    if not match:
        return math.nan
    number = float(match.group())
    if re.search(r"\bmb\b", text):
        number *= 1_000_000
    elif re.search(r"\bkb\b", text):
        number *= 1_000
    return number


def _read_text(path: Path, limit: int = 1_000_000) -> str:
    try:
        if not path.is_file() or path.stat().st_size > limit:
            return ""
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""


ALIASES = {
    "chrom": {"chrom", "chromosome", "chr", "contig"},
    "start": {"start_100kb", "start", "start_bp", "breakpoint_start", "deletion_start"},
    "end": {"end_100kb", "end", "end_bp", "breakpoint_end", "deletion_end"},
    "size": {"size_bp", "size", "length", "length_bp", "deletion_size"},
    "support": {"supporting_signals", "supporting_signal", "support", "signals", "evidence"},
}


def _deletion_rows(path: Path) -> list[dict[str, str]]:
    text = _read_text(path)
    if not text.strip():
        return []
    try:
        dialect = csv.Sniffer().sniff(text[:4096], delimiters="\t,;")
        reader = csv.DictReader(text.splitlines(), dialect=dialect)
        headers = {_norm(h): h for h in (reader.fieldnames or []) if h}
        selected = {key: next((headers[a] for a in names if a in headers), None)
                    for key, names in ALIASES.items()}
        if not all(selected[key] for key in ("chrom", "start", "end", "size")):
            return []
        rows = []
        for raw in reader:
            row = {key: str(raw.get(header, "") or "").strip()
                   for key, header in selected.items() if header}
            if any(row.get(key) for key in ("chrom", "start", "end", "size")):
                rows.append(row)
        return rows
    except (csv.Error, TypeError):
        return []


def _chrom(value: object) -> str:
    return re.sub(r"^chr", "", str(value).strip().lower())


def _close(value: object, expected: float, tolerance: float = TOLERANCE_BP) -> bool:
    number = _number(value)
    return math.isfinite(number) and abs(number - expected) <= tolerance


def _json(path: Path) -> object:
    try:
        if not path.is_file() or path.stat().st_size > 1_000_000:
            return {}
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def _candidate_ratios(obj: object) -> list[float]:
    found: list[float] = []
    if isinstance(obj, dict):
        for key, value in obj.items():
            name = _norm(key)
            if "ratio" in name and any(token in name for token in ("depth", "coverage", "flank")):
                try:
                    number = float(value)
                    if 1 < number <= 100:
                        number /= 100
                    found.append(number)
                except (TypeError, ValueError):
                    found.append(math.nan)
            elif isinstance(value, (dict, list)):
                found.extend(_candidate_ratios(value))
    elif isinstance(obj, list):
        for value in obj:
            found.extend(_candidate_ratios(value))
    return found


def _text_ratios(text: str) -> list[float]:
    values = []
    patterns = (
        r"(?:depth|coverage)(?:[-_ ]to[-_ ]flank)?[-_ ]ratio\s*(?:is|=|:|of)?\s*([0-9]+(?:\.[0-9]+)?%?)",
        r"ratio\s+(?:of\s+)?([0-9]+(?:\.[0-9]+)?%?)\s*(?:for\s+)?(?:depth|coverage)",
    )
    for pattern in patterns:
        for raw in re.findall(pattern, text, flags=re.I):
            number = float(raw.rstrip("%"))
            if raw.endswith("%") or 1 < number <= 100:
                number /= 100
            values.append(number)
    return values


def _evidence_state(text: str, qc: object) -> tuple[bool, bool, bool]:
    lowered = re.sub(r"\s+", " ", text.lower())
    positive = bool(
        re.search(r"(?:depth|coverage).{0,70}(?:deplet|drop|reduc|decreas|deficit|loss|low|lower|half)", lowered)
        or re.search(r"(?:deplet|drop|reduc|decreas|deficit|loss|low|lower|half).{0,70}(?:depth|coverage)", lowered)
    )
    contradicted = bool(
        re.search(r"(?:no|not|without|absent|lacking).{0,30}(?:depth|coverage).{0,35}(?:deplet|drop|reduc|loss|support|evidence|change)", lowered)
        or re.search(r"(?:no|not|without|absent|lacking).{0,55}(?:deplet|drop|reduc|decreas|deficit|loss|low|lower).{0,25}(?:depth|coverage)", lowered)
        or re.search(r"(?:depth|coverage).{0,30}(?:normal|unchanged|not reduced|not depleted|no drop)", lowered)
    )
    ratios = _candidate_ratios(qc) + _text_ratios(lowered)
    numeric_ok = all(math.isfinite(v) and abs(v - TRUTH["depth_ratio"]) <= 0.20 for v in ratios)
    return positive, contradicted, numeric_ok


def check(workspace: Path):
    output = Path(workspace) / "output"
    rows = _deletion_rows(output / "deletion.tsv")
    singleton = len(rows) == 1
    row = rows[0] if singleton else {}

    chrom_ok = singleton and _chrom(row.get("chrom", "")) == TRUTH["chrom"]
    start_ok = singleton and _close(row.get("start", ""), TRUTH["start"])
    end_ok = singleton and _close(row.get("end", ""), TRUTH["end"])
    size_ok = singleton and _close(row.get("size", ""), TRUTH["size"])
    grounded_locus = chrom_ok and start_ok and end_ok and size_ok

    report = _read_text(output / "report.md")
    combined = " ".join((row.get("support", ""), report))
    qc = _json(output / "qc.json")
    positive, contradicted, numeric_ok = _evidence_state(combined, qc)
    manifest_depth_support = TRUTH["depth_ratio"] < 0.65
    grounded_evidence = grounded_locus and manifest_depth_support and positive and not contradicted and numeric_ok

    deletion_claim = bool(re.search(
        r"\bdelet(?:ion|ed)\b|\bheterozygous\s+(?:cnv\s+)?loss\b|\b(?:copy[- ]number|cnv)\s+loss\b",
        combined, flags=re.I,
    ))
    duplication_claim = bool(re.search(r"\bduplication\b|\bcopy[- ]number\s+gain\b", combined, flags=re.I))
    direction_class = grounded_evidence and deletion_claim and not duplication_claim

    precision_ok = bool(re.search(
        r"(?:nearest|rounded|resolution|precision|approx(?:imate|imately)?).{0,45}100\s*[- ]?kb|"
        r"100\s*[- ]?kb.{0,45}(?:nearest|rounded|resolution|precision|approx(?:imate|imately)?)",
        report, flags=re.I,
    ))
    exact_precision_phrase = bool(re.search(
        r"(?:exact|base[- ]pair|single[- ]base).{0,35}(?:breakpoint|precision|resolution)", report, flags=re.I,
    ))
    exact_precision_negated = bool(re.search(
        r"(?:no|not|cannot|can't|do not|does not|without|insufficient).{0,45}(?:exact|base[- ]pair|single[- ]base)",
        report, flags=re.I,
    ))
    false_precision = exact_precision_phrase and not exact_precision_negated
    precision_ok = precision_ok and not false_precision

    core = (8 if chrom_ok else 0) + (8 if start_ok else 0) + (8 if end_ok else 0)
    core += (6 if size_ok else 0) + (10 if grounded_evidence else 0)
    direction = (8 if direction_class else 0) + (7 if grounded_locus and precision_ok else 0)

    compact_report = report.lower().replace(",", "")
    report_locus = "chr22" in compact_report and bool(
        re.search(r"20(?:\.0)?\s*(?:[-–—]|to)\s*21(?:\.0)?\s*mb", compact_report)
        or ("20000000" in compact_report and "21000000" in compact_report)
    )
    report_size = bool(re.search(r"\b1(?:\.0)?\s*mb\b|\b1000000\s*bp\b", compact_report))
    report_consistent = report_locus and report_size and direction_class and precision_ok
    summary = 5 if report_consistent else 0

    criteria = {
        "singleton_deletion_record": singleton,
        "chromosome_matches_ground_truth": chrom_ok,
        "start_within_100kb": start_ok,
        "end_within_100kb": end_ok,
        "size_within_100kb": size_ok,
        "depth_evidence_positive": positive,
        "depth_evidence_not_negated": not contradicted,
        "reported_depth_ratios_concordant": numeric_ok,
        "deletion_direction": direction_class,
        "precision_limit_100kb": precision_ok,
        "report_consistent": report_consistent,
        "fatal_gates": {
            FATAL_GATES[0]: grounded_locus,
            FATAL_GATES[1]: grounded_evidence,
        },
    }
    failures = []
    if not grounded_locus:
        failures.append(FATAL_GATES[0])
    if not grounded_evidence:
        failures.append(FATAL_GATES[1])
    if not direction_class:
        failures.append("DELETION_DIRECTION_UNSUPPORTED")
    if not precision_ok:
        failures.append("PRECISION_LIMIT_NOT_ESTABLISHED")
    if not report_consistent:
        failures.append("REPORT_SCIENCE_MISMATCH")

    return {
        "core_science": core,
        "direction": direction,
        "summary": summary,
        "hardgate_pass": grounded_locus and grounded_evidence,
        "criteria": criteria,
        "failure_codes": failures,
    }
