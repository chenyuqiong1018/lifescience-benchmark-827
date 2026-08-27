from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

ACCEPTED = True

CANONICAL_ID = "R-HSA-6791312"
EXPECTED_P = 0.00014600657788625928  # Hypergeometric tail: N=10489, K=336, M=49, x=8.
REQUIRED_COLUMNS = ("pathway_id", "pathway_name", "overlap", "p_value", "padj", "direction")
ALIASES = {
    "pathway_id": ("pathway_id", "reactome_id", "term_id", "id"),
    "pathway_name": ("pathway_name", "pathway", "term", "term_name"),
    "overlap": ("overlap", "overlap_fraction", "gene_ratio", "genes_ratio"),
    "p_value": ("p_value", "pvalue", "p_val", "p"),
    "padj": ("padj", "adjusted_p_value", "adj_p", "q_value", "fdr"),
    "direction": ("direction", "regulation", "effect_direction"),
}


def _key(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def _read_rows(path: Path) -> list[dict[str, str]]:
    try:
        sample = path.read_text(encoding="utf-8-sig", errors="strict")
        dialect = csv.Sniffer().sniff(sample[:4096], delimiters=",\t;")
        raw = list(csv.DictReader(sample.splitlines(), dialect=dialect))
    except Exception:
        return []
    rows: list[dict[str, str]] = []
    for source in raw:
        normalized = {_key(k): ("" if v is None else str(v).strip()) for k, v in source.items() if k}
        row = {}
        for target, aliases in ALIASES.items():
            row[target] = next((normalized[a] for a in aliases if a in normalized), "")
        rows.append(row)
    return rows


def _flatten(value: object) -> str:
    if isinstance(value, dict):
        return " ".join(_flatten(v) for v in value.values())
    if isinstance(value, list):
        return " ".join(_flatten(v) for v in value)
    return str(value)


def _mechanism_match(text: str) -> bool:
    t = text.lower()
    positive = bool(re.search(r"\b(?:tp53|p53)\b", t) and re.search(r"cell[\s_-]*cycle", t))
    negated = bool(
        re.search(r"\b(?:not|no|without)\s+(?:an?\s+)?(?:tp53|p53)\b", t)
        or re.search(r"\b(?:tp53|p53)\b.{0,35}\b(?:is|was|are|were)?\s*not\s+(?:enriched|supported|the\s+mechanism)", t)
    )
    return positive and not negated


def _number(value: object) -> float | None:
    try:
        result = float(str(value).strip())
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def _overlap_is_official(value: str) -> bool:
    match = re.fullmatch(r"\s*(\d+)\s*/\s*(\d+)\s*", value)
    return bool(match and (int(match.group(1)), int(match.group(2))) == (8, 49))


def _direction_supported(text: str) -> bool:
    t = text.lower()
    supported = bool(re.search(r"\b(mixed|bidirectional|dysregulat|cell[\s_-]*cycle\s+arrest|arrest|repress|suppress|inhibit|downregulat|decreas(?:e|ed|ing))\b", t))
    proliferation_claim = bool(
        re.search(r"\b(?:increase[sd]?|promote[sd]?|activate[sd]?|upregulat\w*)\s+(?:the\s+)?(?:cell[\s_-]*cycle|proliferation)\b", t)
        or re.search(r"\b(?:cell[\s_-]*cycle|proliferation)\s+(?:is\s+)?(?:increase[sd]?|promote[sd]?|activate[sd]?|upregulat\w*)\b", t)
    )
    return supported and not proliferation_claim


def _causal_overclaim(text: str) -> bool:
    for sentence in re.split(r"[.!?\n]+", text.lower()):
        if not re.search(r"\b(?:prove[sd]?|demonstrate[sd]?|establish(?:es|ed)?|confirm(?:s|ed)?)\b", sentence):
            continue
        if re.search(r"\b(?:does|do|did|can|cannot|can't|is|was)\s+not\s+(?:prove|demonstrate|establish|confirm)|\b(?:cannot|can't)\s+(?:prove|demonstrate|establish|confirm)", sentence):
            continue
        if re.search(r"\b(?:caus\w*|mechanism|mediate[sd]?|drive[sn]?|responsible)\b", sentence):
            return True
    return False


def _causal_caveat(text: str) -> bool:
    t = text.lower()
    return bool(
        re.search(r"\b(?:does|do|did|can)\s+not\s+(?:prove|demonstrate|establish)\b", t)
        or re.search(r"\b(?:cannot|can't)\s+(?:prove|demonstrate|establish)\b", t)
        or re.search(r"\b(?:association|hypothesis|consistent\s+with|supports?\s+but|not\s+causation|not\s+causal)\b", t)
    )


def check(workspace: Path):
    out = Path(workspace) / "output"
    rows = _read_rows(out / "pathway_enrichment.csv")
    try:
        call = json.loads((out / "mechanism_call.json").read_text(encoding="utf-8"))
    except Exception:
        call = {}
    try:
        report = (out / "report.md").read_text(encoding="utf-8", errors="replace")
    except Exception:
        report = ""

    call_text = _flatten(call)
    canonical = next(
        (
            row for row in rows
            if CANONICAL_ID.lower() in (row["pathway_id"] + " " + row["pathway_name"]).lower()
            and "tp53" in row["pathway_name"].lower()
            and bool(re.search(r"cell[\s_-]*cycle", row["pathway_name"].lower()))
        ),
        None,
    )
    overlap_ok = bool(canonical and _overlap_is_official(canonical["overlap"]))
    p = _number(canonical["p_value"]) if canonical else None
    q = _number(canonical["padj"]) if canonical else None
    stats_ok = bool(
        p is not None and q is not None
        and math.isclose(p, EXPECTED_P, rel_tol=0.10, abs_tol=1e-8)
        and p <= q <= 0.05
    )
    mechanism_ok = _mechanism_match(call_text)
    direction_ok = bool(canonical and _direction_supported(canonical["direction"] + " " + call_text))
    overclaim = _causal_overclaim(report + " " + call_text)
    report_ok = bool(
        report
        and _mechanism_match(report)
        and "enrich" in report.lower()
        and _causal_caveat(report)
        and not overclaim
    )

    table_parseable = bool(rows)
    canonical_ok = canonical is not None
    core = (10 if canonical_ok else 0) + (10 if overlap_ok else 0) + (8 if stats_ok else 0) + (12 if mechanism_ok else 0)
    direction = 15 if canonical_ok and overlap_ok and direction_ok else 0
    summary = 5 if report_ok else 0

    gates = {
        "FATAL_GROUNDED_ENRICHMENT_EVIDENCE": bool(canonical_ok and overlap_ok and stats_ok),
        "FATAL_PRIMARY_MECHANISM_TRUTH": mechanism_ok,
        "FATAL_NO_CAUSAL_OVERCLAIM": not overclaim,
    }
    criteria = {
        "table_parseable": table_parseable,
        "canonical_reactome_pathway": canonical_ok,
        "official_overlap_8_of_49": overlap_ok,
        "grounded_hypergeometric_p_and_valid_fdr": stats_ok,
        "primary_tp53_cell_cycle_mechanism": mechanism_ok,
        "mixed_or_cell_cycle_repressive_direction": direction_ok,
        "report_scientifically_consistent": report_ok,
        "fatal_gates": gates,
    }
    failures = []
    if not table_parseable: failures.append("ENRICHMENT_TABLE_UNPARSEABLE")
    if not canonical_ok: failures.append("CANONICAL_PATHWAY_NOT_FOUND")
    if canonical_ok and not overlap_ok: failures.append("OFFICIAL_OVERLAP_MISMATCH")
    if canonical_ok and not stats_ok: failures.append("ENRICHMENT_STATISTICS_INVALID")
    if not mechanism_ok: failures.append("PRIMARY_MECHANISM_MISMATCH")
    if not direction_ok: failures.append("DIRECTION_UNSUPPORTED")
    if overclaim: failures.append("CAUSAL_OVERCLAIM")
    if not report_ok: failures.append("REPORT_SCIENCE_INCOMPLETE")
    return {
        "core_science": core,
        "direction": direction,
        "summary": summary,
        "hardgate_pass": all(gates.values()),
        "criteria": criteria,
        "failure_codes": failures,
    }
