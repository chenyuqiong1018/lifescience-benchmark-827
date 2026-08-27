from __future__ import annotations

import csv
import math
import re
from pathlib import Path

ACCEPTED = True

# Independent truth source: active-set non-negative least squares of the raw
# Spot_710-1 expression vector against the six cell-type mean profiles in the
# authorized spatial.sim.tar.gz (archive SHA-256 f290570ad3d3230c0672e3050b54b37f9d89e5673f2a10216ccd359233aa417a).
# The six coefficients below are the NNLS coefficients normalized to sum to one.
EXPECTED = {
    "B_Cell": 0.323040493130,
    "Endothelial": 0.346842122837,
    "Fibroblast_Stroma": 0.0,
    "Macrophage": 0.318705728025,
    "T_Cell": 0.0114116560073,
    "Tumor_Core": 0.0,
}
MAJOR = ("B_Cell", "Endothelial", "Macrophage")
OTHER = ("Fibroblast_Stroma", "T_Cell", "Tumor_Core")

ALIASES = {
    "B_Cell": (r"\bb[-_ ]?cells?\b", r"\bb[-_ ]?lymphocytes?\b"),
    "Endothelial": (r"\bendothelial(?:[-_ ]?cells?)?\b",),
    "Fibroblast_Stroma": (r"\bfibroblasts?\b", r"\bstroma(?:l)?(?:[-_ ]?cells?)?\b"),
    "Macrophage": (r"\bmacrophages?\b",),
    "T_Cell": (r"\bt[-_ ]?cells?\b", r"\bt[-_ ]?lymphocytes?\b"),
    "Tumor_Core": (r"\btumou?r(?:[-_ ]?core|[-_ ]?cells?)?\b", r"\bmalignant(?:[-_ ]?cells?)?\b"),
}


def _canon(value: object) -> str | None:
    s = re.sub(r"[^a-z0-9]+", "", str(value).lower())
    if s.startswith(("bcell", "blymph")):
        return "B_Cell"
    if s.startswith("endothelial"):
        return "Endothelial"
    if s.startswith(("fibro", "stroma")):
        return "Fibroblast_Stroma"
    if s.startswith("macroph"):
        return "Macrophage"
    if s.startswith(("tcell", "tlymph")):
        return "T_Cell"
    if s.startswith(("tumor", "tumour", "malignant", "cancer")):
        return "Tumor_Core"
    return None


def _header(fieldnames: list[str] | None, choices: set[str]) -> str | None:
    for name in fieldnames or []:
        if re.sub(r"[^a-z0-9]+", "", name.lower()) in choices:
            return name
    return None


def _composition(path: Path):
    raw = {k: 0.0 for k in EXPECTED}
    unknown = 0.0
    problems: list[str] = []
    rows = 0
    try:
        text = path.read_text(encoding="utf-8-sig", errors="replace")
        dialect = csv.Sniffer().sniff(text[:4096], delimiters=",;\t")
        reader = csv.DictReader(text.splitlines(), dialect=dialect)
        type_col = _header(reader.fieldnames, {"celltype", "type", "cellidentity", "label"})
        weight_col = _header(reader.fieldnames, {"weight", "proportion", "fraction", "composition", "abundance"})
        if not type_col or not weight_col:
            raise ValueError("required semantic columns not found")
        for row in reader:
            if not any(str(v).strip() for v in row.values() if v is not None):
                continue
            rows += 1
            try:
                value = float(row.get(weight_col, ""))
            except (TypeError, ValueError):
                problems.append("non-numeric weight")
                continue
            if not math.isfinite(value) or value < 0:
                problems.append("non-finite or negative weight")
                continue
            label = _canon(row.get(type_col, ""))
            if label is None:
                unknown += value
            else:
                raw[label] += value
    except Exception as exc:
        problems.append(type(exc).__name__)

    total = sum(raw.values()) + unknown
    usable = rows > 0 and not problems and math.isfinite(total) and total > 0
    vector = {k: (v / total if usable else 0.0) for k, v in raw.items()}
    unknown_fraction = unknown / total if usable else 0.0
    normalized = usable and abs(total - 1.0) <= 0.010000001
    recognized = usable and unknown_fraction <= 0.01
    return vector, total, unknown_fraction, usable, normalized and recognized, problems


def _mention_state(text: str, label: str) -> tuple[bool, bool]:
    positive = negative = False
    for pattern in ALIASES[label]:
        for match in re.finditer(pattern, text, flags=re.I):
            before = text[max(0, match.start() - 40):match.start()].lower()
            after = text[match.end():match.end() + 35].lower()
            pre_neg = bool(re.search(r"(?:\bno\b|\bnot\b|\bwithout\b|\babsent\b|\bmissing\b|\black(?:s|ed|ing)?\b)(?:(?![.!?;]).){0,28}$", before))
            if re.search(r"\bnot\s+only\s*$", before):
                pre_neg = False
            post_neg = bool(re.match(r"[ \t]*(?:(?:is|are|was|were|seems?|appears?)[ \t]*)?(?:not\b|absent\b|missing\b|excluded\b|unsupported\b|negligible\b|zero\b)", after))
            negative |= pre_neg or post_neg
            positive |= not (pre_neg or post_neg)
    return positive, negative


def _reported_number_conflict(text: str, label: str, submitted: float) -> bool:
    number = r"(?P<n>(?:0(?:\.\d+)?|1(?:\.0+)?|\d{1,3}(?:\.\d+)?)\s*%?)"
    for alias in ALIASES[label]:
        patterns = (alias + r"[^\n.!?;]{0,18}?" + number, number + r"[^\n.!?;]{0,18}?" + alias)
        for pattern in patterns:
            for match in re.finditer(pattern, text, flags=re.I):
                token = match.group("n").replace(" ", "")
                value = float(token.rstrip("%")) / (100.0 if token.endswith("%") else 1.0)
                if 0 <= value <= 1 and abs(value - submitted) > 0.15:
                    return True
    return False


def check(workspace: Path):
    vector, total, unknown, usable, valid, problems = _composition(
        Path(workspace) / "output" / "spot_710_composition.csv"
    )
    background = sum(vector[k] for k in OTHER) + unknown
    major_total = sum(vector[k] for k in MAJOR)
    l1 = sum(abs(vector[k] - EXPECTED[k]) for k in EXPECTED) + unknown if usable else 2.0

    # Exactly two explicitly named fatal scientific gates.
    gate_valid = valid
    gate_major = (
        valid
        and all(vector[k] >= 0.18 for k in MAJOR)
        and major_total >= 0.80
        and background <= 0.20
    )
    hardgate = gate_valid and gate_major

    structure_points = 8.0 if valid else 0.0
    major_points = sum(6.0 * max(0.0, 1.0 - abs(vector[k] - EXPECTED[k]) / 0.25) for k in MAJOR) if usable else 0.0
    global_points = 8.0 * max(0.0, 1.0 - l1 / 0.80) if usable else 0.0
    expected_background = sum(EXPECTED[k] for k in OTHER)
    background_points = 6.0 * max(0.0, 1.0 - abs(background - expected_background) / 0.30) if usable else 0.0
    core = max(0, min(40, round(structure_points + major_points + global_points + background_points)))

    recovered = sum(2 for k in MAJOR if vector[k] >= 0.15) if usable else 0
    balanced = 4 if usable and all(0.20 <= vector[k] <= 0.50 for k in MAJOR) else 0
    mixed = 2 if usable and sum(v >= 0.20 for v in vector.values()) >= 2 and max(vector.values()) <= 0.70 else 0
    clean_background = 3 if usable and background <= 0.10 else 0
    direction = recovered + balanced + mixed + clean_background

    report_path = Path(workspace) / "output" / "report.md"
    report = report_path.read_text(encoding="utf-8", errors="replace") if report_path.is_file() else ""
    states = {k: _mention_state(report, k) for k in EXPECTED}
    mixture_stated = bool(re.search(r"\b(?:mix(?:ed|ture)?|admix(?:ed|ture)?|combination|compos(?:ed|ition)|contribut(?:e|es|ion))\b", report, re.I))
    types_supported = all(states[k][0] and not states[k][1] for k in MAJOR)
    unsupported_major_claim = any(
        states[k][0]
        and any(re.search(r"\b(?:major|dominant|substantial|primary|main)\b", report[max(0, m.start()-35):m.end()+35], re.I)
                for alias in ALIASES[k] for m in re.finditer(alias, report, re.I))
        for k in OTHER
    )
    numbers_consistent = not any(_reported_number_conflict(report, k, vector[k]) for k in MAJOR)
    summary = 0
    if hardgate:
        summary += 2 if mixture_stated else 0
        summary += 3 if types_supported and not unsupported_major_claim and numbers_consistent else 0

    criteria = {
        "fatal_gate_valid_normalized_composition": gate_valid,
        "fatal_gate_supported_major_mixture": gate_major,
        "composition_usable": usable,
        "weight_sum": total,
        "unknown_weight_fraction": unknown,
        "submitted_normalized_weights": vector,
        "truth_weights": EXPECTED,
        "l1_distance_to_truth": l1,
        "major_type_points": round(major_points, 4),
        "global_fidelity_points": round(global_points, 4),
        "background_points": round(background_points, 4),
        "direction_major_recovery": recovered,
        "direction_balanced_three_way_mixture": balanced,
        "direction_not_forced_single_type": mixed,
        "direction_unsupported_types_minor": clean_background,
        "report_mixture_stated": mixture_stated,
        "report_expected_types_affirmed": types_supported,
        "report_no_unsupported_major_claim": not unsupported_major_claim,
        "report_numbers_consistent_with_submission": numbers_consistent,
        "parse_problems": problems,
    }
    failures: list[str] = []
    if not gate_valid:
        failures.append("FATAL_INVALID_NORMALIZED_COMPOSITION")
    if not gate_major:
        failures.append("FATAL_SUPPORTED_MAJOR_MIXTURE_NOT_RECOVERED")
    if usable and l1 > 0.30:
        failures.append("QUANTITATIVE_COMPOSITION_MISMATCH")
    if direction < 15:
        failures.append("MIXTURE_DIRECTION_INCOMPLETE")
    if summary < 5:
        failures.append("REPORT_SUMMARY_NOT_SCIENTIFICALLY_SUPPORTED")
    return {
        "core_science": core,
        "direction": direction,
        "summary": summary,
        "hardgate_pass": hardgate,
        "criteria": criteria,
        "failure_codes": failures,
    }
