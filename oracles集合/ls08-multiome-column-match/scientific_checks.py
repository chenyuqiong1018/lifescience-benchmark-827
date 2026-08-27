from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

ACCEPTED = True

# Frozen, provenance-bearing truth in grounding-manifest.json.  The digest makes
# accidental or candidate-side replacement of the truth source fail closed.
GROUNDING_SHA256 = "8abe24df451ce434314b2b09128358c203006a8668a12168d8022f2f14f4d6c9"
INPUT_SHA256 = {
    "multiome.match.atac.rna.q1.atac.tsv.gz": "2b72f84d4a9c75dc78611b3b6ff727a0d2190e3d98f49f2fd08e8d1ebfe00d16",
    "multiome.match.atac.rna.q1.rna.tsv.gz": "caceed54a983fb7c0b58cb5d34a61a999b0e7dec9efa2dd0175450fa0f080784",
}
EXPECTED = {0: 5, 1: 1, 2: 4, 3: 0, 4: 6, 5: 3, 6: 7, 7: 2}
FATAL_GATES = ("grounding_integrity", "exact_bijective_mapping")


def _grounded() -> bool:
    path = Path(__file__).with_name("grounding-manifest.json")
    try:
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != GROUNDING_SHA256:
            return False
        manifest = json.loads(raw)
        integrity = {x["name"]: x["sha256"] for x in manifest["authorized_input_integrity"]}
        truth = manifest["recomputed_truth"]["rna_to_atac_mapping"]
        rule = manifest["frozen_matching_rule"]
        return (
            manifest.get("task_id") == "ls08-multiome-column-match"
            and manifest.get("source", {}).get("dataset_revision")
            == "c673f0855fce09d320f1677f168f7864eec52c1a"
            and integrity == INPUT_SHA256
            and truth == [EXPECTED[i] for i in range(8)]
            and all(k in rule for k in ("gene_activity", "transform", "similarity", "assignment"))
        )
    except Exception:
        return False


def _index(value: object, modality: str) -> int | None:
    text = str(value).strip()
    prefixes = (
        r"(?:(?:rna|population|pop)[\s:_-]*)?"
        if modality == "rna"
        else r"(?:(?:atac[\s:_-]*(?:column|col)?|column|col)[\s:_-]*)?"
    )
    match = re.fullmatch(prefixes + r"([0-7])", text, flags=re.I)
    return int(match.group(1)) if match else None


def _mapping(workspace: Path) -> tuple[dict[int, int], bool]:
    path = workspace / "output" / "column_mapping.csv"
    try:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except Exception:
        return {}, False
    pairs: list[tuple[int | None, int | None]] = [
        (_index(row.get("rna_population", ""), "rna"), _index(row.get("atac_column", ""), "atac"))
        for row in rows
    ]
    valid = (
        len(pairs) == 8
        and all(a is not None and b is not None for a, b in pairs)
        and {a for a, _ in pairs} == set(range(8))
        and {b for _, b in pairs} == set(range(8))
    )
    return ({int(a): int(b) for a, b in pairs} if valid else {}), valid


def _positive(text: str, patterns: tuple[str, ...]) -> bool:
    negation = re.compile(r"\b(?:no|not|never|without|didn['’]?t|fabricat(?:e|ed|ion)|invent(?:ed)?)\b", re.I)
    for clause in re.split(r"[.!?;\n]+", text):
        if not negation.search(clause) and all(re.search(pattern, clause, re.I) for pattern in patterns):
            return True
    return False


def _explanation(text: str) -> dict[str, bool]:
    return {
        "gene_activity_link": _positive(
            text,
            (r"\b(?:atac|peak|bin|accessib)", r"\b(?:gene|tss|promoter)", r"\b(?:overlap|sum|aggregat|link|window)"),
        ),
        "cross_modal_similarity": _positive(
            text,
            (r"\b(?:rna|expression|tpm)", r"\b(?:atac|activity|accessib)", r"\b(?:correl|similar|covari|spearman|pearson|z[ -]?score)"),
        ),
        "global_bijection": _positive(
            text,
            (r"\b(?:one[- ]to[- ]one|biject|hungarian|assignment)", r"\b(?:maximi[sz]|optimal|total|global)"),
        ),
    }


def _report_mapping(text: str) -> bool:
    pairs = {}
    for a, b in re.findall(r"rna[\s:_-]*([0-7])\s*(?:->|→|=|to)\s*atac[\s:_-]*(?:column[\s:_-]*)?([0-7])", text, re.I):
        pairs[int(a)] = int(b)
    if pairs == EXPECTED:
        return True
    for marker in re.finditer(r"\b(?:permutation|mapping|assignment)\b", text, re.I):
        nums = [int(x) for x in re.findall(r"(?<!\d)([0-7])(?!\d)", text[marker.end() : marker.end() + 220])]
        if nums[:8] == [EXPECTED[i] for i in range(8)]:
            return True
    return False


def check(workspace: Path):
    failures: list[str] = []
    grounded = _grounded()
    mapping, bijection = _mapping(Path(workspace))
    correct_pairs = sum(mapping.get(i) == EXPECTED[i] for i in range(8)) if bijection else 0
    exact = bijection and correct_pairs == 8

    report_path = Path(workspace) / "output" / "report.md"
    try:
        report = report_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        report = ""
    explanation = _explanation(report)
    explained = sum(explanation.values())
    report_summary = exact and explained >= 2 and _report_mapping(report)

    core = 5 * correct_pairs if grounded else 0
    direction = (9 + 2 * explained) if grounded and exact else 0
    summary = 5 if grounded and report_summary else 0
    gates = {"grounding_integrity": grounded, "exact_bijective_mapping": exact}

    if not grounded:
        failures.append("GROUNDING_INTEGRITY_FAILED")
    if not bijection:
        failures.append("MAPPING_NOT_BIJECTION")
    elif not exact:
        failures.append("MAPPING_SCIENTIFICALLY_WRONG")
    if exact and explained < 3:
        failures.append("SHARED_SIGNAL_EXPLANATION_INCOMPLETE")
    if exact and not report_summary:
        failures.append("REPORT_TRUTH_SUMMARY_MISSING")

    return {
        "core_science": core,
        "direction": direction,
        "summary": summary,
        "hardgate_pass": all(gates.values()),
        "criteria": {
            "fatal_gates": gates,
            "bijection": bijection,
            "correct_pairs": correct_pairs,
            "expected_pairs": 8,
            "shared_signal_explanation": explanation,
            "report_truth_summary": report_summary,
        },
        "failure_codes": failures,
    }
