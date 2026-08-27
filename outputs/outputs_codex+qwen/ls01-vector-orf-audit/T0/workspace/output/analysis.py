#!/usr/bin/env python3
"""Construct audit for start/stop codons, reading frame, and tag compatibility.

Implements the frozen rule in inputs/AUDIT_RULE.md:

  * start_ok : insert begins with ATG
  * stop_ok  : insert ends in an in-frame stop codon (TAA / TAG / TGA)
  * frame_ok : insert length divisible by 3 AND claimed_frame == "in_frame"
  * tag_ok   : C_terminal_* fusions must not carry an in-frame stop codon
               before the downstream tag; any other tag string is
               unsupported and fails closed
  * overall_status = "pass" only when all four checks are true, else "fail"
    with every failed check listed in issues (labels START/STOP/FRAME/TAG)

Only the fields present in inputs/constructs.csv are examined; no vector
backbone or other sequence features are inferred.
"""

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INPUT_CSV = ROOT / "inputs" / "constructs.csv"
OUTPUT_CSV = ROOT / "output" / "construct_audit.csv"

START_CODON = "ATG"
STOP_CODONS = ("TAA", "TAG", "TGA")
CHECK_ORDER = ("START", "STOP", "FRAME", "TAG")
ISSUE_SEP = ";"
FIELDNAMES = [
    "construct_id",
    "frame_ok",
    "start_ok",
    "stop_ok",
    "tag_ok",
    "overall_status",
    "issues",
]


def normalize_seq(seq: str) -> str:
    """Strip whitespace and upper-case the sequence (input is already clean)."""
    return "".join(seq.split()).upper()


def check_start(insert: str) -> bool:
    """The insert begins with ATG."""
    return insert.startswith(START_CODON)


def check_stop(insert: str) -> bool:
    """The insert ends in an in-frame stop codon.

    The terminal triplet is only in-frame when the insert length is
    divisible by three; otherwise the sequence does not end on a codon
    boundary and cannot end in an in-frame stop codon.
    """
    return len(insert) % 3 == 0 and insert[-3:] in STOP_CODONS


def check_frame(insert: str, claimed_frame: str) -> bool:
    """Insert length divisible by three AND claimed_frame is 'in_frame'."""
    return len(insert) % 3 == 0 and claimed_frame.strip() == "in_frame"


def in_frame_stop_codons(insert: str):
    """Return (codon index, triplet) for each in-frame stop among complete codons."""
    return [
        (i // 3, insert[i:i + 3])
        for i in range(0, len(insert) - 2, 3)
        if insert[i:i + 3] in STOP_CODONS
    ]


def check_tag(insert: str, tag: str) -> bool:
    """Tag/linker compatibility for the tag strings represented in the input.

    For a C_terminal_* fusion, translation must read through the insert
    into the downstream tag, so no in-frame (terminal) stop codon may
    occur inside the insert. Any other tag string is unsupported and
    fails closed.
    """
    tag = tag.strip()
    if tag.startswith("C_terminal_"):
        return not in_frame_stop_codons(insert)
    return False


def audit_row(row):
    """Audit one construct; return (csv_record, diagnostics_dict)."""
    insert = normalize_seq(row["insert_sequence"])
    claimed_frame = row["claimed_frame"].strip()
    tag = row["tag"].strip()

    passed = {
        "START": check_start(insert),
        "STOP": check_stop(insert),
        "FRAME": check_frame(insert, claimed_frame),
        "TAG": check_tag(insert, tag),
    }
    issues = [name for name in CHECK_ORDER if not passed[name]]
    record = {
        "construct_id": row["construct_id"].strip(),
        "frame_ok": str(passed["FRAME"]).lower(),
        "start_ok": str(passed["START"]).lower(),
        "stop_ok": str(passed["STOP"]).lower(),
        "tag_ok": str(passed["TAG"]).lower(),
        "overall_status": "pass" if not issues else "fail",
        "issues": ISSUE_SEP.join(issues),
    }
    diag = {
        "insert_len": len(insert),
        "claimed_frame": claimed_frame,
        "tag": tag,
        "stops": in_frame_stop_codons(insert),
    }
    return record, diag


def main():
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with INPUT_CSV.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    records = []
    for row in rows:
        record, diag = audit_row(row)
        records.append(record)
        stops = ", ".join(f"{triplet}@codon{idx + 1}" for idx, triplet in diag["stops"]) or "-"
        print(
            f"{record['construct_id']}: len={diag['insert_len']} "
            f"claimed_frame={diag['claimed_frame']} tag={diag['tag']} stops=[{stops}] -> "
            f"{record['overall_status']} "
            f"(start_ok={record['start_ok']}, stop_ok={record['stop_ok']}, "
            f"frame_ok={record['frame_ok']}, tag_ok={record['tag_ok']}) "
            f"issues={record['issues'] or '-'}"
        )

    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(records)

    print(f"\nWrote {OUTPUT_CSV} ({len(records)} constructs)")


if __name__ == "__main__":
    main()