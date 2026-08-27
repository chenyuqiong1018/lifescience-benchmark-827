"""Construct audit analysis.

Audits each construct in inputs/constructs.csv against the frozen rules in
inputs/AUDIT_RULE.md:

  - start_ok : the insert begins with ATG.
  - stop_ok  : the insert ends in an in-frame stop codon (TAA, TAG, or TGA).
  - frame_ok : insert length is divisible by three AND claimed_frame is "in_frame".
  - tag_ok   : for a C_terminal_* fusion, the insert must not contain a terminal
               stop codon before the downstream tag; other tag strings are
               unsupported and fail closed.
  - overall_status = pass only when all four checks are true; otherwise fail and
    list every failed check in issues using labels START, STOP, FRAME, TAG.

Only fields present in the input fixture are checked. No vector or tag sequence
that was not supplied is inferred.
"""

from __future__ import annotations

import csv
from pathlib import Path

START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}
C_TERMINAL_TAG_PREFIX = "C_terminal_"

# Canonical label order for failed checks.
CHECK_LABELS = [("start_ok", "START"), ("stop_ok", "STOP"),
                ("frame_ok", "FRAME"), ("tag_ok", "TAG")]

INPUT_CSV = Path(__file__).resolve().parent.parent / "inputs" / "constructs.csv"
OUTPUT_CSV = Path(__file__).resolve().parent / "construct_audit.csv"


def normalize(seq: str) -> str:
    """Upper-case the sequence and strip whitespace/newlines."""
    return "".join(seq.split()).upper()


def check_start(seq: str) -> bool:
    """Insert begins with ATG."""
    return seq.startswith(START_CODON)


def check_stop(seq: str) -> bool:
    """Insert ends in an in-frame stop codon.

    The terminal triplet is only in-frame when the insert length is a multiple
    of three, so that the last three bases form a complete codon aligned with
    the frame that starts at the first base.
    """
    if len(seq) % 3 != 0:
        return False
    return seq[-3:] in STOP_CODONS


def check_frame(seq: str, claimed_frame: str) -> bool:
    """Insert length divisible by three AND claimed_frame is in_frame."""
    return len(seq) % 3 == 0 and claimed_frame.strip() == "in_frame"


def check_tag(seq: str, tag: str) -> bool:
    """Tag/linker compatibility.

    For a C_terminal_* fusion the insert must not carry a terminal stop codon
    before the downstream tag (a stop there would truncate translation before
    the tag). Any tag string that is not a C_terminal_* fusion is unsupported
    and fails closed.
    """
    tag = tag.strip()
    if not tag.startswith(C_TERMINAL_TAG_PREFIX):
        return False  # unsupported tag string -> fail closed
    # Terminal stop codon check: last three bases form a stop codon.
    terminal_stop = len(seq) >= 3 and seq[-3:] in STOP_CODONS
    return not terminal_stop


def audit_construct(row: dict) -> dict:
    seq = normalize(row["insert_sequence"])
    claimed_frame = row.get("claimed_frame", "")
    tag = row.get("tag", "")

    checks = {
        "start_ok": check_start(seq),
        "stop_ok": check_stop(seq),
        "frame_ok": check_frame(seq, claimed_frame),
        "tag_ok": check_tag(seq, tag),
    }

    issues = [label for key, label in CHECK_LABELS if not checks[key]]
    overall = "pass" if not issues else "fail"

    return {
        "construct_id": row["construct_id"],
        "frame_ok": str(checks["frame_ok"]).lower(),
        "start_ok": str(checks["start_ok"]).lower(),
        "stop_ok": str(checks["stop_ok"]).lower(),
        "tag_ok": str(checks["tag_ok"]).lower(),
        "overall_status": overall,
        "issues": ";".join(issues),
    }


def main() -> None:
    with INPUT_CSV.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    results = [audit_construct(row) for row in rows]

    fieldnames = ["construct_id", "frame_ok", "start_ok", "stop_ok",
                  "tag_ok", "overall_status", "issues"]
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"Wrote {OUTPUT_CSV} with {len(results)} rows")
    for r in results:
        print(r)


if __name__ == "__main__":
    main()
