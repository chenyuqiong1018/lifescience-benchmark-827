from __future__ import annotations

import csv
from pathlib import Path


OUT = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[6]
INPUT = ROOT / "inputs" / "ls01-vector-orf-audit" / "constructs.csv"
STOPS = {"TAA", "TAG", "TGA"}


def main() -> None:
    with INPUT.open(newline="", encoding="utf-8") as handle:
        constructs = list(csv.DictReader(handle))

    rows: list[dict[str, str]] = []
    for construct in constructs:
        sequence = construct["insert_sequence"].upper()
        start_ok = sequence.startswith("ATG")
        stop_ok = len(sequence) >= 3 and sequence[-3:] in STOPS
        frame_ok = len(sequence) % 3 == 0 and construct["claimed_frame"] == "in_frame"
        tag = construct["tag"]
        tag_ok = tag.startswith("C_terminal_") and not stop_ok

        checks = [
            ("START", start_ok),
            ("STOP", stop_ok),
            ("FRAME", frame_ok),
            ("TAG", tag_ok),
        ]
        issues = [label for label, passed in checks if not passed]
        rows.append(
            {
                "construct_id": construct["construct_id"],
                "frame_ok": str(frame_ok).lower(),
                "start_ok": str(start_ok).lower(),
                "stop_ok": str(stop_ok).lower(),
                "tag_ok": str(tag_ok).lower(),
                "overall_status": "pass" if not issues else "fail",
                "issues": ";".join(issues),
            }
        )

    fields = ["construct_id", "frame_ok", "start_ok", "stop_ok", "tag_ok", "overall_status", "issues"]
    with (OUT / "construct_audit.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
