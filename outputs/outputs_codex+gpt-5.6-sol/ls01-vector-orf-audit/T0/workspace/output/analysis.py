from __future__ import annotations

import csv
from pathlib import Path


OUT = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[6]
SOURCE = ROOT / "inputs" / "ls01-vector-orf-audit" / "constructs.csv"
STOPS = {"TAA", "TAG", "TGA"}
DNA = set("ACGT")


def main() -> None:
    with SOURCE.open(newline="", encoding="utf-8") as handle:
        records = list(csv.DictReader(handle))
    if len({row["construct_id"] for row in records}) != len(records):
        raise ValueError("construct_id values must be unique")

    output: list[dict[str, str]] = []
    for row in records:
        sequence = row["insert_sequence"].upper()
        if not sequence or set(sequence) - DNA:
            raise ValueError(f"invalid DNA insert for {row['construct_id']}")
        start_ok = sequence[:3] == "ATG"
        stop_ok = sequence[-3:] in STOPS
        frame_ok = len(sequence) % 3 == 0 and row["claimed_frame"] == "in_frame"
        tag_supported = row["tag"].startswith("C_terminal_")
        tag_ok = tag_supported and not stop_ok
        checks = [("START", start_ok), ("STOP", stop_ok), ("FRAME", frame_ok), ("TAG", tag_ok)]
        issues = [name for name, passed in checks if not passed]
        output.append(
            {
                "construct_id": row["construct_id"],
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
        writer.writerows(output)


if __name__ == "__main__":
    main()
