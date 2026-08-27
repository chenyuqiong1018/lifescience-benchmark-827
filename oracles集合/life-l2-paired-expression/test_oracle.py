#!/usr/bin/env python3
"""Self-contained tests for the paired-expression oracle."""

from __future__ import annotations

import csv
import importlib.util
import json
import shutil
import struct
import tempfile
import unittest
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUT_FIXTURE = HERE.parents[1] / "inputs" / "life-l2-paired-expression" / "paired_expression.csv"
SPEC = importlib.util.spec_from_file_location("paired_oracle", HERE / "oracle.py")
ORACLE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(ORACLE)


def chunk(kind: bytes, data: bytes) -> bytes:
    return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)


def tiny_png() -> bytes:
    raw = b"\x00\xff\xff\xff"
    return b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", 1, 1, 8, 2, 0, 0, 0)) + chunk(b"IDAT", zlib.compress(raw)) + chunk(b"IEND", b"")


class OracleTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.workspace = Path(self.temp.name)
        (self.workspace / "inputs").mkdir()
        (self.workspace / "output").mkdir()
        shutil.copy2(INPUT_FIXTURE, self.workspace / "inputs" / "paired_expression.csv")

    def tearDown(self) -> None:
        self.temp.cleanup()

    def write_correct(self) -> None:
        with INPUT_FIXTURE.open(newline="", encoding="utf-8") as handle:
            source = list(csv.DictReader(handle))
        changes = []
        with (self.workspace / "output" / "paired_expression_results.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=ORACLE.EXPECTED_COLUMNS)
            writer.writeheader()
            for row in source:
                change = float(row["post_expression"]) - float(row["pre_expression"])
                changes.append(change)
                writer.writerow({**row, "paired_change": change, "direction": ORACLE.direction(change)})
        directions = [ORACLE.direction(value) for value in changes]
        summary = {
            "task_id": ORACLE.TASK_ID,
            "pair_count": len(changes),
            "mean_paired_change": sum(changes) / len(changes),
            "up_count": directions.count("up"),
            "down_count": directions.count("down"),
            "unchanged_count": directions.count("unchanged"),
        }
        (self.workspace / "output" / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
        (self.workspace / "output" / "paired_expression.png").write_bytes(tiny_png())

    def test_correct_output_scores_80(self) -> None:
        self.write_correct()
        result = ORACLE.grade(self.workspace)
        self.assertEqual(result["deterministic_score"], 80)
        self.assertTrue(result["hardgate_pass"])

    def test_empty_output_fails(self) -> None:
        result = ORACLE.grade(self.workspace)
        self.assertEqual(result["deterministic_score"], 0)
        self.assertFalse(result["hardgate_pass"])

    def test_scientifically_wrong_change_fails(self) -> None:
        self.write_correct()
        result_path = self.workspace / "output" / "paired_expression_results.csv"
        text = result_path.read_text(encoding="utf-8").replace("S01,8.0,10.0,2.0,up", "S01,8.0,10.0,-2.0,down")
        result_path.write_text(text, encoding="utf-8")
        result = ORACLE.grade(self.workspace)
        self.assertLess(result["deterministic_score"], 80)
        self.assertIn("PAIRED_CHANGE_INCORRECT", result["failure_codes"])


if __name__ == "__main__":
    unittest.main()

