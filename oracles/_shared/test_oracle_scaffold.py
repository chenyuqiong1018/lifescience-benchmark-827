import json
import subprocess
import sys
import unittest
import uuid
from pathlib import Path

from oracle_common import REQUIRED_OUTPUTS


class OracleScaffoldTests(unittest.TestCase):
    def test_all_25_tasks_have_entry_points(self):
        root = Path(__file__).resolve().parent.parent
        self.assertEqual(25, len(REQUIRED_OUTPUTS))
        for task_id in REQUIRED_OUTPUTS:
            self.assertTrue((root / task_id / "oracle.py").is_file(), task_id)

    def test_unaccepted_oracle_fails_closed(self):
        root = Path(__file__).resolve().parent.parent
        task_id = "ls01-grna-offtarget-rank"
        repo_root = root.parent.parent
        workspace = repo_root / ".tmp_tests" / f"oracle-{uuid.uuid4().hex}"
        (workspace / "output").mkdir(parents=True)
        completed = subprocess.run(
            [sys.executable, str(root / task_id / "oracle.py"), "--workspace", str(workspace)],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(2, completed.returncode)
        result = json.loads(completed.stdout)
        self.assertEqual("blocked", result["grader_status"])
        self.assertIsNone(result["deterministic_score"])
        self.assertIn("ORACLE_NOT_ACCEPTED", result["failure_codes"])


if __name__ == "__main__":
    unittest.main()
