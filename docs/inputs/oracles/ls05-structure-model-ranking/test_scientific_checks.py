import importlib.util
import unittest
from pathlib import Path

HERE = Path(__file__).parent
SPEC = importlib.util.spec_from_file_location("checker_model_ranking", HERE / "scientific_checks.py")
CHECKER = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(CHECKER)


class AcceptanceTests(unittest.TestCase):
    def workspace(self, mode="reference"):
        return HERE / "testdata" / mode

    def test_reference(self):
        result = CHECKER.check(self.workspace()); self.assertEqual(result["core_science"], 40); self.assertEqual(result["direction"], 15); self.assertTrue(result["hardgate_pass"])
    def test_empty(self):
        result = CHECKER.check(self.workspace("empty")); self.assertFalse(result["hardgate_pass"]); self.assertLess(result["core_science"], 40)
    def test_wrong(self):
        result = CHECKER.check(self.workspace("wrong")); self.assertIn("ORDER_MISMATCH", result["failure_codes"]); self.assertLess(result["core_science"], 40)


if __name__ == "__main__": unittest.main()
