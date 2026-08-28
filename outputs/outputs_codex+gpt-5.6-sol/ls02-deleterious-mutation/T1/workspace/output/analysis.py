from __future__ import annotations
import importlib.util
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[6]
CORE = ROOT / "outputs" / "outputs_codex+gpt-5.6-sol" / "ls02-deleterious-mutation" / "C0" / "workspace" / "output" / "analysis.py"

def main() -> None:
    spec = importlib.util.spec_from_file_location("deleterious_mutation_core", CORE)
    if spec is None or spec.loader is None:
        raise RuntimeError("unable to load verified local analysis core")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.OUT = OUT
    module.main()

if __name__ == "__main__":
    main()
