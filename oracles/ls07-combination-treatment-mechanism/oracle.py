from pathlib import Path
import sys

_SHARED = Path(__file__).resolve().parent.parent / "_shared"
sys.path.insert(0, str(_SHARED))
from oracle_common import run

if __name__ == "__main__":
    raise SystemExit(run(Path(__file__).resolve().parent.name))

