from __future__ import annotations

import importlib.util
import json
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("ls05_scientific_checks", HERE / "scientific_checks.py")
CHECKER = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(CHECKER)


def _run(payload=None):
    with tempfile.TemporaryDirectory() as tmp:
        workspace = Path(tmp)
        if payload is not None:
            out = workspace / "output"
            out.mkdir()
            (out / "shape_call.json").write_text(json.dumps(payload), encoding="utf-8")
        # Candidate code is never imported, invoked, or inspected.
        return CHECKER.check(workspace)


def _case(name, fn):
    try:
        fn()
        return {"name": name, "passed": True}
    except AssertionError as exc:
        return {"name": name, "passed": False, "detail": str(exc)}


def reference_like_correct():
    result = _run({
        "letter": "F",
        "confidence": 0.86,
        "orientation_notes": (
            "Viewed orthographically down the Z axis (the XY projection), then rotated 90 degrees "
            "in-plane. The long body becomes a vertical spine; the terminal upper arm and thinner "
            "middle arm project rightward on the same side."
        ),
    })
    assert result["hardgate_pass"] is True, result
    assert (result["core_science"], result["direction"], result["summary"]) == (40, 15, 5), result


def empty_or_missing():
    missing = _run()
    empty = _run({})
    for result in (missing, empty):
        assert result["hardgate_pass"] is False, result
        assert (result["core_science"], result["direction"], result["summary"]) == (0, 0, 0), result
        assert "FATAL_SHAPE_IDENTITY_MISMATCH" in result["failure_codes"], result


def scientifically_wrong():
    result = _run({
        "letter": "T",
        "confidence": 0.91,
        "orientation_notes": (
            "XY projection viewed down Z and rotated by 90 degrees. A long vertical spine and "
            "two rightward arms, one upper and one middle, are visible."
        ),
    })
    assert result["hardgate_pass"] is False, result
    assert result["core_science"] == 0, result
    assert "FATAL_SHAPE_IDENTITY_MISMATCH" in result["failure_codes"], result


def plausible_but_unsupported():
    vague = _run({"letter": "F", "confidence": 0.99, "orientation_notes": "It looks convincingly like F."})
    assert vague["hardgate_pass"] is False, vague
    assert vague["direction"] == 0 and vague["summary"] == 0, vague

    # Negated keywords must not be mistaken for positive scientific evidence.
    negated = _run({
        "letter": "F",
        "confidence": 0.9,
        "orientation_notes": "This is not an XY projection, was not rotated 90 degrees, and has no vertical spine or two arms.",
    })
    assert negated["hardgate_pass"] is False, negated
    assert negated["direction"] == 0, negated

    # A fabricated/out-of-range number earns no summary credit even when geometry is stated correctly.
    fabricated = _run({
        "letter": "F",
        "confidence": 7.3,
        "orientation_notes": (
            "XY-plane projection down Z, rotated 90 degrees: the long vertical spine has upper "
            "and middle rightward arms."
        ),
    })
    assert fabricated["hardgate_pass"] is True, fabricated
    assert fabricated["summary"] == 0 and not fabricated["criteria"]["confidence_bounded"], fabricated


def valid_alternative_implementation():
    result = _run({
        "letter": "f",
        "confidence": 0.72,
        "orientation_notes": (
            "Use camera view vector [0, 0, 1] with up vector [1, 0, 0]. In that orthographic "
            "layout, the two low-Y lobes become rightward lateral protrusions at the terminal "
            "and medial positions of the long main axis."
        ),
        "method": "alternative renderer and JSON layout are intentionally allowed",
    })
    assert result["hardgate_pass"] is True, result
    assert (result["core_science"], result["direction"], result["summary"]) == (40, 15, 5), result


CASES = (
    ("reference_like_correct", reference_like_correct),
    ("empty_or_missing", empty_or_missing),
    ("scientifically_wrong", scientifically_wrong),
    ("plausible_but_unsupported", plausible_but_unsupported),
    ("valid_alternative_implementation", valid_alternative_implementation),
)


if __name__ == "__main__":
    results = [_case(name, fn) for name, fn in CASES]
    payload = {
        "schema_version": 1,
        "all_passed": all(item["passed"] for item in results),
        "candidate_code_executed": False,
        "cases": results,
    }
    (HERE / "acceptance-result.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    raise SystemExit(0 if payload["all_passed"] else 1)
