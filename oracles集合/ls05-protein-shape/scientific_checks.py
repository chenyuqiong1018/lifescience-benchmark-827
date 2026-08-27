from __future__ import annotations

import json
import math
import re
from pathlib import Path


ACCEPTED = True

# Independently derived from protein.shape.q1.pdb (SHA-256 recorded in the
# author result): an XY projection, viewed along Z and rotated a quarter turn,
# has a long vertical spine with terminal and medial arms on the same side.
EXPECTED_LETTER = "F"
ALLOWED_LETTERS = frozenset("BDFHJL NPRTVXZ".replace(" ", ""))

_NEGATION = re.compile(
    r"\b(?:not|no|never|without|neither|isn['’]?t|doesn['’]?t|cannot|can['’]?t|fails?\s+to)\b",
    re.IGNORECASE,
)


def _affirmed(text: str, patterns: tuple[str, ...]) -> bool:
    """Find a positive assertion, excluding locally negated evidence."""
    for pattern in patterns:
        for match in re.finditer(pattern, text, re.IGNORECASE):
            prefix = text[max(0, match.start() - 48) : match.start()]
            # A sentence boundary makes an earlier negation irrelevant.
            prefix = re.split(r"[.;!?\n]", prefix)[-1]
            if not _NEGATION.search(prefix):
                return True
    return False


def _orientation_evidence(notes: str) -> tuple[bool, bool, bool]:
    projection = _affirmed(
        notes,
        (
            r"\bxy[-\s]?plane\b",
            r"\bx[-\s]?y\s+projection\b",
            r"\bproject(?:ed|ion)?\b.{0,28}\b(?:along|down|onto)\b.{0,12}\b[+\-±]?z(?:[-\s]?axis)?\b",
            r"\bview(?:ing)?\s+(?:vector|direction)\b.{0,20}\b0\s*[, ]\s*0\s*[, ]\s*[+\-]?1\b",
            r"\b[+\-±]?z(?:[-\s]?axis)?\b.{0,24}\b(?:line of sight|view direction|orthographic)\b",
        ),
    )

    explicit_turn = _affirmed(
        notes,
        (
            r"\b(?:rotat(?:e|ed|ion)|turn(?:ed)?)\b.{0,20}\b(?:90\s*(?:°|degrees?)|quarter[-\s]?turn)\b",
            r"\b(?:90\s*(?:°|degrees?)|quarter[-\s]?turn)\b.{0,20}\b(?:rotat(?:e|ed|ion)|turn(?:ed)?)\b",
        ),
    )
    axis_mapping = (
        _affirmed(notes, (r"\b(?:pdb\s+)?x(?:[-\s]?axis|\s+extent)?\b.{0,28}\bvertical\b",))
        and _affirmed(notes, (r"\bnegative[-\s]?y|\blow[-\s]?y\b",))
        and _affirmed(notes, (r"\bright(?:ward)?\b|\blateral\b",))
    )
    camera_mapping = (
        projection
        and _affirmed(notes, (r"\bup\s+vector\b.{0,20}\b[+]?1\s*[, ]\s*0\s*[, ]\s*0\b",))
    )
    quarter_turn = explicit_turn or axis_mapping or camera_mapping

    spine = _affirmed(
        notes,
        (r"\b(?:long|main|vertical)\b.{0,22}\b(?:vertical\s+)?(?:spine|stem|backbone|shaft|axis)\b",),
    )
    two_named_arms = (
        _affirmed(notes, (r"\b(?:top|upper|terminal)\b.{0,20}\b(?:arm|branch|lobe|protrusion)\b",))
        and _affirmed(notes, (r"\b(?:mid(?:dle)?|medial|central)\b.{0,20}\b(?:arm|branch|lobe|protrusion)\b",))
    )
    two_counted_arms = _affirmed(
        notes,
        (r"\b(?:two|2)\b.{0,28}\b(?:arms|branches|lobes|protrusions)\b",),
    )
    same_side = _affirmed(
        notes,
        (r"\bsame[-\s]?side\b", r"\bright(?:ward)?\b.{0,20}\b(?:arms?|branches?|lobes?|protrusions?)\b",
         r"\b(?:arms?|branches?|lobes?|protrusions?)\b.{0,20}\bright(?:ward)?\b"),
    )
    topology = (spine and (two_named_arms or two_counted_arms)) or ((two_named_arms or two_counted_arms) and same_side)
    return projection, quarter_turn, topology


def check(workspace: Path):
    criteria: dict[str, bool] = {}
    failures: list[str] = []
    try:
        data = json.loads((Path(workspace) / "output" / "shape_call.json").read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            data = {}
    except Exception:
        data = {}

    raw_letter = data.get("letter")
    letter = str(raw_letter).strip().upper() if isinstance(raw_letter, str) else ""
    notes = data.get("orientation_notes")
    notes = notes.strip() if isinstance(notes, str) else ""

    projection, quarter_turn, topology = _orientation_evidence(notes)
    evidence_count = sum((projection, quarter_turn, topology))

    identity = letter == EXPECTED_LETTER
    allowed_call = len(letter) == 1 and letter in ALLOWED_LETTERS
    grounded_orientation = evidence_count >= 2

    # Exactly two explicitly named fatal scientific gates.
    gate_shape_identity = identity
    gate_grounded_orientation = grounded_orientation

    try:
        confidence = float(data.get("confidence"))
        confidence_valid = math.isfinite(confidence) and 0.0 <= confidence <= 1.0
    except (TypeError, ValueError):
        confidence_valid = False

    criteria.update(
        allowed_shape_call=allowed_call,
        expected_shape_F=identity,
        xy_projection_grounded=projection,
        quarter_turn_grounded=quarter_turn,
        f_topology_grounded=topology,
        confidence_bounded=confidence_valid,
        orientation_notes_present=bool(notes),
        GATE_SHAPE_IDENTITY=gate_shape_identity,
        GATE_GROUNDED_ORIENTATION=gate_grounded_orientation,
    )

    core_science = (32 if identity else 0) + (8 if identity and topology else 0)
    direction = (6 if projection else 0) + (5 if quarter_turn else 0) + (4 if topology else 0)
    summary = 5 if identity and grounded_orientation and confidence_valid and bool(notes) else 0
    hardgate_pass = gate_shape_identity and gate_grounded_orientation

    if not allowed_call:
        failures.append("SHAPE_CALL_INVALID")
    if not gate_shape_identity:
        failures.append("FATAL_SHAPE_IDENTITY_MISMATCH")
    if not gate_grounded_orientation:
        failures.append("FATAL_ORIENTATION_UNGROUNDED")
    if not confidence_valid:
        failures.append("CONFIDENCE_INVALID")
    if direction < 15:
        failures.append("ORIENTATION_INCOMPLETE")

    return {
        "core_science": core_science,
        "direction": direction,
        "summary": summary,
        "hardgate_pass": hardgate_pass,
        "criteria": criteria,
        "failure_codes": failures,
    }
