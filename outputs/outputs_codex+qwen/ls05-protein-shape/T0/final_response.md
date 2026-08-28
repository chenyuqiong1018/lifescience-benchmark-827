Both deliverables are complete and validated.

**Result: the structure is a "Z".**

- **Evidence (primary method):** CA-backbone trace projected onto the XY plane (the structure is near-planar, best-fit plane z ≈ −0.19·y with ~5 Å residual). The smoothed trace shows exactly three strokes: a top horizontal bar (x −175→75 at y≈+55), a diagonal from the top bar's right end to the bottom bar's left end, and a bottom horizontal bar (x −25→120 at y≈−30) — the canonical Z topology (diagonal joining opposite bar ends, ruling out T/H/F junctions; diagonal direction ruling out S).
- **Validation:** multi-view inspection (XY/XZ/YZ + 3D) confirmed XY is the letter plane; topology/junction checks and letter-template comparisons confirmed Z over N (bars horizontal, not vertical) and over T/S.

**Artifacts:**
- `output/shape_call.json` — `{"letter": "Z", "confidence": 0.92, "orientation_notes": ...}` (schema asserted: exact keys, allowed letter, 0–1 confidence)
- `output/shape_view.png` — backbone projection showing the Z.