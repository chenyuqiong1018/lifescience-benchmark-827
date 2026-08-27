# Frozen construct-audit rule

- `start_ok`: the insert begins with `ATG`.
- `stop_ok`: the insert ends in an in-frame stop codon (`TAA`, `TAG`, or `TGA`).
- `frame_ok`: insert length is divisible by three **and** `claimed_frame` is `in_frame`.
- `tag_ok`: for a `C_terminal_*` fusion, the insert must not contain a terminal stop codon before the downstream tag; other tag strings are unsupported and fail closed.
- `overall_status=pass` only when all four checks are true; otherwise use `fail` and list every failed check in `issues` using the labels `START`, `STOP`, `FRAME`, and `TAG`.

This rule checks only fields represented by the fixture and does not infer vector sequence that was not supplied.
