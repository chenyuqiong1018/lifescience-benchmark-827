# Oracle contract

Each task has an `oracle.py` entry point. The shared runner performs only static artifact parsing and static Python AST checks; it never imports or executes the submission's script.

To activate scoring, add a task-local `scientific_checks.py` with `ACCEPTED = True` and `check(workspace)`, but only after reference, empty-output and scientifically-wrong fixtures each pass/fail as intended in 3/3 isolated runs. `check` returns `core_science` (0–40), `direction` (0–15), `summary` (0–5), `hardgate_pass`, `failure_codes`, and criterion details. Until then the entry point fails closed with exit code 2 and `grader_status=blocked`.

Usage: `python docs/oracles/<task-id>/oracle.py --workspace <frozen-workspace> --json-out result.json`.
