All three deliverables are complete and verified. Here's a summary of the audit:

## Deliverables (in `output/`)
- **`construct_audit.csv`** — audit results in the requested column order
- **`analysis.py`** — self-contained script that reads `inputs/constructs.csv`, applies the frozen rules from `inputs/AUDIT_RULE.md`, and regenerates both the CSV and report
- **`report.md`** — full report with rules applied, results table, per-construct explanations, and scope notes

## Audit results

| construct | frame_ok | start_ok | stop_ok | tag_ok | status | issues |
|---|---|---|---|---|---|---|
| c01 | True | True | True | False | fail | TAG |
| c02 | False | True | False | True | fail | STOP;FRAME |
| c03 | False | True | True | False | fail | FRAME;TAG |

**Reasoning per the frozen rules:**
- **c01** (33 nt, in-frame, ends `TAA`): passes start/stop/frame, but as a `C_terminal_FLAG` fusion the terminal `TAA` stop would terminate translation before the downstream tag → **TAG** failure.
- **c02** (32 nt): length not divisible by 3, so no in-frame terminal stop codon (**STOP**) and frame check fails (**FRAME**); its tail (`...AAGAA`) contains no terminal stop, so tag_ok is True.
- **c03** (33 nt, ends `TGA`): `claimed_frame=out_of_frame` fails the frame check (**FRAME**), and the terminal `TGA` blocks the C-terminal FLAG tag (**TAG**).

Per the instruction not to infer absent features: the fixture supplies no linker or cloning-flag fields, so none were evaluated; `promoter` is recorded but has no audit rule. The report also notes that under these frozen rules, `stop_ok` and `tag_ok` are mutually exclusive for `C_terminal_*` fusions, so no construct in this fixture can pass — all checks were nonetheless reported independently as specified.