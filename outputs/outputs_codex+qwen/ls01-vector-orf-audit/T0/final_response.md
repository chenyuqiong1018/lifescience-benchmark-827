Audit complete. All three deliverables are in `output/`, built strictly from the frozen rule in `inputs/AUDIT_RULE.md` and the fields present in `inputs/constructs.csv` (no vector/linker features inferred).

**Results (`output/construct_audit.csv`)** — 0 of 3 constructs pass:

| construct_id | frame_ok | start_ok | stop_ok | tag_ok | overall_status | issues |
|---|---|---|---|---|---|---|
| c01 | true | true | true | false | fail | TAG |
| c02 | false | true | false | true | fail | STOP;FRAME |
| c03 | false | true | true | false | fail | FRAME;TAG |

**Key findings:**
- **c01** (33 nt, in-frame): starts with ATG, ends in-frame `TAA`, frame consistent — but as a `C_terminal_FLAG` fusion, the terminal `TAA` (codon 11) would stop translation before the FLAG tag → TAG failure.
- **c02** (32 nt): starts with ATG and has no in-frame stop inside the insert (so tag_ok is true), but 32 nt is not divisible by 3 — it ends mid-codon (`...AAG|AA`), so it neither ends in an in-frame stop nor satisfies the frame check despite the `in_frame` claim → STOP + FRAME failures.
- **c03** (33 nt): starts with ATG and ends in-frame `TGA`, but `claimed_frame=out_of_frame` fails the frame rule, and the `TGA` before the C-terminal FLAG tag fails the tag rule → FRAME + TAG failures.

**Deliverables:**
- `output/construct_audit.csv` — machine-readable audit with the exact requested columns
- `output/analysis.py` — reproducible auditor (`python output/analysis.py`); stdlib-only, prints per-construct diagnostics
- `output/report.md` — rule summary, per-construct codon-level findings, and interpretation notes (e.g., stop_ok requires the terminal triplet to be in-frame; only `C_terminal_*` tags are supported, others fail closed)