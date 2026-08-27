# Construct Audit Report

**Date:** 2026-08-27
**Inputs:** `inputs/constructs.csv` (3 constructs), `inputs/AUDIT_RULE.md` (frozen audit rule)
**Outputs:** `output/construct_audit.csv`, `output/analysis.py`

## Scope and method

Each construct was audited strictly against the frozen rule supplied in
`inputs/AUDIT_RULE.md`. Only fields represented in the fixture
(`insert_sequence`, `tag`, `claimed_frame`) were evaluated. No vector backbone,
tag DNA sequence, linker sequence, or other feature that is absent from the
input was inferred. The `promoter` column is present in the fixture but is not
referenced by any check in the frozen rule, so it is not used in scoring.

Checks applied (as defined in the frozen rule):

| Check | Criterion | Failure label |
|---|---|---|
| `start_ok` | Insert begins with `ATG` | `START` |
| `stop_ok` | Insert ends in an in-frame stop codon (`TAA`, `TAG`, or `TGA`). The terminal triplet is in frame only when the insert length is divisible by 3. | `STOP` |
| `frame_ok` | Insert length divisible by 3 **and** `claimed_frame` is `in_frame` | `FRAME` |
| `tag_ok` | For a `C_terminal_*` fusion: insert must not carry a terminal stop codon before the downstream tag. Any other tag string is unsupported and fails closed. | `TAG` |

`overall_status = pass` only when all four checks are true; otherwise `fail`,
with every failed check listed in `issues` (canonical order START, STOP, FRAME, TAG).

## Per-construct findings

All three inserts share the same coding body (`ATG GCC GAC TAC AAA GAC GAT GAC
GAC AAG ...`, a DYKDDDDK-like ORF) and differ only in their 3'-terminal bases,
`claimed_frame`, and promoter.

### c01 — `C_terminal_FLAG`, claimed `in_frame`
- Insert: `ATGGCCGACTACAAAGACGATGACGACAAGTAA`, length **33** (divisible by 3).
- `start_ok`: begins with `ATG` → **true**
- `stop_ok`: terminal in-frame codon is `TAA` (stop) → **true**
- `frame_ok`: 33 % 3 = 0 and `claimed_frame = in_frame` → **true**
- `tag_ok`: `C_terminal_*` fusion; the insert ends in a terminal stop codon
  (`TAA`) before the downstream FLAG tag, which would truncate translation
  before the tag → **false**
- **Overall: fail — issues: TAG**

### c02 — `C_terminal_FLAG`, claimed `in_frame`
- Insert: `ATGGCCGACTACAAAGACGATGACGACAAGAA`, length **32** (not divisible by 3).
- `start_ok`: begins with `ATG` → **true**
- `stop_ok`: length is not a multiple of 3, so there is no in-frame terminal
  codon (last three bases are `GAA`, not a stop) → **false**
- `frame_ok`: 32 % 3 = 2 → **false** (despite `claimed_frame = in_frame`)
- `tag_ok`: `C_terminal_*` fusion; the insert does not end in a stop codon, so
  no terminal stop stands before the tag → **true**
- **Overall: fail — issues: STOP; FRAME**

### c03 — `C_terminal_FLAG`, claimed `out_of_frame`
- Insert: `ATGGCCGACTACAAAGACGATGACGACAAGTGA`, length **33** (divisible by 3).
- `start_ok`: begins with `ATG` → **true**
- `stop_ok`: terminal in-frame codon is `TGA` (stop) → **true**
- `frame_ok`: length divisible by 3 but `claimed_frame = out_of_frame` → **false**
- `tag_ok`: `C_terminal_*` fusion; the insert ends in a terminal stop codon
  (`TGA`) before the downstream FLAG tag → **false**
- **Overall: fail — issues: FRAME; TAG**

## Results summary

| construct_id | frame_ok | start_ok | stop_ok | tag_ok | overall_status | issues |
|---|---|---|---|---|---|---|
| c01 | true | true | true | false | fail | TAG |
| c02 | false | true | false | true | fail | STOP;FRAME |
| c03 | false | true | true | false | fail | FRAME;TAG |

**0 of 3 constructs pass.** Each construct starts correctly with `ATG`, but each
fails at least one other check:

- **c01** fails only tag compatibility: a terminal stop codon (`TAA`) sits before
  the C-terminal FLAG tag.
- **c02** fails frame and stop: the 32 bp insert is not a multiple of three, so
  it is out of frame and lacks an in-frame terminal stop codon.
- **c03** fails frame and tag: the claimed frame is `out_of_frame`, and a
  terminal stop codon (`TGA`) sits before the C-terminal FLAG tag.

Note on rule interaction: under the frozen rule, `stop_ok` (insert must end in a
stop codon) and `tag_ok` for `C_terminal_*` fusions (insert must not end in a
stop codon before the tag) are evaluated independently. Consequently a
`C_terminal_*` construct cannot simultaneously satisfy both; this tension is
inherent to the supplied rule and is reported as-is without inferring any
unprovided vector or tag sequence.

## Reproducibility

`output/analysis.py` reads `inputs/constructs.csv`, applies exactly the checks
above, and writes `output/construct_audit.csv`. Run from the workspace root:

```
python output/analysis.py
```
