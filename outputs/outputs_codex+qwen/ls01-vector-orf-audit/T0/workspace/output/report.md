# Construct Audit Report

**Date:** 2026-08-27
**Inputs:** `inputs/constructs.csv` (3 constructs); frozen rule in `inputs/AUDIT_RULE.md`
**Deliverables:** `output/construct_audit.csv`, `output/analysis.py`, this report

## Scope

Each construct was audited strictly on the fields represented in the input
fixture (`construct_id`, `promoter`, `insert_sequence`, `tag`,
`claimed_frame`). No vector backbone, linker sequence, or any other sequence
feature was inferred. The `promoter` field is carried in the input but does
not participate in any audit check.

## Checks applied (frozen rule)

| Check | Criterion |
|---|---|
| `start_ok` | Insert begins with `ATG`. |
| `stop_ok` | Insert ends in an in-frame stop codon (`TAA`, `TAG`, or `TGA`). The terminal triplet is only in-frame when the insert length is divisible by 3, i.e. the sequence ends on a codon boundary. |
| `frame_ok` | Insert length divisible by 3 **and** `claimed_frame` is `in_frame`. |
| `tag_ok` | For a `C_terminal_*` fusion, the insert must not contain a terminal (in-frame) stop codon before the downstream tag, since it would terminate translation before the tag. Any other tag string is unsupported and fails closed. |

`overall_status` is `pass` only when all four checks are true; otherwise it
is `fail`, and `issues` lists every failed check using the labels `START`,
`STOP`, `FRAME`, `TAG` (joined with `;`, in that canonical order).

## Results summary

| construct_id | frame_ok | start_ok | stop_ok | tag_ok | overall_status | issues |
|---|---|---|---|---|---|---|
| c01 | true | true | true | false | fail | TAG |
| c02 | false | true | false | true | fail | STOP;FRAME |
| c03 | false | true | true | false | fail | FRAME;TAG |

**0 of 3 constructs pass.**

## Per-construct findings

### c01 — EF1a / C_terminal_FLAG / claimed in_frame
- Insert length 33 nt (11 codons), divisible by 3.
- Begins with `ATG` → start_ok.
- Ends in in-frame `TAA` (codon 11) → stop_ok.
- Length % 3 = 0 and claimed frame `in_frame` → frame_ok.
- C-terminal FLAG fusion: the in-frame `TAA` at codon 11 terminates
  translation before the downstream FLAG tag → tag_ok = false.
- **Status: fail — issues: TAG.**

### c02 — CMV / C_terminal_FLAG / claimed in_frame
- Insert length 32 nt — not divisible by 3 (10 complete codons + 2 trailing nt).
- Begins with `ATG` → start_ok.
- Does not end in an in-frame stop codon: the sequence does not end on a
  codon boundary and its terminal triplet `GAA` is not a stop → stop_ok = false.
- Length not divisible by 3 → frame_ok = false despite the `in_frame` claim.
- C-terminal FLAG fusion: no in-frame stop codon occurs among the insert's
  complete codons, so nothing represented in the input blocks read-through
  into the tag → tag_ok = true.
- **Status: fail — issues: STOP;FRAME.**

### c03 — EF1a / C_terminal_FLAG / claimed out_of_frame
- Insert length 33 nt (11 codons), divisible by 3.
- Begins with `ATG` → start_ok.
- Ends in in-frame `TGA` (codon 11) → stop_ok.
- `claimed_frame` is `out_of_frame`; the rule requires both divisibility by
  3 and an `in_frame` claim → frame_ok = false.
- C-terminal FLAG fusion: the in-frame `TGA` at codon 11 terminates
  translation before the downstream tag → tag_ok = false.
- **Status: fail — issues: FRAME;TAG.**

## Interpretation notes

- "Terminal stop codon before the downstream tag" is evaluated on the
  complete in-frame codons of the insert as supplied. Any in-frame stop
  codon inside the insert (including its final codon) terminates translation
  before a C-terminal tag; in this fixture the offending stops are the final
  codons of c01 (`TAA`) and c03 (`TGA`).
- `stop_ok` requires the terminal triplet to be in-frame. c02 (length 32)
  does not end on a codon boundary, so it cannot end in an in-frame stop.
- Only `C_terminal_*` tag strings are supported by the frozen rule; any other
  tag value fails closed. All constructs in this fixture carry
  `C_terminal_FLAG`.
- No features beyond the fixture columns were inferred or assumed (no linker
  sequence, vector backbone, or expression-context checks).

## Reproduction

```
python output/analysis.py
```

The script reads `inputs/constructs.csv`, applies the frozen rule, prints
per-construct diagnostics, and writes `output/construct_audit.csv`.