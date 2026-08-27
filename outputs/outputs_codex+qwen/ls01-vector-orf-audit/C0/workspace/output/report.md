# Construct Audit Report

Audit of the constructs in `inputs/constructs.csv` against the frozen rules in `inputs/AUDIT_RULE.md`. Only fields represented in the input fixture are checked; no vector sequence or other sequence features are inferred.

## Rules applied

| Check | Rule |
| --- | --- |
| `start_ok` | Insert begins with `ATG`. |
| `stop_ok` | Insert ends in an in-frame stop codon (`TAA`/`TAG`/`TGA`). |
| `frame_ok` | Insert length divisible by 3 **and** `claimed_frame` is `in_frame`. |
| `tag_ok` | For `C_terminal_*` fusions the insert must not contain a terminal stop codon before the downstream tag; other tag strings are unsupported and fail closed. |
| `overall_status` | `pass` only when all four checks are true; otherwise `fail`, listing every failed check in `issues` (labels START/STOP/FRAME/TAG). |

## Summary

- Constructs audited: 3
- Pass: 0
- Fail: 3

## Results

| construct_id | frame_ok | start_ok | stop_ok | tag_ok | overall_status | issues |
| --- | --- | --- | --- | --- | --- | --- |
| c01 | True | True | True | False | fail | TAG |
| c02 | False | True | False | True | fail | STOP;FRAME |
| c03 | False | True | True | False | fail | FRAME;TAG |

## Per-construct detail

### c01

- promoter: `EF1a`
- tag: `C_terminal_FLAG`
- claimed_frame: `in_frame`
- insert_sequence: `ATGGCCGACTACAAAGACGATGACGACAAGTAA`
- insert length: 33 nt (length mod 3 = 0)
- start_ok: True
- stop_ok: True
- frame_ok: True
- tag_ok: False
- overall_status: **fail**
- issues: TAG

Explanation:

- TAG: `C_terminal_FLAG` fusion but the insert ends in stop codon `TAA`, which would terminate translation before the downstream tag.

### c02

- promoter: `CMV`
- tag: `C_terminal_FLAG`
- claimed_frame: `in_frame`
- insert_sequence: `ATGGCCGACTACAAAGACGATGACGACAAGAA`
- insert length: 32 nt (length mod 3 = 2)
- start_ok: True
- stop_ok: False
- frame_ok: False
- tag_ok: True
- overall_status: **fail**
- issues: STOP;FRAME

Explanation:

- STOP: insert length 32 is not divisible by 3, so the insert does not end in an in-frame stop codon.
- FRAME: insert length 32 is not divisible by 3.

### c03

- promoter: `EF1a`
- tag: `C_terminal_FLAG`
- claimed_frame: `out_of_frame`
- insert_sequence: `ATGGCCGACTACAAAGACGATGACGACAAGTGA`
- insert length: 33 nt (length mod 3 = 0)
- start_ok: True
- stop_ok: True
- frame_ok: False
- tag_ok: False
- overall_status: **fail**
- issues: FRAME;TAG

Explanation:

- FRAME: claimed_frame is `out_of_frame`, not `in_frame`.
- TAG: `C_terminal_FLAG` fusion but the insert ends in stop codon `TGA`, which would terminate translation before the downstream tag.

## Scope notes

- The fixture supplies only `construct_id`, `promoter`, `insert_sequence`, `tag`, and `claimed_frame`. The frozen rules define checks over the insert sequence, tag, and claimed frame; `promoter` is recorded but has no audit rule attached.
- No linker sequence and no cloning flags (e.g., restriction sites) are represented in the input, so none were evaluated or inferred, per the rule that only fields represented by the fixture are checked.
- For `C_terminal_*` fusions the rules make `stop_ok` (insert ends in an in-frame stop) and `tag_ok` (no terminal stop before the downstream tag) mutually exclusive, so a C-terminal fusion construct cannot pass both; each check is reported independently as specified.
