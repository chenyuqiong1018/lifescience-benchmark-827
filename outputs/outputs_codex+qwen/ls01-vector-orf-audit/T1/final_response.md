Audit complete. All three deliverables are in `output/`:

## Results — 0 of 3 constructs pass

| construct_id | frame_ok | start_ok | stop_ok | tag_ok | overall_status | issues |
|---|---|---|---|---|---|---|
| c01 | true | true | true | false | fail | TAG |
| c02 | false | true | false | true | fail | STOP;FRAME |
| c03 | false | true | true | false | fail | FRAME;TAG |

**Per-construct findings** (all inserts are the same DYKDDDDK-like ORF with different 3′ ends):
- **c01** (33 bp, `in_frame`): starts with ATG, ends in in-frame `TAA` — but that terminal stop sits *before* the C-terminal FLAG tag, so `tag_ok` fails.
- **c02** (32 bp): length not divisible by 3 → no in-frame terminal codon (ends `GAA`), so `stop_ok` and `frame_ok` both fail; no terminal stop, so `tag_ok` passes.
- **c03** (33 bp, claimed `out_of_frame`): ends in in-frame `TGA`, but `claimed_frame ≠ in_frame` fails `frame_ok`, and the terminal stop before the tag fails `tag_ok`.

**Method notes:**
- Implemented the frozen rule from `inputs/AUDIT_RULE.md` exactly: `stop_ok` requires the terminal triplet to be in-frame (length % 3 == 0); `tag_ok` fails closed for non-`C_terminal_*` tags; unsupported features were never inferred (promoter column unused, no vector/tag DNA assumed).
- The inherent tension in the rule (C-terminal fusions can't satisfy `stop_ok` and `tag_ok` simultaneously) is reported as-is per the frozen rule.

**Deliverables:**
- `output/construct_audit.csv` — required columns, verified against hand computation
- `output/analysis.py` — deterministic, path-robust script that regenerates the CSV from `inputs/constructs.csv`
- `output/report.md` — full per-construct rationale and summary