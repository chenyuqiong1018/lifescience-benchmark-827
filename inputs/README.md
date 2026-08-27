# Input index

This directory contains inputs for exactly the 25 task IDs listed in `docs/contracts/selected-tasks-v1.tsv`.

For a run, copy the contents of exactly one `inputs/<task-id>/` directory into a fresh workspace's `inputs/` directory. Keep file names and subdirectories unchanged and never edit the repository copy.

The authoritative file inventory, byte size and SHA-256 for every participant-visible input is [`SHA256SUMS.tsv`](SHA256SUMS.tsv). The matching individual task card under `task-card/` contains the frozen exact input list. Those two sources supersede prose summaries elsewhere.

Several task pairs intentionally repeat byte-identical source files:

- the two LS06 tasks share the same two workbooks;
- the two LS07 tasks share counts, sample layout and gene-name mapping;
- two sequence tasks share the same frozen chr9 reference.

These copies are intentional: every task directory must remain independently uploadable. Git deduplicates identical blob content internally.

## Integrity rule

`SHA256SUMS.tsv` must cover every file under the 25 task directories and no stale path. A missing file, extra file, byte-size mismatch or hash mismatch blocks release.

## Adding or changing inputs

Record provenance and redistribution terms in the task card or a task-local README. Do not add secrets, personal data, controlled clinical data, answer-bearing notebooks or files whose license prohibits redistribution. After an intentional change, run:

```bash
python3 scripts/refresh-input-manifest.py
python3 scripts/audit-c0-t1-campaign.py
```

Update the matching standalone task card's input inventory in the same change; the audit fails if a packaged input is absent from the card.
