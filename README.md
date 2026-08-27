# Life-science evaluation inputs and oracles

This repository contains the participant-visible inputs and deterministic graders for 25 life-science agent evaluation tasks.

The files are frozen from [`Nono111-dot/lifescience`](https://github.com/Nono111-dot/lifescience) commit [`c14f6974804e09bff9fd1b435d9d888582c5a41c`](https://github.com/Nono111-dot/lifescience/tree/c14f6974804e09bff9fd1b435d9d888582c5a41c/docs/inputs).

- Input directories: `inputs/<task-id>/`
- Authoritative inventory and checksums: `inputs/SHA256SUMS.tsv`
- Task count: 25
- Participant-visible files: 89

Keep file names, directory structure, and file contents unchanged. For one evaluation run, copy the contents of exactly one task directory into a fresh workspace's `inputs/` directory.

Provenance and redistribution terms remain governed by the source repository and its task-local documentation.

## Oracle directories

- `oracles/<task-id>/oracle.py`: self-contained grader used for scoring; each task directory contains exactly one Python file.
- `oracles集合/`: archived upstream oracle materials, including shared code, scientific checks, acceptance evidence, gold data, manifests, and tests.

Run a grader with `python oracles/<task-id>/oracle.py --workspace <workspace> --out result.json`. Keep `oracles集合/` private from evaluated agents.
