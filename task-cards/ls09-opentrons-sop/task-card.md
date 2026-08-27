# 任务卡： `ls09-opentrons-sop`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## `ls09-opentrons-sop`

| 字段 | 内容 |
| --- | --- |
| ID | ls09-opentrons-sop |
| Domain / sub-domain | laboratory automation / liquid handling |
| Level / time | L2, 45 min |
| Anchor / related capabilities | [protocol-planning] |
| 来源思想 | custom / new-fixture |

### Inputs (authoritative packaged inventory)
- `inputs/README.md` — 1,699 bytes
- `inputs/instrument.csv` — 123 bytes
- `inputs/labware.csv` — 329 bytes
- `inputs/reagent_map.csv` — 200 bytes
- `inputs/sample_map.csv` — 274 bytes
- `inputs/simulator_contract.json` — 1,108 bytes
- `inputs/sop.md` — 1,481 bytes

**Total:** 5,214 bytes (0.00 MiB).

### Prompt

> Read every file under `inputs/`. Translate the frozen 24-sample magnetic-bead cleanup SOP into an auditable OT-2 Opentrons protocol. Use no external data and do not alter inputs. Write `output/protocol.py`, `output/transfer_plan.csv` with exactly the columns `step,source,destination,volume_uL,pipette,tip_policy`, `output/simulation.txt`, and `output/report.md`. `transfer_plan.csv` is a net liquid-transfer stage table, not a command log: write exactly one row per SOP net-transfer stage per sample (`lysis`, `beads`, `supernatant`, `wash1_add`, `wash1_remove`, `wash2_add`, `wash2_remove`, `elution`), for 8 × 24 = 192 rows. Do not add one row per mix stroke, aspirate command, dispense command, delay, magnet action, or tip action; represent those operations in `protocol.py` and summarize them in `report.md`. Identify each row as `<stage>:<well>` and use the frozen role/well identifiers. Respect the declared deck, Magnetic Module compatibility, API level, pipette range, well capacity, reagent dead volumes, liquid balance, and tip policy. Run the evaluator-pinned simulator using the supplied invocation and record its unedited outcome in `simulation.txt`; if that simulator or invocation is unavailable, record the failure and abort rather than claiming success.
### Deliverables

static Opentrons protocol; 192-row net-transfer plan; verbatim simulator record; report. All wells/volumes/pipettes/tip policies must be explicit.

### Hard gates

exact net-transfer contract and liquid balance; valid deck/labware/wells/pipette range and contamination-safe tip policy; static protocol contract; pinned simulation success.

### DeterministicArtifactScore（0–80）

coverage/schema 10; transfer contract, balance, pipette and tip policy 40; protocol/simulation decision 15; report consistency 5; static protocol plus isolated pinned simulation 10.

### JudgeScore（0–20）

| Dimension | 0 | 3 | 5 |
| --- | --- | --- | --- |
| Evidence | No task-specific evidence, fabricated/untraceable support, or contradiction with machine-readable artifacts. | At least one traceable task-specific input/result supports the main claim, but coverage, linkage or uncertainty is incomplete. | Every main claim links to the relevant input/result artifact; decisive measurements/rows and evidence-bound uncertainty are explicit. |
| Method | Missing or materially invalid method. | Broadly valid method with incomplete choices, direction, units or limitations. | Correct, auditable method with all consequential choices, direction, units and limitations stated. |
| Restraint | Unsupported causal/clinical/experimental claim. | Mostly calibrated but one important boundary is vague. | Claims stay within the supplied evidence and all important limitations are explicit. |
| Readability | Unusable or internally confusing report. | Understandable with notable ambiguity or clutter. | Concise, internally consistent and easy to audit against formal artifacts. |

The judge records main-claim count, supported-claim count, evidence locations and a short rationale. Judge points cannot repair a deterministic hard-gate failure.
- Readiness closure: the reference protocol completed three real Opentrons 7.1.0 simulations and wrong/empty/legacy-tip controls fail; both C0 and T1 receive the same pre-provisioned simulator. The macOS-arm64 versus production Linux-x86_64 platform difference is retained as a disclosed campaign deviation.

### Ablation（不进入 Prompt）

Skill: C0=`NONE`; T0=`AUTO`; T1=`synthetic_biology_design`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

