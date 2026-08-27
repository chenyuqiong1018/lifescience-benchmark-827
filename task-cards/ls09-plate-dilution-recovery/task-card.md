# 任务卡： `ls09-plate-dilution-recovery`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## `ls09-plate-dilution-recovery`

| 字段 | 内容 |
| --- | --- |
| ID | ls09-plate-dilution-recovery |
| Domain / sub-domain | laboratory automation / dilution recovery |
| Level / time | L2, 40 min |
| Anchor / related capabilities | [liquid-handling] |
| 来源思想 | custom / new-fixture |

### Inputs (authoritative packaged inventory)
- `inputs/README.md` — 1,204 bytes
- `inputs/dilution_request.csv` — 94 bytes
- `inputs/pipettes.csv` — 62 bytes
- `inputs/plate_map.csv` — 187 bytes
- `inputs/run_log.csv` — 587 bytes
- `inputs/source_inventory.csv` — 182 bytes

**Total:** 2,316 bytes (0.00 MiB).

### Prompt

> Read every file under `inputs/`. Diagnose the stopped dilution run and generate only the recovery work that remains physically necessary. Use no external data and do not alter inputs. Write `output/root_cause.json` with `failed_well,failure_mode,liquid_moved,completed_wells,recovery_wells`; `failure_mode` may be a concise controlled phrase that states the failed operation and whether it occurred before aspiration. Write `output/recovery_plan.csv` with exactly `step,source,destination,transfer_uL,transfer_pipette,diluent_source,diluent_uL,diluent_pipette,final_concentration,final_volume_uL`; the two pipette fields identify the physical instrument for each distinct liquid movement. Also write rerunnable `output/analysis.py` and `output/report.md`. Enforce the plate map, event ordering, `C_source*V_transfer=C_final*V_final`, source inventory, solvent identity, and frozen pipette ranges. Do not redo completed wells. Abort explicitly rather than inventing a missing transfer.
### Deliverables

structured root cause; one row per required recovery destination; separate pipettes for solute and diluent transfers; report and rerunnable script.

### Hard gates

root cause traceable to the run log; only failed/requested wells are recovered; dilution mass balance and both pipette ranges pass; no source overdraw.

### DeterministicArtifactScore（0–80）

coverage/schema 10; root cause 14, recovery plan 10, concentration/mass balance 10, pipette plus inventory feasibility 6; recover/abort decision 15; report consistency 5; static/rerunnable script 10.

### JudgeScore（0–20）

| Dimension | 0 | 3 | 5 |
| --- | --- | --- | --- |
| Evidence | No task-specific evidence, fabricated/untraceable support, or contradiction with machine-readable artifacts. | At least one traceable task-specific input/result supports the main claim, but coverage, linkage or uncertainty is incomplete. | Every main claim links to the relevant input/result artifact; decisive measurements/rows and evidence-bound uncertainty are explicit. |
| Method | Missing or materially invalid method. | Broadly valid method with incomplete choices, direction, units or limitations. | Correct, auditable method with all consequential choices, direction, units and limitations stated. |
| Restraint | Unsupported causal/clinical/experimental claim. | Mostly calibrated but one important boundary is vague. | Claims stay within the supplied evidence and all important limitations are explicit. |
| Readability | Unusable or internally confusing report. | Understandable with notable ambiguity or clutter. | Concise, internally consistent and easy to audit against formal artifacts. |

The judge records main-claim count, supported-claim count, evidence locations and a short rationale. Judge points cannot repair a deterministic hard-gate failure.

### Ablation（不进入 Prompt）

Skill: C0=`NONE`; T0=`AUTO`; T1=`bioassay_analysis`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

