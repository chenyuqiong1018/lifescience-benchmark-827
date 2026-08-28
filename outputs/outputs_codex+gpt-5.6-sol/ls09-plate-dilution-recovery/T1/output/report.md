# Dilution-run recovery report

## Root cause

The run stopped at well **B2** because tip pickup failed before aspiration. The event log explicitly states that no liquid moved, so the failed well remains physically unfilled. This is a categorical mechanical failure, not evidence of a transferred-volume measurement error.

Completed final-dispense wells, in numeric event order: **A1, A2, A3, B1**. They are excluded from recovery. The only wells still requiring work are **B2, B3**.

## Recovery

For each of B2, B3, transfer **2 uL** from **source:A2** with **P20**, then add **98 uL** from **diluent:R1** with **P300**. Each well finishes at **100 uL of 0.5 uM CMPD_B**. No completed well is repeated.

The concentration check is `25 uM * 2 uL = 0.5 uM * 100 uL = 50 uM*uL`. Frozen ranges select exactly one physical instrument for each movement: P20 for 2 uL and P300 for 98 uL. Total remaining demand is **4 uL** from source:A2 and **196 uL** from diluent:R1, below the recorded current inventories of 1998 uL and 5000 uL. The recorded source solvent (`DMSO_0.5pct_in_media`) and diluent solvent (`media`) are preserved without substitution.

## Skill-assisted audit

The controlled protocol-generation result agreed with the two fixed liquid movements and introduced no extra operation. The code-execution capability returned the requested exact-arithmetic snippet but did not execute it, so this rerunnable script independently performs the Decimal calculation and all validations. Measurement-error guidance was applied by treating the pre-aspiration tip-pickup failure as a discrete run-state event rather than inventing an uncertainty correction. Bioassay database lookups and unrelated unit-conversion endpoints were not used because the task forbids external data and already fixes concentration and volume units.
