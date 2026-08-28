# OT-2 magnetic-bead cleanup translation

## Status

**Not execution-ready.** The evaluator-specified invocation `python -m opentrons.simulate output/protocol.py` was run from the experiment-arm root. It exited with code 1 because the active Python environment has no `opentrons` module. `simulation.txt` records the command, empty stdout, unedited stderr, and exit code. In accordance with the frozen failure policy, no package was installed, no alternative simulator was substituted, and no successful simulation is claimed.

## Protocol and transfer-plan audit

`protocol.py` targets OT-2 Protocol API 2.16, loads the P300 single-channel pipette on the right mount, two 96-tip racks in slots 4 and 7, the reagent reservoir in slot 1, liquid waste in slot 6, and the deep-well processing plate on a Magnetic Module Gen2 in slot 5. Samples follow the supplied column-major order A1–D6.

The protocol performs batch lysis, bead binding, magnetic separation, supernatant removal, two washes while engaged, air drying, disengagement, and elution. Each wash tip remains attached during its 30-second dwell and is reused only for removal from the same sample. Other stages use one fresh tip per sample. The run therefore requires 144 of 192 tips.

`transfer_plan.csv` is a net-transfer table with exactly the required columns and 192 unique rows: 24 rows for each of `lysis`, `beads`, `supernatant`, `wash1_add`, `wash1_remove`, `wash2_add`, `wash2_remove`, and `elution`. Mix strokes, waits, magnet operations, and tip commands appear only in the protocol.

Static balances pass: lysis uses 1,920 µL leaving 580 µL (300 µL dead); beads use 2,880 µL leaving 620 µL (300 µL dead); wash uses 8,640 µL leaving 360 µL (300 µL dead); elution uses 960 µL leaving 540 µL (300 µL dead). Waste is 14,640 µL against 195,000 µL capacity. Every commanded liquid volume is within the declared 20–300 µL P300 range, and the maximum sample-well liquid load is 250 µL against 2,000 µL capacity.
