# OT-2 magnetic-bead cleanup translation

## Status

**Not execution-ready.** From this experiment-arm root, the frozen invocation `python -m opentrons.simulate output/protocol.py` exited 1 because `opentrons` is absent from the active Python environment. `simulation.txt` preserves the empty stdout, exact stderr, and exit code. The simulator was not installed or replaced, so no successful simulation is claimed.

## Frozen-procedure translation

`protocol.py` declares OT-2/API 2.16 and the supplied deck: reagent reservoir in slot 1, tip racks in 4 and 7, Magnetic Module Gen2 with the deep-well plate in 5, and liquid waste in 6. The P300 single-channel pipette is on the right mount. Samples are processed in supplied column-major order A1–D6.

The code preflights reagent dead volumes, 144-tip demand, P300 range, 2,000 µL well capacity, and 195,000 µL waste capacity. It batch-executes lysis/mixing, the 5-minute incubation, bead-source resuspension and sample mixing, the second 5-minute incubation, magnet engagement and 7-minute separation, supernatant removal, two engaged washes, 2-minute air dry, magnet disengagement, and elution. Each wash uses one fresh tip that stays attached for its 30-second dwell and is reused only to remove liquid from that same sample.

`transfer_plan.csv` has exactly the requested six columns and 192 unique `<stage>:<well>` rows: eight net-transfer stages for each of 24 wells. It excludes all mix strokes, delays, magnet actions, and tip commands.

Static balances pass. Remaining volumes after planned consumption are lysis 580 µL, beads 620 µL, wash 360 µL, and elution 540 µL, each at or above its 300 µL dead-volume boundary. Planned waste is 14,640 µL and maximum processing-well load is 250 µL.

## Skill boundary

The freshly installed `protocol-generation-from-description` skill was used to check procedural coverage. Its response introduced non-frozen labware, temperature, centrifugation, storage, and tip suggestions; every such addition was rejected. Only the supplied inputs govern these artifacts.
