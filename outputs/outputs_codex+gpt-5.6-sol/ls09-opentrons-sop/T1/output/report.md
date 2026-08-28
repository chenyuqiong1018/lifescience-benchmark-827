# Audited OT-2 magnetic-bead cleanup

## Execution status

**Not execution-ready.** The evaluator-pinned command `python -m opentrons.simulate output/protocol.py` was run unchanged from this experiment-arm root. It returned exit code 1 because `opentrons` is unavailable in the active Python environment. `simulation.txt` contains the command, empty stdout, unedited stderr, and exit code. No dependency was installed, no alternate invocation was used, and no simulation success is claimed.

## Protocol coverage

`protocol.py` declares Protocol API 2.16 and uses only the frozen OT-2 resources: P300 single-channel on the right mount, two 96-tip racks in slots 4 and 7, reagent reservoir in slot 1, Magnetic Module Gen2 plus deep-well processing plate in slot 5, and liquid waste in slot 6. The sample tuple exactly follows A1, B1, C1, D1 through A6, B6, C6, D6.

The operation order is batch lysis/addition/mixing; 5-minute incubation; bead-source resuspension, bead addition and mixing; 5-minute incubation; magnet engagement and 7-minute wait; supernatant removal; two engaged wash cycles; 2-minute air dry; disengagement; and elution/mixing. Each wash loop picks up one fresh tip, adds wash, holds that tip during the 30-second dwell, removes wash from the same sample, and discards it. The other four transfer stages each use a separate fresh tip per sample, totaling 144 tips.

## Transfer and balance audit

`transfer_plan.csv` contains exactly the requested columns and 192 unique `<stage>:<well>` rows: 24 for each of eight net-transfer stages. It deliberately excludes mixes, delays, magnet actions, and tip commands.

Preflight arithmetic is embedded in the protocol. Consumption/remaining volumes are lysis 1,920/580 µL, beads 2,880/620 µL, wash 8,640/360 µL, and elution 960/540 µL; all remain above the 300 µL dead-volume boundary. Waste is 14,640 µL (14.64 cm³) versus 195,000 µL capacity. Peak processing-well load is 250 µL versus 2,000 µL. All 30–250 µL liquid actions are within the declared 20–300 µL pipette range.

## Skill audit boundary

The prescribed skills were freshly installed and opened. Protocol generation and executable-JSON conversion were used only for coverage checking; their repeated per-sample incubation and any non-frozen suggestions were rejected. Unit arithmetic and code checks were reproduced locally because the exposed code endpoint only echoed the audit snippet and no compatible unit-conversion endpoint was callable. Synthetic-biology sequence tools and external dataset/literature searches were irrelevant and not invoked.
