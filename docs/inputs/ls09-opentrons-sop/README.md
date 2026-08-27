# Input contract: ls09-opentrons-sop

## Files and roles

- `sop.md`: frozen 24-sample magnetic-bead cleanup workflow and contamination-safe tip policy.
- `instrument.csv`: OT-2, Protocol API 2.16, P300, mount, range, channels, and tip inventory.
- `labware.csv`: deck slots, API load names, module placement, roles, and well capacities.
- `reagent_map.csv`: reagent wells, starting volumes, and dead volumes.
- `sample_map.csv`: 24 sample identifiers, wells, and starting volumes.
- `simulator_contract.json`: evaluator-owned simulator version, platform contract, invocation, and success criteria.

## Provenance and license

This is a locally authored synthetic fixture grounded in the official Opentrons Python Protocol API, instrument, module, and labware documentation. It contains no patient data. The fixture is intended for CC0-1.0 release after scientific review. Opentrons software is Apache-2.0 and is not embedded in this input directory.

## Schema, units, and missing values

CSV identifiers are unique within their declared role. Volumes are microlitres. Empty module capacity is not a numeric zero. The Protocol API level and API load names are exact strings.

## Integrity and readiness boundary

`../SHA256SUMS.tsv` and `docs/oracles/ls09-opentrons-sop/input-sha256.tsv` freeze the same files. The simulator dependency lock is under `docs/environments/opentrons-api-2.16/`. The reference protocol completed three clean Opentrons 7.1.0 simulations on the campaign host; the formal runner prepends the same pre-provisioned environment for both C0 and T1. The original Linux lock remains the release environment authority, while host-platform use is recorded as a campaign deviation.
