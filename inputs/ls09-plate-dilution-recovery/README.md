# Input contract: ls09-plate-dilution-recovery

## Files and roles

- `dilution_request.csv`: compounds, stock concentrations, requested final concentrations/volumes, and replicate counts.
- `plate_map.csv`: requested wells and replicate mapping.
- `run_log.csv`: completed operations and the failed-before-aspirate event.
- `source_inventory.csv`: intermediate sources, concentrations, available volumes, and solvent identity.
- `pipettes.csv`: P20/P300 usable ranges and increments.

## Provenance and license

This is a locally authored, source-supported synthetic fixture. It contains no patient data and is intended for CC0-1.0 release after scientific review.

## Schema, units, and missing values

Concentrations are micromolar unless a column explicitly says millimolar; volumes are microlitres. Empty compound identifiers for pure diluent remain empty. Event order is the numeric `event_id` order and must not be inferred from row sorting by well.

## Integrity

`../SHA256SUMS.tsv` and `docs/oracles/ls09-plate-dilution-recovery/input-sha256.tsv` freeze the same five data files. The accepted checker requires the reference to pass 3/3 and empty, wrong, and ambiguous-pipette controls to fail.
