# Stopped dilution run: recovery-only plan

## Root cause classification

The numeric event sequence records final dispenses to A1, A2, A3, B1, followed by a tip-pickup failure at B2 before aspiration. The log explicitly says no liquid moved. This is a discrete pre-transfer mechanical failure, not a noisy volume measurement: there is no partial dose to estimate or correct. The only physically necessary wells are B2, B3; completed wells remain untouched.

## Recovery

Both remaining wells are CMPD_B at 0.5 uM in 100 uL. The completed preparation record and inventory identify `source:A2` as the unique 25 uM intermediate in `DMSO_0.5pct_in_media`. Exact dilution gives 25 x 2 = 0.5 x 100, so each well receives 2 uL from `source:A2` with P20 and 98 uL from `diluent:R1` (`media`) with P300. The frozen ranges and increments make each pipette choice unique.

The recovery consumes 4 uL intermediate and 196 uL diluent versus 1998 uL and 5000 uL available. `analysis.py` validates event order, request/plate consistency, source identity and concentration, solvent-preserving source selection, the dilution equation, pipette range/increment uniqueness, and inventory. It aborts rather than estimating an unrecorded movement.

The freshly installed `measurement-error-analysis` skill sharpened the classification boundary: because failure occurred before aspiration, uncertainty propagation is inappropriate; the correct action is a complete transfer only for the uncompleted wells.
