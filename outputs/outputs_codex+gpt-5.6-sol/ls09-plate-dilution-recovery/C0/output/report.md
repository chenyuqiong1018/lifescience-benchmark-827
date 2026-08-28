# Dilution-run recovery

## Diagnosis

Numeric event order shows completed final dispenses to A1, A2, A3, B1. Event 7 then records a tip-pickup failure at B2 before aspiration and explicitly states that no liquid moved. The run stopped there. Consequently, only B2, B3 remain; completed wells are not repeated.

## Recovery calculation

Both remaining wells request 0.5 uM in 100 uL for CMPD_B. The completed intermediate preparation and inventory agree on source:A2 at 25 uM in `DMSO_0.5pct_in_media`. For each well, C_source x V_transfer = C_final x V_final gives 25 uM x 2 uL = 0.5 uM x 100 uL. Add 2 uL intermediate with P20, then 98 uL from diluent:R1 (`media`) with P300. These instruments are unique under the frozen ranges/increments: P20 2-20 uL by 0.1 uL; P300 20-300 uL by 1 uL.

Total new demand is 4 uL from source:A2 and 196 uL from diluent:R1, below the recorded available volumes of 1998 uL and 5000 uL. No intermediate preparation or completed final dispense is repeated. `analysis.py` aborts on absent/ambiguous sources, unexplained liquid movement, inconsistent requests, invalid or ambiguous pipette choice, equation failure, solvent/source mismatch, or insufficient inventory rather than inventing work.
