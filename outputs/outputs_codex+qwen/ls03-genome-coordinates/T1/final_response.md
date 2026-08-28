All deliverables are complete, validated against the frozen rules, and independently cross-checked.

## Summary

**Data:** 150,000 observations = 600 cells × 250 consecutive timepoints (0–249). Enhancer coordinates move over time; promoter coordinates are fixed per cell.

**Deliverables in `output/`:**
- **`cell_metrics.csv`** — 600 rows: `cell_id,n_timepoints,mean_distance_nm,contact_fraction,transcription_fraction`. Euclidean distance in nm; contact = distance ≤ 260 nm; fractions over all 250 rows per cell.
- **`lag_analysis.csv`** — 41 rows (every integer lag −20…+20): `lag,association,n_observations`. Pooled Pearson correlation of contact at *t* vs transcription at *t+lag*, paired only within cells (never across cell boundaries).
- **`analysis.py`** — reproducible script (`python output/analysis.py`).
- **`report.md`** — methods, results, and the association-vs-causation discussion.

**Validation:** an independent brute-force pairing (by explicit time values) reproduced every association and observation count; cell-level metrics spot-checked for cells 0/299/599; structure asserts all pass.

**Key findings:**
- Mean E–P distance ≈ 526 nm; contact (≤260 nm) is sparse (5.3% of frames); transcription active in 15.4% of frames.
- Lag 0 association ≈ 0 (no contemporaneous coupling).
- **Strongest |association| at lag −1 (r = +0.078):** contact at *t* tracks transcription at *t−1* — contact tends to *follow* transcription by one step (P(contact | prior transcription) = 0.095 vs 0.046).
- Weaker secondary signal at lag +1 (r = +0.032); all |lag| ≥ 2 are at noise (max |r| = 0.0066).
- Transcription has no one-step persistence, so the lag −1 peak is not autocorrelation leakage.

**Association ≠ causation (as required):** the report states explicitly that this observational record cannot establish causal direction, cannot rule out confounding (TFs, chromatin state, cell cycle, polymer dynamics) or reverse causation, yields no causal effect size or mechanism, and cannot confirm the one-step lag is a true mechanistic delay. Establishing causation would require intervention (enhancer deletion/silencing, forced tethering, factor depletion), which the dataset lacks.