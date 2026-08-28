All deliverables are complete and independently verified. 

## Summary

**Data:** 600 cells × 250 timepoints (150,000 rows) of enhancer/promoter 3D coordinates + binary transcription.

**Deliverables** (all in `output/`):

| File | Content |
|---|---|
| `cell_metrics.csv` | 600 rows, `cell_id,n_timepoints,mean_distance_nm,contact_fraction,transcription_fraction` — Euclidean distance in nm, contact = d ≤ 260 nm, fractions over all 250 rows per cell |
| `lag_analysis.csv` | 41 rows (lags −20…+20), `lag,association,n_observations` — pooled Pearson correlation between contact(t) and transcription(t+lag), pairs never joined across cell boundaries (n = 600 × (250 − \|lag\|)) |
| `analysis.py` | Deterministic, validated pipeline reproducing both CSVs |
| `report.md` | Full methods, results, and association-vs-causation discussion |

**Key findings:**
- Mean enhancer–promoter distance ≈ 526 nm; contact fraction ≈ 5.3%; transcription fraction ≈ 15.4%.
- Lag profile is flat (|r| ≤ 0.007) except at ±1 step: **strongest association at lag −1 (r ≈ 0.078)** — transcription one step earlier coincides with contact now — and a weaker secondary peak at **lag +1 (r ≈ 0.032)**, where contact leads transcription. Contemporaneous (lag 0) association ≈ 0.
- A within-cell sensitivity check (mean of per-cell correlations) reproduces the pooled profile, ruling out between-cell heterogeneity as the source.

**Association vs. causation:** The report states explicitly that the lag −1 peak is a temporal association only. The observational data cannot establish: the direction of causation; necessity/sufficiency of contact for transcription (requires perturbation experiments); absence of unmeasured confounders (TF binding, chromatin state, cell cycle); the underlying mechanism or true time scale; nor treat temporal precedence as proof — lagged correlation can arise from autocorrelation and common upstream drivers.