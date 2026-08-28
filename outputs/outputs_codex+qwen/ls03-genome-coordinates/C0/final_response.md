All deliverables are complete and verified. Summary of the analysis:

## Deliverables (in `output/`)

| File | Content |
|---|---|
| `cell_metrics.csv` | 600 cells × `cell_id,n_timepoints,mean_distance_nm,contact_fraction,transcription_fraction` (fractions over all 250 rows/cell, contact = distance ≤ 260 nm) |
| `lag_analysis.csv` | 41 rows, every integer lag −20…+20: `lag,association,n_observations` (pooled within-cell Pearson r of contact at *t* vs transcription at *t*+lag) |
| `analysis.py` | Self-contained, deterministic (seeded bootstrap) script that regenerates all outputs |
| `report.md` | Full report with methods, results, and association-vs-causation analysis |

Both CSVs were cross-validated against an independent implementation (max deviation ~5×10⁻⁷, i.e., 6-decimal rounding only).

## Key findings

- **600 cells × 250 timepoints**; mean E–P distance 526 nm, contact is rare (5.3% of timepoints at the supplied **260 nm** threshold), transcription active in 15.4%.
- **Contemporaneous coupling ≈ 0**: r = +0.0007 at lag 0 — current contact says nothing about current transcription.
- **Strongest association at lag −1**: r = +0.0778 (block-bootstrap 95% CI [+0.069, +0.086]; 100% of bootstrap replicates select lag −1). Transcription at *t*−1 precedes contact at *t*: P(tx=1 at *t*−1 | contact) = 0.272 vs 0.154 baseline. A weaker mirror association exists at lag +1 (r = +0.032). Neither signal is autocorrelated at lag 1, so this is genuine adjacency structure, not slow-dynamics leakage.
- **Between-cell** contact vs transcription rates are uncorrelated (r = −0.002), so the signal is purely within-cell temporal structure. The lag −1 peak is stable across thresholds 220–320 nm.

## Causation caveat (detailed in report §4)

The peak is a **temporal association only**: with no intervention, the data cannot establish that contact causes transcription (or vice versa), the direction of influence, that the link is direct rather than driven by a common upstream process (TF binding, chromatin remodeling, cell cycle) with differential delays, or that lag −1 is a true mechanistic delay. Causal claims would require perturbations (optogenetic tethering, enhancer deletion/inversion, degron depletion of bridging factors) against matched controls.