# Enhancer–Promoter 3D Distance and Transcription Dynamics

## 1. Objective and data

Quantify, across cells and over time, how the 3D distance between an enhancer
and a promoter relates to transcriptional activity, and characterize the
temporal (lagged) association between enhancer–promoter **contact** and
**transcription**.

- Input: `inputs/single_cell_dynamics_question.csv`
  - 150,000 rows = **600 cells × 250 timepoints** per cell.
  - Columns: `cell_id, time, enh_x/y/z, prom_x/y/z, transcription`.
  - Coordinates are in nm; `transcription` is binary (0/1).
- Deliverables: `output/cell_metrics.csv`, `output/lag_analysis.csv`,
  `output/analysis.py`, and this report.

## 2. Methods

**Distance and contact.** For each row, the Euclidean enhancer–promoter
distance is computed from the supplied x/y/z coordinates (nm):
`d = sqrt((enh_x−prom_x)² + (enh_y−prom_y)² + (enh_z−prom_z)²)`.
A **contact** is defined with the supplied threshold: `contact = (d ≤ 260 nm)`.

**Per-cell metrics** (`output/cell_metrics.csv`, one row per cell, columns
`cell_id, n_timepoints, mean_distance_nm, contact_fraction,
transcription_fraction`): fractions are computed over **all 250 rows** of the
cell.

**Lag analysis** (`output/lag_analysis.csv`, columns
`lag, association, n_observations`): for every integer lag from **−20 to +20**,
the **pooled Pearson correlation across all cells** between `contact(t)` and
`transcription(t + lag)` is computed. Positive lag means contact **leads** the
later transcription value. Pairs are formed strictly within each cell's own
time series (no joins across cell boundaries), giving
`n_observations = 600 × (250 − |lag|)`.

Per the frozen rule, the strongest absolute association is interpreted as a
**temporal association only, never proof of causation** (Section 5).

## 3. Results — per-cell distance and contact summary

| Quantity | Value |
|---|---|
| Cells | 600 (250 timepoints each) |
| Mean enhancer–promoter distance | 526.07 nm (cell means: 497.7–560.2 nm) |
| Contact fraction (d ≤ 260 nm), pooled | 0.0534 |
| Transcription fraction, pooled | 0.1539 |
| Contact fraction across cells | 0.012–0.100 (mean 0.053) |
| Transcription fraction across cells | 0.100–0.224 (mean 0.154) |

Contacts are rare relative to the time series (≈5.3% of frames) and the mean
distance (~526 nm) sits well above the 260 nm threshold, so contact events are
brief, intermittent excursions.

## 4. Results — lagged contact–transcription association

Full profile in `output/lag_analysis.csv`. Key features:

| lag | association (pooled r) | n_observations | reading |
|---:|---:|---:|---|
| −1 | **+0.0778** | 149,400 | transcription at t−1 associates with contact at t (transcription leads contact by 1 step) |
| +1 | +0.0316 | 149,400 | contact at t associates with transcription at t+1 (contact leads transcription by 1 step) |
| 0 | +0.0007 | 150,000 | essentially no contemporaneous association |
| all other lags | \|r\| ≤ 0.0066 | — | indistinguishable from noise |

- **Strongest absolute association: lag = −1 (r ≈ 0.078).** Under the frozen
  sign convention (positive lag = contact leads), a *negative* lag means the
  transcription value comes from the *earlier* timepoint: transcription one
  step before is weakly but consistently more likely when the enhancer and
  promoter are in contact now.
- A secondary, weaker positive peak appears at **lag = +1** (r ≈ 0.032), the
  direction in which contact precedes transcription. Both peaks are confined to
  ±1 time step; the profile is flat elsewhere.
- **Sensitivity check (not part of the frozen output):** averaging per-cell
  within-cell Pearson correlations instead of pooling gives an almost identical
  profile (peak at lag −1, mean r ≈ 0.074; lag +1 ≈ 0.033). The pooled signal
  is therefore not an artifact of between-cell heterogeneity (cells differing
  in mean contact/transcription rates).
- Effect sizes are small: even the strongest association explains well under
  1% of variance (r² ≈ 0.006). With ~150,000 non-independent observations,
  p-values would be misleadingly extreme; effect size and the flatness of the
  rest of the profile are the meaningful evidence here.

## 5. Temporal association ≠ causation

The lag structure is a **temporal association only**. Concretely, this
observational dataset **cannot** establish:

1. **Direction of causation.** The strongest peak (lag −1) is compatible with
   multiple stories: (a) transcriptional activity (or machinery recruited with
   it) helping to establish or stabilize enhancer–promoter contact one time
   step later; (b) contact enabling subsequent transcription (consistent with
   the weaker lag +1 peak); (c) no direct link at all, with a common upstream
   driver (e.g., transcription-factor binding, chromatin opening, bursting
   kinetics, cell-cycle stage) producing both signals with different temporal
   signatures. The data cannot choose among these.
2. **Necessity or sufficiency of contact for transcription** (or vice versa).
   Only interventions — e.g., acute enhancer deletion/silencing, degron-based
   disruption of looping factors, or optogenetic forcing of contact — can test
   whether contact is required for, or capable of triggering, transcription.
3. **Absence of confounding.** Unmeasured variables (TF occupancy, Pol II
   loading, local chromatin state, nuclear position, cell-cycle phase) may
   drive both distance and transcription. The within-cell design avoids joining
   observations across cells, but cannot remove shared upstream drivers inside a
   cell.
4. **Mechanism or time scale.** A ±1-time-step peak only says the two binary
   signals co-vary across adjacent frames; it does not identify the molecular
   sequence of events, and the discrete time resolution may alias faster or
   slower dynamics.
5. **Temporal precedence as proof.** Lagged correlation can arise from
   autocorrelation within each series and from a common cause; "contact
   precedes transcription" (or the reverse) in a cross-correlation sense is not
   evidence that one event makes the other happen.

**Bottom line:** the only defensible claim from these data is descriptive —
contact and transcription are weakly temporally coupled at ±1 time step
(strongest at lag −1, pooled r ≈ 0.078), with no detectable contemporaneous
association and none at longer lags. Establishing whether enhancer–promoter
contact causes, follows, or merely accompanies transcription requires
perturbation experiments, not observational time series.

## 6. Reproducibility

Run from the workspace root:

```powershell
python output/analysis.py
```

The script validates input shape (600 cells × 250 unique, increasing
timepoints, no missing values), recomputes both CSVs deterministically, and
prints the summary statistics quoted above.
