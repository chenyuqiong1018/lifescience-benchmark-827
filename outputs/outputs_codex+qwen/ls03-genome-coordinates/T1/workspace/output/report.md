# Enhancer–Promoter 3D Distance and Transcription Dynamics

## 1. Objective

Quantify, across single cells and over time, the 3D spatial relationship between an
enhancer and a promoter and the transcriptional state of the locus, and test whether
enhancer–promoter **contact** is temporally associated with **transcription** at a lag.
Per the supplied rule, the contact threshold is **260 nm**, and any lagged association is
interpreted strictly as a **temporal association, never as evidence of causation**.

## 2. Data

- Input: `inputs/single_cell_dynamics_question.csv`
- 150,000 observations = **600 cells × 250 timepoints** (consecutive integer times 0–249 per cell).
- Columns: `cell_id, time, enh_x/y/z, prom_x/y/z, transcription`.
- Coordinates are in nm; `transcription` is binary (0/1). Promoter coordinates are
  fixed within a cell; the enhancer moves over time.

## 3. Methods

All analyses are implemented in `output/analysis.py` (Python 3, pandas/numpy).

1. **Distance.** Euclidean enhancer–promoter distance per observation:
   `d(t) = sqrt((enh_x−prom_x)² + (enh_y−prom_y)² + (enh_z−prom_z)²)` (nm).
2. **Contact.** `contact(t) = 1` iff `d(t) <= 260 nm` (supplied threshold).
3. **Per-cell metrics** (`output/cell_metrics.csv`,
   `cell_id,n_timepoints,mean_distance_nm,contact_fraction,transcription_fraction`):
   fractions are computed over **all 250 rows** of each cell.
4. **Lag analysis** (`output/lag_analysis.csv`, `lag,association,n_observations`):
   for every integer lag in [−20, +20], the **pooled Pearson correlation** across all
   cells between `contact(t)` and `transcription(t+lag)`. Positive lag means contact
   leads the later transcription value. Pairs are formed **only within each cell**
   (never across cell boundaries); for lag L each cell contributes 250 − |L| pairs, so
   `n_observations = 600 × (250 − |L|)`.
5. **Validation.** An independent brute-force implementation that pairs observations
   by explicit time values (instead of positional shifts) reproduced every association
   and observation count in `lag_analysis.csv` (max absolute difference < 5×10⁻⁶).

## 4. Results

### 4.1 Cell-level metrics (600 cells)

| Quantity | Value |
|---|---|
| Mean E–P distance (all frames) | **526.1 nm** (per-cell means 497.7–560.2 nm; median 525.9) |
| Distance distribution | median 536.5 nm, SD 139.1 nm, range 34.6–1071.5 nm |
| Overall contact fraction (d ≤ 260 nm) | **0.0534** (8,004 / 150,000 frames) |
| Per-cell contact fraction | min 0.0120, median 0.0520, max 0.1000 |
| Overall transcription fraction | **0.1539** |
| Per-cell transcription fraction | min 0.1000, median 0.1520, max 0.2240 |
| Mean distance in contact frames | 206.5 nm (vs 544.1 nm in non-contact frames) |

Contact is sparse and transient: the locus spends ~95% of frames beyond 260 nm.

### 4.2 Lagged contact–transcription association

Pooled Pearson correlation of `contact(t)` with `transcription(t+lag)`
(sign convention: positive lag = contact leads transcription):

| lag | association (r) | n_observations |
|---:|---:|---:|
| −20 | −0.000006 | 138,000 |
| −10 | −0.002507 | 144,000 |
| −5 | +0.002508 | 147,000 |
| −2 | −0.002529 | 148,800 |
| **−1** | **+0.077843** | **149,400** |
| 0 | +0.000668 | 150,000 |
| **+1** | **+0.031616** | **149,400** |
| +2 | −0.005052 | 148,800 |
| +5 | −0.001803 | 147,000 |
| +10 | −0.001419 | 144,000 |
| +20 | −0.004519 | 138,000 |

Key findings:

- **Strongest absolute association: lag = −1, r = +0.0778.** With the defined sign
  convention, lag −1 correlates contact at time *t* with transcription at time *t−1*:
  contact tends to occur **one timepoint after** a transcription event
  (equivalently, transcription leads contact by one step).
  Consistent conditional rates: P(contact(t) | transcription(t−1)=1) = **0.0945** vs
  P(contact(t) | transcription(t−1)=0) = **0.0459** (≈2.1-fold enrichment).
- **Secondary association: lag = +1, r = +0.0316.** Contact at *t* is weakly associated
  with transcription at *t+1*: P(transcription(t+1) | contact(t)) = **0.202** vs
  **0.151** when not in contact (≈1.34-fold).
- **No contemporaneous association:** lag 0 gives r ≈ 0.0007; contact and transcription
  do not co-occur more than chance within the same frame.
- Beyond |lag| ≥ 2 all associations are at the noise floor (max |r| = 0.0066).
- Transcription shows **no one-step persistence**
  (P(trans(t+1) | trans(t)=1) = 0.154 ≈ P(trans(t+1) | trans(t)=0) = 0.154;
  mean burst length ≈ 1.18 frames), so the lag −1 peak is not an artifact of
  transcription autocorrelation.
- Across cells, `contact_fraction` and `transcription_fraction` are uncorrelated
  (r = −0.0024): cells that contact more often are not cells that transcribe more often.

### 4.3 Interpretation of the temporal pattern

The only structured associations in the data sit at lags −1 and +1, with the maximum at
**lag −1**: in this observational record, enhancer–promoter contact most strongly follows,
rather than precedes, transcription. The weaker positive association at lag +1 means
contact also slightly raises the probability of transcription in the next frame. Together
the pattern is symmetric-ish around contact frames (transcription at t−1, contact at t,
slightly elevated transcription at t+1) with no same-frame coupling.

## 5. Temporal association is not causation

Per the analysis rule, the strongest absolute association (lag −1, r ≈ 0.078) is reported
as a **temporal association only**. This study is purely observational, and the data
**cannot** establish the following:

1. **Causal direction.** Whether transcription causes subsequent enhancer–promoter
   contact (e.g., transcription machinery, cofactors, or supercoiling recruiting or
   stabilizing the loop), whether contact causes transcription, or whether neither
   causes the other. The sign and location of the correlation peak describe *ordering
   in time*, not *direction of influence*.
2. **Absence of confounding.** A shared upstream driver could produce both signals with
   the observed one-step offset — e.g., transcription-factor binding/coactivator
   recruitment, local chromatin state or accessibility, cell-cycle stage, nuclear
   microenvironment, or global polymer dynamics. None of these are measured here, so
   confounding cannot be ruled out.
3. **Reverse causation.** The data cannot exclude that the apparent "effect" runs in the
   opposite direction of any proposed mechanism.
4. **Causal effect sizes.** A pooled correlation of 0.078 explains <1% of variance; it
   is not an estimate of any causal effect, and no causal effect is identified.
5. **Mechanism or mediation.** Nothing in the data distinguishes direct physical
   coupling from indirect coordination through unobserved intermediates.
6. **That the one-step lag equals a true mechanistic delay.** The sampling interval is
   one frame; dynamics faster than the frame rate are aliased, and the true delay (if
   any causal link exists) could be sub-frame or distributed.

Establishing causation would require **intervention** — e.g., acute enhancer deletion or
silencing, degron-mediated depletion of the bridging factor, forced tethering of the
enhancer and promoter, or live-cell perturbation with matched temporal resolution — none
of which are present in this dataset.

## 6. Limitations

- Binary contact calls depend on the supplied 260 nm threshold; results at other
  thresholds may differ.
- Pooled correlations treat all frames as exchangeable observations; within-cell
  dependence is respected for pairing but not modeled (no significance testing was
  requested or performed).
- Coordinates are observed, not perturbed; measurement noise and (if applicable)
  fixation artifacts are not quantifiable from this file alone.
- One locus, one cell type, one time resolution; generalizability is unknown.

## 7. Deliverables

| File | Description |
|---|---|
| `output/cell_metrics.csv` | 600 rows: `cell_id,n_timepoints,mean_distance_nm,contact_fraction,transcription_fraction` |
| `output/lag_analysis.csv` | 41 rows (lag −20…+20): `lag,association,n_observations` |
| `output/analysis.py` | Reproducible analysis script (`python output/analysis.py` from the workspace root) |
| `output/report.md` | This report |
