# Enhancer-Promoter 3D Distance and Transcription Dynamics

## 1. Data and scope

- Source: `inputs/single_cell_dynamics_question.csv`
- Design: **600 cells**, each observed for **250 consecutive time points** (t = 0..249);
  150,000 observations in total. One enhancer and one promoter locus per cell.
- Coordinates are supplied in **nm**. The promoter position is static within each cell;
  the enhancer position moves over time, so the enhancer-promoter distance is driven
  entirely by enhancer motion.
- `transcription` is binary (0/1) per time point.
- Contact threshold: **260 nm (supplied)** - a locus pair is "in contact"
  when its Euclidean distance is <= 260 nm.

## 2. Methods

**Distance.** Euclidean distance between enhancer and promoter coordinates (nm) at each
time point.

**Per-cell metrics** (`output/cell_metrics.csv`): `n_timepoints`, mean distance,
contact fraction (share of the 250 time points with distance <= 260 nm), and
transcription fraction (share of the 250 time points with transcription = 1). All
fractions use all 250 rows per cell.

**Lag-resolved association** (`output/lag_analysis.csv`): for every integer lag
L in [-20, +20], the pair (contact at time t, transcription at time t + L) is formed
**within each cell only** (no observation is ever joined across cell boundaries; a
positive lag means contact leads the later transcription value). All pooled pairs are
then correlated with a single Pearson correlation (for two binary variables this is the
phi coefficient). `n_observations` is the number of pooled within-cell pairs
(600 x (250 - |L|)).

**Uncertainty.** Because observations within a cell are dependent, the naive Pearson
p-value (which treats all pooled pairs as independent) is reported only for reference.
Primary inference uses a **cell-level block bootstrap**: whole cells are resampled with
replacement (500 replicates, seed 20260828), the pooled correlation is recomputed per
replicate, and 2.5/97.5 percentiles form 95% CIs. This respects the "no joins across
cell boundaries" rule.

## 3. Results

### 3.1 Global and per-cell summary

| quantity | value |
|---|---|
| mean enhancer-promoter distance | 526.1 nm |
| median distance | 536.5 nm |
| distance range | 34.6 - 1071.5 nm |
| pooled contact fraction (<= 260 nm) | 0.0534 |
| pooled transcription fraction | 0.1539 |

Contact is brief/rare relative to the full trajectory: a cell is in contact in only a
few percent of time points on average (median contact fraction
0.052, range 0.012-0.100).

Per-cell distributions (across 600 cells):

| metric | min | p25 | median | p75 | max |
|---|---:|---:|---:|---:|---:|
| mean_distance_nm | 497.683 | 520.440 | 525.905 | 532.011 | 560.246 |
| contact_fraction | 0.012 | 0.044 | 0.052 | 0.064 | 0.100 |
| transcription_fraction | 0.100 | 0.140 | 0.152 | 0.168 | 0.224 |

Across cells, contact fraction and transcription fraction are essentially uncorrelated (Pearson r = -0.0024, p = 0.953, n = 600 cells): cells that spend more time in contact do not have systematically higher (or lower) overall transcription levels. This between-cell null is informative in its own right: the lag-resolved association below is within-cell temporal structure, not an artifact of some cells being both more contact-prone and more transcriptionally active.

### 3.2 Lag-resolved contact-transcription association

Full profile (pooled Pearson r of contact at t vs transcription at t + lag). The
bootstrap CI is from the cell-level block bootstrap; the peak |r| lag is in bold.

| lag | association (pooled r) | n_observations | bootstrap 95% CI |
|----:|-----------------------:|---------------:|------------------|
| -20 | -0.000006 | 138000 | [-0.005378, +0.005704] |
| -19 | -0.001956 | 138600 | [-0.007122, +0.003544] |
| -18 | -0.000493 | 139200 | [-0.005492, +0.004229] |
| -17 | +0.001002 | 139800 | [-0.003790, +0.006107] |
| -16 | -0.001828 | 140400 | [-0.006828, +0.002971] |
| -15 | -0.000573 | 141000 | [-0.005370, +0.004547] |
| -14 | -0.003980 | 141600 | [-0.008907, +0.001648] |
| -13 | -0.000647 | 142200 | [-0.005772, +0.004863] |
| -12 | +0.003623 | 142800 | [-0.002045, +0.009433] |
| -11 | -0.001560 | 143400 | [-0.006516, +0.003518] |
| -10 | -0.002507 | 144000 | [-0.006624, +0.002459] |
| -9 | -0.006580 | 144600 | [-0.011390, -0.001216] |
| -8 | +0.003870 | 145200 | [-0.001682, +0.008687] |
| -7 | -0.000921 | 145800 | [-0.005588, +0.004240] |
| -6 | -0.001430 | 146400 | [-0.006627, +0.003940] |
| -5 | +0.002508 | 147000 | [-0.001904, +0.007158] |
| -4 | -0.000400 | 147600 | [-0.005162, +0.004356] |
| -3 | +0.004165 | 148200 | [-0.001381, +0.009453] |
| -2 | -0.002529 | 148800 | [-0.007507, +0.002524] |
| -1 | **+0.077843** | 149400 | [+0.069287, +0.086326] |
| +0 | +0.000668 | 150000 | [-0.004395, +0.005937] |
| +1 | +0.031616 | 149400 | [+0.025230, +0.037210] |
| +2 | -0.005052 | 148800 | [-0.009842, +0.000728] |
| +3 | -0.000073 | 148200 | [-0.005382, +0.005116] |
| +4 | -0.000306 | 147600 | [-0.006082, +0.004879] |
| +5 | -0.001803 | 147000 | [-0.006961, +0.003346] |
| +6 | -0.002144 | 146400 | [-0.007140, +0.002926] |
| +7 | +0.000941 | 145800 | [-0.003921, +0.005998] |
| +8 | -0.001233 | 145200 | [-0.006357, +0.003828] |
| +9 | -0.001469 | 144600 | [-0.006755, +0.003138] |
| +10 | -0.001419 | 144000 | [-0.007646, +0.003703] |
| +11 | +0.001703 | 143400 | [-0.003509, +0.006643] |
| +12 | -0.004338 | 142800 | [-0.009742, +0.000444] |
| +13 | -0.001235 | 142200 | [-0.006269, +0.003517] |
| +14 | -0.001397 | 141600 | [-0.006398, +0.003888] |
| +15 | -0.006582 | 141000 | [-0.012051, -0.000996] |
| +16 | -0.001798 | 140400 | [-0.007354, +0.003173] |
| +17 | -0.001192 | 139800 | [-0.005849, +0.004182] |
| +18 | -0.002676 | 139200 | [-0.007701, +0.002702] |
| +19 | -0.002253 | 138600 | [-0.007150, +0.003215] |
| +20 | -0.004519 | 138000 | [-0.009050, +0.000936] |

ASCII profile (each '#' ~ 0.0025 in |r|):

```
 -20                                #                                -0.0000
 -19                               ##                                -0.0020
 -18                                #                                -0.0005
 -17                                #                                +0.0010
 -16                               ##                                -0.0018
 -15                                #                                -0.0006
 -14                              ###                                -0.0040
 -13                                #                                -0.0006
 -12                                ##                               +0.0036
 -11                               ##                                -0.0016
 -10                               ##                                -0.0025
  -9                             ####                                -0.0066
  -8                                ###                              +0.0039
  -7                                #                                -0.0009
  -6                               ##                                -0.0014
  -5                                ##                               +0.0025
  -4                                #                                -0.0004
  -3                                ###                              +0.0042
  -2                               ##                                -0.0025
  -1                                ################################ +0.0778
   0                                #                                +0.0007
   1                                ##############                   +0.0316
   2                              ###                                -0.0051
   3                                #                                -0.0001
   4                                #                                -0.0003
   5                               ##                                -0.0018
   6                               ##                                -0.0021
   7                                #                                +0.0009
   8                                #                                -0.0012
   9                               ##                                -0.0015
  10                               ##                                -0.0014
  11                                ##                               +0.0017
  12                              ###                                -0.0043
  13                                #                                -0.0012
  14                               ##                                -0.0014
  15                             ####                                -0.0066
  16                               ##                                -0.0018
  17                                #                                -0.0012
  18                               ##                                -0.0027
  19                               ##                                -0.0023
  20                              ###                                -0.0045
```

Key findings:

- **Contemporaneous coupling is essentially absent.** At lag 0, r = +0.000668
  (n = 150,000); whether the loci are in contact right now says almost nothing
  about whether the gene is transcribed right now.
- **The strongest association is at lag -1**: r = +0.077843
  (bootstrap 95% CI [+0.069287, +0.086326]; naive p = 1.76e-199,
  n = 149,400 pooled pairs). Because this lag is
  negative, the dominant temporal ordering in the data is
  transcription at t-1 followed by contact at t: transcription tends to *precede* contact, while the reverse ordering (contact leading transcription) is present but weaker (lag +1, r = +0.031616).
- On the positive-lag side (contact before transcription), the strongest association is
  at lag +1 with r = +0.031616; on the negative-lag side the strongest is at
  lag -1 with r = +0.077843. The profile is therefore **asymmetric**: the
  negative-lag association is stronger than the positive-lag one.
- **Autocorrelation context:** neither signal is persistent one step apart (pooled lag-1 autocorrelation: contact +0.002, transcription +0.000; lag-2: contact -0.001, transcription -0.002). The lag +/-1 cross-associations are therefore not artifacts of slow, autocorrelated dynamics: contact at time t is specifically associated with transcription in the two adjacent steps (t-1 and t+1) but not with transcription at t itself.

Conditional rates at informative lags (pooled):

| lag | P(contact at t) | P(tx=1 at t+lag) overall | P(tx=1 at t+lag \| contact at t) | P(tx=1 at t+lag \| no contact at t) |
|----:|----------------:|-------------------------:|---------------------------------:|-------------------------------------:|
| -1 | 0.0534 | 0.1539 | 0.2722 | 0.1472 |
| +0 | 0.0534 | 0.1539 | 0.1549 | 0.1538 |
| +1 | 0.0534 | 0.1539 | 0.2020 | 0.1512 |

At the peak lag (-1), P(transcription = 1 at t-1 | contact at t) =
0.272 versus the overall rate 0.154 and the no-contact rate
0.147 - a modest enrichment consistent with the small but reliable
correlation above.

**Lag-selection uncertainty.** In the cell-level block bootstrap, the lag with the largest |r| was again lag -1 in 100.0% of replicates. The location and sign of the peak are robustly identified.

### 3.3 Sensitivity to the contact threshold

The supplied threshold is 260 nm. Re-running the lag scan with alternative
distance cutoffs:

| threshold (nm) | contact fraction | best lag (max abs r) | r at best lag | r at lag 0 |
|---------------:|-----------------:|---------------------:|--------------:|-----------:|
| 220 | 0.0302 | -1 | +0.058868 | +0.001605 |
| 240 | 0.0413 | -1 | +0.068790 | +0.002036 |
| 260 (supplied) | 0.0534 | -1 | +0.077843 | +0.000668 |
| 280 | 0.0653 | -1 | +0.083069 | +0.002815 |
| 300 | 0.0767 | -1 | +0.084320 | +0.000828 |
| 320 | 0.0884 | -1 | +0.083644 | +0.000694 |

The peak-lag pattern (strongest association at a small negative lag, near-zero
contemporaneous association) is stable across this threshold range; only the magnitude
of the association and the contact fraction change.

## 4. Temporal association is not causation

The lag analysis shows a **temporal association**: contact and transcription are not
independently distributed over time, and the association is asymmetric in lag, with the
maximum at lag -1. This is a statement about predictive temporal structure in
observational data. It is **not** evidence of a causal mechanism, for several reasons:

1. **No intervention.** Nothing in the data manipulates contact or transcription.
   Causal claims require comparing outcomes under interventions (e.g., forcing or
   preventing contact), which these observations never provide.
2. **Common-cause (confounding) explanations are not ruled out.** An upstream process -
   for example transcription-factor binding, local chromatin state remodeling, polymerase
   loading, cell-cycle phase, or nuclear microenvironment - could drive *both* contact and
   transcription with different response delays. Such a common cause can produce a
   lagged correlation whose peak sits at a negative, zero, or positive lag even if
   contact and transcription never interact directly.
3. **Reverse causation is not ruled out.** The peak at lag -1
   (transcription preceding contact) is equally compatible with
   the transcriptional machinery recruiting or stabilizing enhancer-promoter contact (transcription -> contact) as with any contact -> transcription story. Observational time ordering alone cannot
   choose between these directions.
4. **The peak lag is not a measured mechanistic delay.** Its location depends on the
   sampling interval, on the autocorrelation of both signals, and on any confounder
   dynamics; it cannot be read off as "transcription causes contact 1 step(s) later".
5. **Pooling can mask heterogeneity.** The pooled correlation averages over cells with
   different contact propensities and transcription rates; a within-cell latent state
   could contribute to the pooled association. (The between-cell correlation in 3.1 is
   near zero, so at least rate-level cell heterogeneity does not explain the signal.)

**What these observational data cannot establish:**

- that enhancer-promoter contact *causes* transcriptional activation (or repression);
- that transcription *causes* changes in contact;
- the *direction* of any causal influence between the two;
- that the association is *direct* rather than mediated by unmeasured molecular
  processes;
- that the observed lag (-1) reflects a true biological delay rather than
  confounder dynamics or sampling;
- the effect size of any hypothetical intervention on contact or transcription.

**What would be needed for causal evidence:** perturbations that exogenously set the
distance or contact state (optogenetic tethering of the enhancer to the promoter,
CRISPR deletion/inversion of the enhancer, acute degron-mediated depletion of bridging
factors), ideally with dose-response and rescue experiments, measured against matched
controls. Even then, causality would attach to the specific intervention, not to the
observational correlation reported here.

## 5. Limitations

- Contact is a binary summary of distance at a single supplied threshold (260 nm);
  the threshold-sensitivity table above shows the main pattern is robust, but magnitudes
  are threshold-dependent.
- The pooled Pearson correlation treats all within-cell pairs as one sample; inference
  therefore relies on the cell-level block bootstrap rather than naive p-values.
- 41 lags were examined; the peak lag was selected from the same data, so
  its exact value carries selection uncertainty (quantified by the bootstrap above).
- Time is in arbitrary observation steps; no physical time calibration is supplied, so
  no delay in physical units can be inferred.

## 6. Reproducibility

Run `python output/analysis.py` from the workspace root (or anywhere; paths are
anchored to the script location). It reads
`inputs/single_cell_dynamics_question.csv` and regenerates `output/cell_metrics.csv`,
`output/lag_analysis.csv`, and this `output/report.md`. Bootstrap seed: 20260828;
bootstrap replicates: 500.
