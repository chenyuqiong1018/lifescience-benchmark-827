# Enhancer-promoter distance and transcription dynamics

## Summary

The dataset contains 600 cells and 150,000
observations, with exactly 250 timepoints per cell. Across cells,
the mean of cell-level mean enhancer-promoter distances is 526.067 nm,
the mean contact fraction at the supplied `≤260 nm` threshold is
0.0534, and the mean transcription fraction is 0.1539.

The strongest absolute pooled lag association occurs at **lag -1**:
Pearson `r = 0.077843` from
149,400 within-cell pairs. Under the supplied sign convention,
transcription leads the later contact measurement. At lag 0, `r = 0.000668` from
150,000 pairs.

## Method

For every row, enhancer-promoter distance was calculated as the three-dimensional
Euclidean distance in nanometres. Contact equals 1 when distance is at most 260 nm.
Cell summaries use all 250 rows for that cell. For each integer lag from -20 through
+20, contact at time `t` was paired with transcription at `t+lag` only when both
indices existed within the same cell. Pearson correlation was then calculated after
pooling those valid within-cell pairs. No pairs cross cell boundaries.

## Interpretation and limitations

The lag peak is a temporal association, not evidence that physical contact causes
transcription. These observational trajectories cannot establish causal direction,
exclude a shared upstream regulator, remove time-dependent or cell-state confounding,
or show that experimentally changing distance would change transcription. Binary
contact also discards variation above and below the 260-nm threshold. A causal claim
would require an intervention (for example targeted perturbation of contact), suitable
controls, and a design that addresses temporal autocorrelation and cell-level effects.

## Reproduction

Run `python analysis.py` from any directory. The script locates the repository from
its own path, validates the 250-point within-cell time grids, and deterministically
rewrites both CSV files and this report.
