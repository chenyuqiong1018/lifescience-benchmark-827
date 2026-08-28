# Enhancer-promoter distance and transcription dynamics

## Methods

The input contains 150,000 observations from 600 cells, with 250 time points per cell and no cross-cell joins. Enhancer-promoter distance is the three-dimensional Euclidean distance in nanometers. A contact is defined by the supplied threshold, distance <= 260 nm. Cell fractions use all 250 observations in each cell.

For every integer lag from -20 through +20, the analysis pools correctly aligned within-cell pairs and computes Pearson's correlation between contact at time *t* and transcription at time *t + lag*. Positive lag therefore means contact precedes transcription. Correlation itself is the requested effect-size measure; all 41 planned lags are reported without selecting only favorable results.

## Results

The largest absolute pooled temporal association occurs at lag -1: r = 0.077843, n = 149,400, r^2 = 0.006060. Under the specified sign convention, this means transcription leads later contact by 1 time step(s). The magnitude is weak: the variables share only about 0.606% of pooled variance at this selected lag. As a descriptive sensitivity check, 600 cells had nonconstant paired values at that lag; their median within-cell correlation was 0.062595, and 73.0% were positive.

The complete lag profile is in `lag_analysis.csv`; cell-level distance, contact, and transcription summaries are in `cell_metrics.csv`.

## Interpretation and limitations

This is an observational temporal association, not evidence that contact causes transcription (or that transcription causes contact). Selecting the strongest absolute value across 41 lags is exploratory, and adjacent time points are repeated, autocorrelated measurements within cells, so the pooled observation count is not a count of independent experimental units. The data alone cannot establish intervention effects, exclude common or time-varying confounding, determine molecular mechanism, or prove temporal direction at a resolution finer than the sampling interval. A causal claim would require an appropriate perturbation or other identification strategy, controls, and uncertainty analysis that respects cell-level clustering and serial dependence.
