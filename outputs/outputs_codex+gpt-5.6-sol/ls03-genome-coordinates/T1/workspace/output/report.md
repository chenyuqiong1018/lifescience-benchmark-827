# Enhancer-promoter contact and transcription dynamics

## Data and analysis

The dataset has 150,000 complete observations: 600 cells with 250 ordered time points each. For every observation, enhancer-promoter distance is the Euclidean norm of the x/y/z coordinate difference in nanometers. Contact is 1 when distance is <= 260 nm, inclusive, and 0 otherwise. Per-cell fractions use all 250 rows.

For each integer lag from -20 through +20, contact(t) is paired only with transcription(t + lag) from the same cell. The requested pooled Pearson correlation and exact number of aligned pairs are reported for every lag. Positive lag means contact leads later transcription. No cross-cell pairs are formed.

The overall mean distance is 526.067 nm, the pooled contact fraction is 0.0534, and the pooled transcription fraction is 0.1539. These pooled summaries do not replace the 600 per-cell values in `cell_metrics.csv`.

## Lag result

The largest absolute association is at lag -1: r = 0.077843, n = 149,400, r^2 = 0.006060. By the supplied sign convention, transcription precedes contact by 1 time step(s). This is a weak association in magnitude. A 2,000-resample whole-cell bootstrap, conditional on examining this selected lag, gives a descriptive 95% interval of [0.069653, 0.086153]. It preserves cell series during resampling, but it does not correct for choosing the maximum absolute value among 41 lags.

The complete numerical profile is in `lag_analysis.csv`. `lag_profile.svg` displays all 41 prespecified lags with a zero reference line, accessible text, and colorblind-safe styling; the plotted point count for each lag is available in the CSV rather than implied to be constant.

## Regulatory annotation scope

The available measurements are microscopy-style Cartesian positions in nanometers. They are not chromosome coordinates and provide no organism, chromosome, assembly, interval, or gene identifier. Consequently, genomic-region overlap, regulatory-element-to-gene, sequence, binding-matrix, or phenotype-region queries are not identifiable from this input and were not fabricated or run. Such annotation would require a real genomic interval plus organism/assembly metadata.

## Association is not causation

The lag pattern is observational. It cannot establish that enhancer-promoter contact causes transcription, that transcription causes subsequent contact, or that either direction is free of shared or time-varying causes. Repeated measurements are serially dependent, and the 149,400 pooled pairs at lag -1 are not 149,400 independent experimental units. The data cannot establish a molecular mechanism, rule out confounding, resolve direction below the sampling interval, or estimate an intervention effect. Causal conclusions require a suitable perturbation or identification design, appropriate controls, and inference that accounts for cells and temporal dependence.
