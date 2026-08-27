# Frozen distance/lag rule

`cell_metrics.csv` must contain `cell_id,n_timepoints,mean_distance_nm,contact_fraction,transcription_fraction`. Euclidean distance uses the supplied x/y/z coordinates in nm; contact is distance `<=260`; fractions use all 250 rows per cell.

`lag_analysis.csv` must contain one row for every integer lag from -20 through +20 with `lag,association,n_observations`. The association is the pooled Pearson correlation across all cells between contact at time `t` and transcription at time `t+lag`; positive lag means contact leads the later transcription value. Do not join observations across cell boundaries. The strongest absolute association is interpreted as temporal association only, never proof of causation.
