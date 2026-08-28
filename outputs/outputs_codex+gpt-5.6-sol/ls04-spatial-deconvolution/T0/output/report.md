# Skill-guided spatial composition of Spot_710-1

## Result

| Cell type | Reported weight | Marker-model MLE | Joint bootstrap 95% CI | Marker-count sensitivity |
|---|---:|---:|---:|---:|
| B_Cell | 0.329049 | 0.327 | [0.294, 0.358] | [0.320, 0.334] |
| Endothelial | 0.360644 | 0.359 | [0.328, 0.389] | [0.350, 0.370] |
| Macrophage | 0.310308 | 0.309 | [0.275, 0.333] | [0.301, 0.320] |

The nonnegative reported weights sum to 1.000001. `Spot_710-1` is a three-way mixture, not a single-type spot.

## Skill-guided method

The autonomously selected whitelist skill was `scvi-tools`. Its spatial-task guidance explicitly routes deconvolution to raw-count reference-signature models such as cell2location, DestVI, or Tangram. Those model runtimes and their GPU environment were unavailable, so the analysis uses a transparent cell2location-inspired multinomial abundance model and does not claim to have trained cell2location or scVI.

Raw nonnegative integer counts were validated. Each reference cell was library-normalized, and the mean of 200 cells defined each of six cell-type signatures. For the primary model, the 100 most specific genes per type were selected, producing 598 unique markers. The target marker counts were fit as a convex mixture of the marker signatures by multinomial maximum likelihood.

## Evidence and ambiguity

The mixture cosine is 0.925, compared with 0.707 for the best pure type on the marker panel. A 200-replicate joint bootstrap resampled both the target spot counts and reference cells within each type, reselected markers, and refit the mixture. All three main-type intervals remain well above zero.

Repeating the fit with 50, 75, 100, 150, and 200 markers per type retained the same three components within the reported sensitivity ranges. The unthresholded model assigned approximately 0.006 to Fibroblast_Stroma and zero to T_Cell and Tumor_Core; components below 0.02 were treated as unstable tails and the retained components were renormalized. The relative ordering is somewhat uncertain because the three weights are close, but their joint presence is robust.
