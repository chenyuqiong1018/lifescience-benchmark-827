# Spot_710-1 spatial composition

## Result

`Spot_710-1` is a supported three-way mixture:

| Cell type | Reported weight | Unthresholded MLE | Parametric bootstrap 95% CI |
|---|---:|---:|---:|
| B_Cell | 0.332577 | 0.331 | [0.306, 0.359] |
| Endothelial | 0.352118 | 0.351 | [0.327, 0.371] |
| Macrophage | 0.315306 | 0.314 | [0.284, 0.336] |

The reported nonnegative weights sum to 1.000001. The spot should not be represented as a single cell type.

## Method

The single-cell count matrix contains 1,200 cells, with 200 cells in each of six reference types. Each cell was converted to a within-cell expression proportion; the mean of those profiles defined each type signature. Gene names and ordering were explicitly checked against all 1,000 Visium features.

The target spot has 3,693 counts across 883 detected genes. Its count vector was fit as a convex combination of the six type signatures using multinomial maximum likelihood (EM). A fixed-seed, 200-replicate parametric bootstrap resampled 3,693 counts from the fitted profile and refit the mixture to quantify sampling sensitivity.

## Mixture evidence and ambiguity

The fitted mixture has cosine similarity 0.915 to the observed expression profile, compared with 0.723 for the best pure reference type. This large fit improvement and the three bootstrap intervals bounded away from zero support a genuine mixture.

The unthresholded fit assigned 0.0039 to Fibroblast_Stroma and effectively zero to T_Cell and Tumor_Core. Components below the prespecified 0.02 reporting threshold were treated as unstable tails; the three retained weights were renormalized. The main components have similar weights, with Endothelial slightly largest, so their precise ordering is less certain than the conclusion that all three are present.
