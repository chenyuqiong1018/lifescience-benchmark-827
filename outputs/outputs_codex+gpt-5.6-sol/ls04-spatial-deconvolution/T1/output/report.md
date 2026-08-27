# Multi-view composition of Spot_710-1

## Result

| Cell type | Ensemble weight | Joint bootstrap 95% CI | All-gene / marker / rank views |
|---|---:|---:|---:|
| B_Cell | 0.335656 | [0.309, 0.366] | 0.331 / 0.327 / 0.345 |
| Endothelial | 0.350943 | [0.326, 0.379] | 0.351 / 0.359 / 0.340 |
| Macrophage | 0.313401 | [0.285, 0.343] | 0.314 / 0.309 / 0.315 |

The nonnegative weights sum to exactly 1.000000. `Spot_710-1` is a supported three-way mixture rather than a pure cell type.

## Controlled-skill analysis

The controlled skills were `scvi-tools`, `scgpt`, `statistical-analysis`, `code_execution_analysis`, and `scientific-visualization`.

- `scvi-tools` routed the spatial task to raw-count reference-signature deconvolution rather than ordinary scVI latent clustering. `scvi`, `anndata`, CUDA, cell2location, and DestVI were unavailable, so the first two views are transparent all-gene and marker-count multinomial proxies; no unavailable model is claimed.
- `scgpt` motivated a scale-robust gene-order view. Both the package/checkpoint and vocabulary-compatible gene symbols were unavailable—the simulation uses synthetic `Gene_N` names—so the third view is an explicit rank-based NNLS proxy rather than claimed scGPT inference.
- `statistical-analysis` motivated the prespecified equal-weight three-view ensemble, joint reference/spot bootstrap, complete intervals, and effect-size reporting. No cell-level p-value was manufactured because the task is composition estimation, not a replicate-level group hypothesis test.
- `code_execution_analysis` supplied an EM unit-test code template. It returned code rather than an executed result; a mismatched toy count vector was caught before local execution, corrected to an exact 25/75 mixture, and then recovered `[0.25, 0.75]` with sum 1.
- `scientific-visualization` guided the diagnostic figure: a zero-baseline bar chart with 95% CI error bars, Okabe–Ito colors, redundant hatching for grayscale accessibility, sans-serif typography, lossless 300-DPI PNG, and vector SVG. The rendered PNG was visually inspected and passed layout/readability review.

## Evidence and uncertainty

Raw nonnegative integer counts, matching cell/metadata indices, and identical spatial/reference gene ordering were validated. Each view used the same six reference types and 1,000 genes; the marker view selected 100 type-specific genes per reference type. The reported estimate is the equal-weight mean of the all-gene count, marker-count, and rank-proxy weights, followed by a 0.02 tail threshold and renormalization.

The mixture has cosine similarity 0.915 to the observed spot profile, compared with 0.723 for the best pure type. Its count log-likelihood exceeds the best pure signature by 2082.816 log units. Across 200 fixed-seed joint bootstraps that resampled target counts and the 200 reference cells within each type, all three 95% intervals remain well above zero. The three view estimates are also close, supporting presence of all three types; their relative ordering is less certain because the proportions are similar.

The diagnostic is available as `composition_diagnostic.png` and `composition_diagnostic.svg`; error bars are joint-bootstrap 95% intervals (`n=200`).
