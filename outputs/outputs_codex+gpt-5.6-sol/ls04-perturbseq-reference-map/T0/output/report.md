# Perturb-seq latent reference mapping

## Result

| Target gene | Query guide | 30D latent score | Runner-up score | Margin | Confidence |
|---|---:|---:|---:|---:|---:|
| PABPC1 | guide18 | 0.843426 | 0.445620 | 0.397806 | 1.000 |
| NUDT21 | guide13 | 0.766347 | 0.610190 | 0.156157 | 1.000 |
| LEO1 | guide14 | 0.582859 | 0.464379 | 0.118480 | 1.000 |

The requested assignments are **PABPC1 → guide18**, **NUDT21 → guide13**, and **LEO1 → guide14**.

## Skill-guided method

The selected skill was `scvi-tools`, chosen for its batch-corrected latent-space and reference label-transfer workflow. The installed guidance requires raw integer UMI counts for an scVI negative-binomial model. The execution environment had `torch` but did not have `scvi` or `anndata`, so no scVI model was claimed or trained.

The auditable CPU fallback preserves and validates the raw counts, then uses log-CP10K only to construct a low-dimensional surrogate representation. Query and reference guide means are independently centered on their own pooled non-targeting cells to remove the dominant dataset/cell-type shift. Gene IDs are intersected and explicitly reordered (11,858 shared genes), and 2,000 features are selected by pooled variance across anonymous guide signatures. An exact Gram-matrix eigendecomposition gives latent guide coordinates; cosine similarity in the preselected 30-dimensional representation supplies the reported scores.

## Leakage and ambiguity controls

The complete score matrices for latent dimensions 10, 20, 30, 40, 50, and 60 are computed before `obs/target_gene` is opened. Target labels therefore only choose report rows and cannot influence preprocessing, feature selection, latent construction, or similarity scores.

`runner_up_score` is the second-highest query-guide score at 30 dimensions. `confidence` is the fraction of the six fixed latent dimensionalities that retain the reported top guide; all three are 1.000. The NUDT21 and LEO1 margins are narrower than the PABPC1 margin, so those mappings remain relatively more ambiguous despite their dimensional stability. This confidence is a sensitivity measure, not a calibrated biological probability.
