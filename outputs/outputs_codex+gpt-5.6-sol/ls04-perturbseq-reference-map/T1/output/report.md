# Multi-view perturb-seq reference mapping

## Result

| Target gene | Query guide | Ensemble score | Runner-up score | Margin | Confidence |
|---|---:|---:|---:|---:|---:|
| PABPC1 | guide18 | 0.665450 | 0.346476 | 0.318974 | 1.000 |
| NUDT21 | guide13 | 0.697352 | 0.573815 | 0.123537 | 1.000 |
| LEO1 | guide14 | 0.439745 | 0.322117 | 0.117628 | 1.000 |

The requested mappings are **PABPC1 → guide18**, **NUDT21 → guide13**, and **LEO1 → guide14**.

## Controlled-skill analysis

The controlled skills were `scvi-tools`, `scgpt`, `statistical-analysis`, and `code_execution_analysis`.

- `scvi-tools` motivated raw-integer-count validation, within-dataset shift removal, and a 30-dimensional latent transfer view. The runtime lacked `scvi`, `anndata`, CUDA, and a GPU, so the analysis uses an explicit Gram-SVD latent proxy and does not claim scVI/scANVI training.
- `scgpt` motivated a complementary gene-order representation. The runtime lacked `scgpt` and a checkpoint, so the analysis uses a transparent per-signature gene-rank correlation proxy and does not claim foundation-model inference.
- `statistical-analysis` motivated a prespecified sensitivity analysis: 200 fixed-seed feature bootstraps, complete reporting of runner-up scores, and separation of stability from absolute probability. No cell-level p-value is manufactured because cells assigned to one guide are not independent perturbation replicates for this matching question.
- `code_execution_analysis` was used to formulate a deterministic score-aggregation unit test. Its returned code was executed locally and passed; the full local script is the auditable execution record. Dataset and literature search tools were intentionally not used.

## Preprocessing and score

Raw count matrices are validated before any transformation. Cells are library-size normalized to 10,000 and log-transformed only for the documented local proxies. Guide means are centered on the pooled non-targeting-cell mean separately in query and reference, reducing the cell-type/dataset shift. Gene IDs are explicitly intersected and reordered, yielding 11,858 shared genes, and the 2,000 highest-variance features are selected without labels.

The reported score is the equal-weight mean of three cosine-scale views: direct perturbation-signature cosine, 30D latent cosine, and centered gene-rank correlation. `runner_up_score` is the next-largest ensemble score among query guides. `confidence` is the fraction of 200 feature bootstraps retaining the reported top guide. All three are 1.000, but the smaller NUDT21 and LEO1 margins show greater relative ambiguity than PABPC1. Bootstrap support is a feature-sensitivity measure, not a calibrated probability of biological correctness.

## Leakage firewall

The full direct, rank, latent, and ensemble score matrices—and bootstrap winners for every reference row—are completed before `obs/target_gene` is opened. Reference target labels only select the three rows to report; they cannot affect normalization, feature selection, representation learning, score construction, or bootstrap outcomes. Query target metadata are absent.
