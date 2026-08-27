# Perturb-seq reference mapping

## Result

| Target gene | Query guide | Score | Runner-up score | Margin | Confidence |
|---|---:|---:|---:|---:|---:|
| PABPC1 | guide18 | 0.628924 | 0.394111 | 0.234813 | 1.000 |
| NUDT21 | guide13 | 0.701070 | 0.593042 | 0.108028 | 1.000 |
| LEO1 | guide14 | 0.438264 | 0.339762 | 0.098502 | 1.000 |

The requested mappings are therefore **PABPC1 → guide18**, **NUDT21 → guide13**, and **LEO1 → guide14**.

## Method

The query and reference matrices were independently library-size normalized to 10,000 counts per cell and transformed with `log1p`. For each guide group, its mean profile was centered against the pooled non-targeting-cell mean from the same dataset. This within-dataset contrast removes much of the cell-type shift before cross-dataset comparison.

Gene IDs were explicitly intersected and reordered, yielding 11,858 shared genes; the raw column orders were not assumed to match. The 2,000 genes with the greatest variance across the combined anonymous guide-level signatures were selected without labels. Reference-to-query scores are cosine similarities on these signatures.

## Leakage and ambiguity controls

`obs/target_gene` is not read until after the complete 33-by-32 similarity matrix and the unsupervised feature set have been computed. Thus target labels can select the three report rows but cannot affect normalization, feature selection, or scoring. Query data contain only anonymous guide IDs.

`runner_up_score` is the second-highest query-guide similarity for the same reference group. `confidence` is the fraction of 500 fixed-seed feature bootstraps for which the reported guide remains top-ranked. All three mappings are feature-stable (1.000), but NUDT21 and LEO1 have smaller raw-score margins than PABPC1 and should be regarded as relatively more ambiguous. Bootstrap support measures stability to feature resampling, not an absolute probability that the biological assignment is correct.
