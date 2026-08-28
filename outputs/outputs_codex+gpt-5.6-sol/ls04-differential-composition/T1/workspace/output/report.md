# Retinal composition comparison

## QC and frozen analysis

Sample 1 contains 36,601 genes x 6,295 cells and 9,246,637 stored count entries; sample 2 contains 36,601 genes x 5,004 cells and 9,953,348 entries. No cell or gene was removed. All libraries were positive, header entry totals matched the parsed entries, and every listed marker symbol occurred exactly once.

Each cell's raw integer counts were divided by its full library size, multiplied by 10,000, and transformed with `log1p`. For each panel row, the arithmetic mean of its listed markers was computed; the largest score assigned the type, with panel order resolving ties. Fractions use the complete matrix column count. The call was restricted to types with sample-1 fraction >= 1% and selected by the minimum sample-2/sample-1 fraction ratio.

Library-size QC: sample 1 range 205–45,431, median 1999.0; sample 2 range 204–48,129, median 3197.0.

## Depleted population

The depleted call is **horizontal cell**. It decreases from 237/6,295 cells (fraction 0.037649) to 49/5,004 (fraction 0.009792); the sample-2/sample-1 fraction ratio is 0.260092. Thus, under the frozen labels, the sample-2 proportion is about 26.0% of the sample-1 proportion.

Conditional on treating cells as independent observations, descriptive 95% Wilson intervals are [0.033222, 0.042640] and [0.007415, 0.012921], and the log-ratio approximation is [0.191656, 0.352964]. These intervals are shown in `composition_profile.svg`; they do not represent donor-level uncertainty.

## Annotation evidence, model scope, and uncertainty

Annotation evidence is restricted to the supplied marker sets. Sample 1 has 172 exact winning-score ties (median winning margin 0.993560); sample 2 has 26 (median margin 2.176783). Ties were retained and resolved as required.

The matrices are raw integer counts, which would be appropriate input to probabilistic single-cell models. However, scVI/scANVI training, scGPT embeddings, vocabulary-based gene dropping, learned clustering, or label transfer would change the explicitly frozen normalization and annotation rule. They were therefore not substituted for the required calculation. No checkpoint, reference labels, donor/batch covariates, or GPU modeling is needed to reproduce this artifact.

Biological uncertainty remains: the rule does not handle doublets, ambient RNA, marker overlap, continuous states, batch effects, or donor replication. With one aggregate matrix per sample, cells may be pseudoreplicates; the analysis cannot establish population-level variance, a causal loss mechanism, or generalization beyond these matrices. The result is a deterministic composition call under the supplied rule.
