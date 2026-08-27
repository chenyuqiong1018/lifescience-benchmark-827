# Frozen enhancer–promoter integration rule

1. For every Hi-C row, take the arithmetic mean of the three contact counts. On background rows only, fit ordinary least squares `log1p(mean_count) ~ log10(distance_bp)`.
2. Convert every residual to a robust z score using the background residual median and `1.4826 × MAD`. This is `contact_evidence`.
3. For each guide, calculate `log2((mean_perturbed + 0.5)/(mean_control + 0.5))`; the pair median is `perturbation_effect`.
4. A pair is physically eligible when `contact_evidence >= 1.645`. For eligible pairs, `combined_support = contact_evidence × abs(perturbation_effect)` and ranks run from strongest (1) to weakest. Ineligible pairs receive support and rank 0 and are not called causal candidates.
5. The least-supported causal pair is the eligible pair with the largest positive rank. Ties use `pair_id`.

The files contain EP1–EP7; there is no EP8 row to invent. This is a transparent local integration policy, not a benchmark-native hidden verifier.
