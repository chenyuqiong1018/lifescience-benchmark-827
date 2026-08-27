# Enhancer-promoter evidence integration

All supplied rows were integrated by `pair_id`, preserving physical contact and CRISPR-expression evidence as distinct measurements. The Hi-C model used 200 background rows: ordinary least squares of `log1p(mean_count)` on `log10(distance_bp)`. Residuals were standardized with the background residual median and 1.4826 times MAD to produce `contact_evidence`.

For each guide, the perturbation effect was `log2((mean_perturbed + 0.5)/(mean_control + 0.5))`; each pair's `perturbation_effect` is the median across 4 guides. Pairs with contact evidence below 1.645 were ineligible and assigned support and rank zero. Eligible pairs were ranked from strongest to weakest by `contact_evidence × abs(perturbation_effect)`, with `pair_id` as the deterministic tie-breaker.

The least-supported eligible causal candidate under this frozen rule is **EP3** (contact evidence 3.1048, perturbation effect -0.0108, combined support 0.0335, rank 4). The files contain exactly EP1-EP7; no EP8 row was invented.

“Causal candidate” is a label defined by the supplied integration rule. Physical proximity and expression response are complementary evidence, but this calculation alone does not establish a biological causal enhancer-promoter link or exclude indirect perturbation effects.
