# Enhancer-promoter evidence integration

All rows were joined only by `pair_id`; Hi-C contact and CRISPR-expression response remained distinct. Mean Hi-C counts were modeled on 200 background rows by ordinary least squares, `log1p(mean_count) ~ log10(distance_bp)`. Residuals were standardized with the background median and 1.4826×MAD. Each guide's expression effect was `log2((mean_perturbed + 0.5)/(mean_control + 0.5))`, and each pair used the median across 4 guides.

Pairs required contact evidence ≥ 1.645. Eligible pairs were ranked by contact evidence times absolute perturbation effect, strongest first; `pair_id` broke ties. Ineligible pairs received support and rank zero. The least-supported eligible pair is **EP3** (contact 3.1048, perturbation -0.0108, support 0.0335, rank 4). EP8 was not invented.

The regulatory-region and IGVF skills require genomic coordinates, gene identities, or biosample context absent from the frozen inputs, so Ensembl, UCSC, and IGVF were not queried. The code endpoint only echoed a logic check; the local script is authoritative.

Critical appraisal: the two modalities and replicated guides strengthen the frozen-rule ranking, but Hi-C proximity is not direct regulation, and CRISPR effects may be indirect or off-target. The product score has no uncertainty interval. Therefore EP3 is the least-supported **rule-defined candidate**, not a demonstrated biological causal link. The seven supplied pairs also limit conclusions about the originally mentioned eighth candidate.
