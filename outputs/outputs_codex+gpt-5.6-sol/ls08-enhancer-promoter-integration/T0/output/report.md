# Enhancer-promoter evidence integration

Following the `regulatory_region_analysis` evidence-separation perspective, all supplied rows were joined only by `pair_id`. Physical Hi-C contact and CRISPR-expression response remained distinct. No Ensembl/UCSC overlap, sequence, binding-matrix, or phenotype endpoint was called: the frozen pair-level files provide no genomic coordinates and the supplied integration rule does not permit external annotation.

For Hi-C, the mean of three counts was transformed with `log1p`; ordinary least squares against `log10(distance_bp)` used 200 background rows. Residuals were standardized by the background median and 1.4826×MAD. For each guide, expression effect was `log2((mean_perturbed + 0.5)/(mean_control + 0.5))`; pair effects are medians across 4 guides.

Pairs required contact evidence ≥ 1.645. Eligible pairs were ranked from strongest to weakest by contact evidence times absolute perturbation effect; ineligible pairs received support and rank zero. The least-supported eligible pair is **EP3** (contact 3.1048, perturbation -0.0108, support 0.0335, rank 4). Inputs contain EP1-EP7 only; EP8 was not invented.

The term “causal candidate” follows the frozen local rule. Contact plus expression response is stronger than either modality alone, but it does not itself prove a direct causal enhancer-promoter relationship or rule out indirect perturbation effects.
