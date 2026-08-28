# Multiome column matching

## Result

The maximum-total-correlation bijection is: RNA 0 → ATAC 5, RNA 1 → ATAC 1, RNA 2 → ATAC 4, RNA 3 → ATAC 0, RNA 4 → ATAC 6, RNA 5 → ATAC 3, RNA 6 → ATAC 7, RNA 7 → ATAC 2. The assigned correlations sum to 3.663500; the smallest assigned-minus-row-runner-up margin is 0.012746.

## Shared biological signal

RNA abundance and chromatin accessibility are different molecular layers, but active regulatory programs tend to make transcription start sites accessible while their linked genes are expressed. I therefore compared populations through matched gene-level profiles: each eligible gene's strand-aware TSS was mapped to its containing 10 kb ATAC bin, and the corresponding RNA TPM and ATAC value were log-transformed. This is the cross-modal biological signal used for matching; it does not imply a causal relationship for individual genes.

## Procedure and safeguards

The Ensembl release 112 annotation was filtered to symbols occurring exactly once in the complete annotation, present in RNA, and located on chromosomes 1–22, X, or Y. Coordinates were treated as one-based and converted to zero-based half-open ATAC bins. No missing value was replaced with zero. Among genes with an observed target bin, the 2,000 with highest variance across log1p RNA populations were selected. All 64 Pearson correlations were computed across the same ordered genes, then an O(n³) Hungarian maximum-weight assignment enforced the required bijection. `runner_up_score` is the second-highest correlation in that RNA row, irrespective of the global assignment.

The freshly installed `multiomics_integration` skill informed the cross-layer integration framing; its external TCGA, UniProt, STRING, and KEGG calls were not relevant to this self-contained column-matching input and were not used.
