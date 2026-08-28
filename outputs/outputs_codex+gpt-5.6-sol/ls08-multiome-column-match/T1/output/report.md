# Multiome column matching

## Result

The recovered bijection is RNA 0 → ATAC 5, RNA 1 → ATAC 1, RNA 2 → ATAC 4, RNA 3 → ATAC 0, RNA 4 → ATAC 6, RNA 5 → ATAC 3, RNA 6 → ATAC 7, RNA 7 → ATAC 2. The maximum total Pearson correlation is 3.663500. The next-best complete assignment totals 3.629567, a global gap of 0.033933. The smallest assigned-minus-row-runner-up margin is 0.012746; a negative row margin can occur because the optimum is global rather than eight independent row choices.

## Shared biological signal and method

The shared signal is coordinated gene activity across molecular layers: accessible chromatin at a gene's transcription start site tends to accompany expression of that gene. Following the supplied rule, Ensembl release 112 symbols were retained only when they occurred exactly once in the complete annotation, appeared in RNA, and lay on chromosomes 1–22, X, or Y. Strand determined whether gene start or end was the TSS. One-based coordinates were converted to zero-based half-open 10 kb ATAC bins. RNA TPM and mapped ATAC values were log1p-transformed; the 2,000 mapped genes with greatest across-RNA-population variance were used for all 64 Pearson correlations. An O(n³) Hungarian maximum-weight solution enforced the bijection, and exhaustive enumeration of all 8! assignments independently confirmed its total.

## Data-quality audit

The inputs contained 29,165 unique RNA genes, 295,420 unique ATAC bins, and 70,611 annotation rows. After exact-symbol/chromosome filtering there were 18,921 eligible genes; 18,875 had an observed TSS bin. Every matrix row had nine fields, all eight measurement columns were unique, and no blank, negative, or non-finite measurement was found. No zero imputation was performed. `runner_up_score` is the second-highest row-wise correlation, not the second-best global assignment. Correlations describe association, not gene-level causation.

The five prescribed skills guided cross-omics framing, regulatory-coordinate handling, integrity checks, correlation interpretation, and executable-code verification; external reference datasets were neither needed nor used.
