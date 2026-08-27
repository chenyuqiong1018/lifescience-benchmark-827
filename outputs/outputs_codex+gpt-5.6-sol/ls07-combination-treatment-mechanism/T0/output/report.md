# Combination-treatment mechanism report

The approved differential-expression analysis used only three `Cisplatin_IC50_CBD_IC50` samples (9-1, 9-2, 9-3) and three `DMSO` samples (3-1, 3-2, 3-3), with the combination as numerator. PyDESeq2 0.5.0 used `~Group`, `refit_cooks=True`, and the standard contrast. Genes were fitted when any selected raw count exceeded 10; 555 of 18,029 retained genes passed adjusted p-value < 0.05, absolute log2 fold change > 0.5, and base mean > 10.

Following the `go_term_analysis` functional-genomics workflow, passing genes were stratified into overall, upregulated, and downregulated sets before annotation. GSEApy 1.1.4 tested them only against the evaluator-supplied, hash-verified `Reactome_2022` GMT and its 10,489-symbol background. The external STRING, GO, and Ensembl endpoints described by the skill were not called because they would introduce non-frozen annotations forbidden by the task.

The best-supported primary mechanism is **TP53 Regulates Transcription Of Cell Cycle Genes** (R-HSA-6791312; overlap 8/49, adjusted p-value 0.0377, odds ratio 7.59). The related TP53 terms among the highest-ranked overall results provide coherent functional context. Selection was deterministic: minimum adjusted p-value, then raw p-value, then maximum combined score among all passing genes.

This Reactome over-representation is statistical association with the treatment contrast, not evidence that the combination caused the pathway change. The design also cannot separate synergy from either component's individual effect.
