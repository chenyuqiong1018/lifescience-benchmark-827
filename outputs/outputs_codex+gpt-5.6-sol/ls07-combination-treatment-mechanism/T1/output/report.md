# Combination-treatment mechanism report

The approved differential-expression analysis compared three `Cisplatin_IC50_CBD_IC50` samples (9-1, 9-2, 9-3) with three `DMSO` samples (3-1, 3-2, 3-3), combination as numerator. PyDESeq2 0.5.0 used `~Group`, `refit_cooks=True`, and the standard contrast. Genes were fitted when any selected raw count exceeded 10; 555 of 18,029 retained genes passed adjusted p-value < 0.05, absolute log2 fold change > 0.5, and base mean > 10.

GSEApy 1.1.4 tested overall, upregulated, and downregulated passing-gene sets only against the evaluator-supplied, hash-verified `Reactome_2022` GMT and its 10,489-symbol background. No current library, alternate mapping, STRING network, GO annotation, or Ensembl annotation was fetched. The code-execution endpoint only echoed a logic check, so the local executable analysis is authoritative.

The best-supported primary mechanism is **TP53 Regulates Transcription Of Cell Cycle Genes** (R-HSA-6791312; overlap 8/49, adjusted p-value 0.0377, odds ratio 7.59). Related TP53 transcription terms among the highest-ranked overall pathways support a coherent TP53/cell-cycle interpretation. Selection was prespecified and deterministic: adjusted p-value, raw p-value, then combined score.

Critical appraisal limits confidence. There are only three replicates per group, the design includes only `Group`, and no single-agent arms enter this contrast; therefore synergy and component attribution cannot be inferred. Over-representation also ignores gene correlation and pathway redundancy. The result is an experiment-specific statistical association, not proof that the combination caused TP53 pathway activation or the observed expression changes.
