# Combination-treatment mechanism report

The approved six-sample differential-expression analysis compared `Cisplatin_IC50_CBD_IC50` (numerator; 9-1, 9-2, 9-3) with `DMSO` (denominator; 3-1, 3-2, 3-3), using only `Group` in the design. PyDESeq2 0.5.0 used `refit_cooks=True`; genes were retained when any selected raw count exceeded 10. Of 18,029 retained genes, 555 passed adjusted p-value < 0.05, absolute log2 fold change > 0.5, and base mean > 10.

GSEApy 1.1.4 tested the mapped passing genes against the evaluator-supplied, byte-verified `Reactome_2022` GMT and its explicit 10,489-gene background. No current pathway library or alternate mapping was downloaded. The all-passing query contained 282 unique background-mapped symbols; the up- and downregulated queries contained 187 and 95, respectively.

The best-supported primary cellular mechanism is **TP53 Regulates Transcription Of Cell Cycle Genes** (R-HSA-6791312; overlap 8/49, adjusted p-value 0.0377, odds ratio 7.59). This call follows a prespecified deterministic rule: among all-passing genes, minimize adjusted p-value and then raw p-value, with combined score as the next discriminator. The complete overall and directional enrichment tables are in `pathway_enrichment.csv`; resource provenance and integrity checks are in `resource_manifest.json`.

This is pathway over-representation associated with the treatment contrast, not proof that the combination treatment caused the pathway change. The experiment also does not isolate pharmacologic interaction from the effects of either component alone.
