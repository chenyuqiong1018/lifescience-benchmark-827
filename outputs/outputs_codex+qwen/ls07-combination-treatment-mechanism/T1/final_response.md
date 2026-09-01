All deliverables complete and validated. Summary:

## What was done

**Environment**: GSEApy 1.1.4 requires Python 3.12 (cp312 wheel); conda failed on the Unicode workspace path, so I bootstrapped an embedded Python 3.12.10 inside the workspace (`py312/`) with `gseapy==1.1.4` + `pydeseq2 0.5.4`. No pathway library or ID mapping was ever downloaded — only Python tooling.

**DE analysis** (`Cisplatin_IC50_CBD_IC50` 9-1/9-2/9-3 vs `DMSO` 3-1/3-2/3-3, PyDESeq2):
- 21,253 genes tested after count filter; 1,168 at padj<0.05; **82 DE genes** (56 up / 26 down) at padj<0.05 & |log2FC|≥1
- Leading hits: AHRR, CYP1A1, ALDH1A3, CDKN1A, TP53I3, FAS (up); DDIT3, NR1D1, FOS (down)

**Enrichment** (GSEApy 1.1.4 `enrichr`, local offline mode, frozen `Reactome_2022` GMT + supplied 10,489-gene background):
- Resource integrity verified (sha256 of GMT/background/mapping match the manifest; GMT gene-union == background)
- 275 of 1,818 terms tested (GSEApy natively tests only terms with ≥1 query-gene overlap; 38/82 DE genes are in the universe)
- **No term reaches FDR<0.05** (min adj-p = 0.062); top by raw p: EET/DHET synthesis (CYP1A1/CYP1B1, p=3.5e-4), HETE synthesis, Cytochrome P450, xenobiotics, PPARα lipid programs, then TP53 cell-death genes

**Mechanism call**: best-supported primary mechanism = **cytochrome P450-mediated xenobiotic/drug and lipid (arachidonic-acid/eicosanoid, PPARα) metabolism**, with weaker secondary TP53 cell-death support. The report explicitly separates enrichment (associative over-representation) from causation.

**Validation**: independent recomputation (scipy hypergeometric + BH, re-derived from the raw GMT/background/DE list) reproduced all 275 p-values, adjusted p-values, overlaps, odds ratios and combined scores exactly; final schema assertions on all artifacts passed.

## Artifacts in `output/`
| File | Content |
|---|---|
| `pathway_enrichment.csv` | 275 tested terms, columns Term/Overlap/P-value/Adjusted P-value/Odds Ratio/Combined Score/Genes |
| `mechanism_call.json` | primary mechanism, evidence counts, top-10 pathways, methods, causation caveat |
| `resource_manifest.json` | input manifest + run verification (sha256, counts, versions, `network_used: false`) |
| `analysis.py` | full deterministic offline pipeline |
| `report.md` | 396 words |
| `de_results.csv` (supplementary) | gene-level DE statistics |