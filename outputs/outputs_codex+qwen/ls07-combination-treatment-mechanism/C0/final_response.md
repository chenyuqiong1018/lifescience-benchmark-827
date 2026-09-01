Task complete. All five required artifacts are in `output/` and validated.

## What was done

**Pipeline** (`output/analysis.py`, run in a workspace venv):
1. **DE**: pydeseq2 0.5.0 (DESeq2, design `~condition`), `Cisplatin_IC50_CBD_IC50` (9-1/9-2/9-3) vs `DMSO` (3-1/3-2/3-3); Ensembl→HGNC mapping via the supplied TSV → 84 significant genes (padj < 0.05, |log2FC| ≥ 1; 59 up / 25 down)
2. **Enrichment**: GSEApy **1.1.4** `enrichr`, offline Fisher exact test against the byte-frozen `Reactome_2022.gmt` (1,818 terms) with the supplied `Reactome_2022.background.txt` universe (10,489 genes). Nothing downloaded or substituted (input SHA-256s re-verified in the manifest).

**Primary mechanism call — induction of AhR-associated xenobiotic and eicosanoid metabolism**: the combination strongly induces `AHRR` (log2FC +2.13, padj 9.8e-26), `CYP1A1` (+2.03, padj 1.1e-19), `CYP1B1` (+1.12), which drive the six top pathways (EET/DHET synthesis P=3.9e-04, HETE synthesis, Cytochrome P450, Xenobiotics, Endogenous Sterols, Phase-I functionalization). Secondary: CDKN1A/p53 cell-cycle-arrest component. Evidence strength is documented as moderate — min adjusted P = 0.070 (no term passes BH 0.05 with n=3 per group), so the call rests on convergent raw-P enrichment plus individually very significant markers.

**Artifacts** (all schema-asserted):
| file | check |
|---|---|
| `pathway_enrichment.csv` | 278 rows × 8 gseapy columns, sorted by Adjusted P-value, values in [0,1] |
| `mechanism_call.json` | parses; comparison, method, mechanism, markers, causation caveat |
| `resource_manifest.json` | parses; all 5 input SHA-256s match; gseapy 1.1.4 recorded |
| `analysis.py` | parses; staged (`full` / `finalize`) |
| `report.md` | 405 words ≤ 500; explicit enrichment-vs-causation section |

Plus supplementary `_de_results.csv` (full DE table). Environment notes: a prebuilt gseapy-1.1.4 wheel required Anaconda Python 3.13; anndata was pinned to 0.11.4 for pandas 2.2.3 compatibility.