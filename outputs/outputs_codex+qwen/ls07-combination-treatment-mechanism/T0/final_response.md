Task complete. All five required artifacts are in `output/`, validated against exact schemas.

**Pipeline (single production run)**
- **Environment**: gseapy has no wheel for the system Python 3.14, so I created a workspace venv on Python 3.13.5 with the prebuilt **gseapy==1.1.4** wheel + pydeseq2 0.5.0, plus a documented anndata-0.13 compatibility shim (keeps 1-D varm masks; no effect on statistics).
- **DE**: PyDESeq2 0.5.0, `~condition`, contrast `Cisplatin_IC50_CBD_IC50` (9-1..9-3) vs `DMSO` (3-1..3-3); ENSG→symbol via supplied mapping (duplicates summed), 50,746 genes tested → **84 DEGs** at padj<0.05 & |log2FC|≥1 (59 up / 25 down).
- **Enrichment**: GSEApy 1.1.4 `enrichr` local mode against the frozen `Reactome_2022.gmt` (1,818 terms; sha256 verified vs supplied manifest) with the supplied 10,489-gene background universe. One-sided hypergeometric + BH. 278 terms had ≥1 overlapping gene and are reported; **no term passed BH<0.05** (best adjusted p = 0.070).
- **Mechanism call**: best-supported primary mechanism = **cytochrome P450–mediated xenobiotic/drug metabolism**, consistent with AhR-pathway activation (CYP1A1, CYP1B1, AHRR are top up-regulated genes and drive the top terms: EET/DHET & HETE synthesis, CYP450, Xenobiotics, Phase I functionalization, PPARα lipid metabolism).

**Independent validation**: all 278 hypergeometric p-values and BH-adjusted values recomputed from the frozen GMT/universe with SciPy — exact match with `pathway_enrichment.csv`; GMT/background/mapping sha256 all match the input manifest.

**Artifacts** (`output/`): `pathway_enrichment.csv` (278 rows × 8 gseapy columns), `mechanism_call.json`, `resource_manifest.json`, `analysis.py`, `report.md` (371 words, explicitly distinguishes enrichment association from demonstrated causation), plus supplementary `deseq2_results.csv`.