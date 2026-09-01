# Combination-treatment mechanism: differential expression and Reactome enrichment

## Setup
Comparison: `Cisplatin_IC50_CBD_IC50` (n=3) vs `DMSO` vehicle (n=3) from the frozen bix-43 raw counts. Differential expression (DE): PyDESeq2 0.5.4 (DESeq2-style negative-binomial GLM), genes with total counts >= 10 across the 6 samples (21,253 tested); DE defined as adjusted p < 0.05 and |log2FC| >= 1. Enrichment: GSEApy 1.1.4 `enrichr` in local offline mode against the evaluator-frozen `Reactome_2022` GMT (1,818 terms) with the supplied background universe (10,489 genes); one-sided hypergeometric (Fisher) test, Enrichr odds ratio and combined score, Benjamini-Hochberg FDR. No pathway library or identifier mapping was downloaded or substituted; GMT, background and mapping sha256 values were verified against the input manifest.

## Differential expression
1,168 genes had adjusted p < 0.05; 82 genes met the full DE criteria (56 up, 26 down). Top up-regulated genes: GDA, AHRR, ALDH1A3, CYP1A1, CDKN1A, ZBED2, FAS, TP53I3; top down-regulated: ANKRD37, DDIT3, NR1D1, FOS. The leading hits combine an aryl-hydrocarbon/xenobiotic-response module (AHRR, CYP1A1) with p53-mediated stress signaling (CDKN1A, TP53I3, FAS). Full statistics: `de_results.csv`.

## Pathway enrichment
38 of the 82 DE symbols fall inside the background universe; 275 of 1,818 Reactome terms share at least one gene with the DE list and are tested (GSEApy's native behaviour). No term reaches FDR < 0.05 (minimum adjusted p = 0.062), so ranking rests on raw p-values:

1. Synthesis of epoxy (EET) and dihydroxyeicosatrienoic acids (2/8, p = 3.5e-4; CYP1A1, CYP1B1)
2. Synthesis of (16-20)-hydroxyeicosatetraenoic acids (2/9, p = 4.5e-4)
3. Cytochrome P450 arranged by substrate type (3/65, p = 1.6e-3)
4. Xenobiotics (2/24); endogenous sterols (2/27); phase I functionalization (3/104)
5. PPARA activates gene expression; regulation of lipid metabolism by PPARalpha
6. TP53 regulates transcription of cell death genes (2/44, p = 0.011)

## Mechanism call
The best-supported primary cellular mechanism is **cytochrome P450-mediated xenobiotic/drug and lipid metabolism**, specifically arachidonic-acid/eicosanoid (EET/HETE) biotransformation and PPARalpha-associated lipid handling, consistent with AHR-linked CYP1A1/CYP1B1 induction. A secondary, weaker signal supports TP53-mediated cell-death transcription.

## Enrichment is not causation
Pathway enrichment is an over-representation statistic computed on a transcriptional snapshot: it shows annotated gene sets are over-represented among changed genes. It does not show those pathways drive or mediate the combination-treatment effect. Here the evidence is further limited: no term survives 5% FDR, overlaps are small (2-7 genes), and only 38 of 82 DE genes lie in the tested universe. Establishing causation would require targeted perturbation (e.g., CYP/PPARalpha/p53 inhibition) and phenotypic readouts.
