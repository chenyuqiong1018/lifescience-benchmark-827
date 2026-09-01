# LS07 Combination Treatment Mechanism Report

## Objective

Identify the primary cellular mechanism of the approved combination treatment
`Cisplatin_IC50_CBD_IC50` relative to the `DMSO` vehicle control, using the
frozen BixBench bix-43 RNA-seq inputs and the evaluator-supplied frozen
`Reactome_2022` pathway library with its explicit background universe.

## Methods

Differential expression (DE) was run with pydeseq2 0.5.0 (DESeq2 model,
design `~condition`, 3 treated vs 3 control replicates) after mapping Ensembl
IDs to HGNC symbols. Significant genes (padj < 0.05, |log2FC| >= 1.0;
59 up, 25 down) were tested for pathway enrichment with GSEApy 1.1.4
(`enrichr`, offline Fisher exact test) against the byte-frozen
`Reactome_2022.gmt` (1818 terms), using the supplied
`Reactome_2022.background.txt` universe (10,489 genes). No library was
downloaded or substituted. 278 terms overlapped the DE list and were scored.

## Key results

- Most strongly induced genes: AHRR, CYP1A1, CYP1B1 (AhR/xenobiotic program),
  GDA, ALDH1A3, CDKN1A (p21).
- Top enriched pathways (by P-value): "Synthesis Of Epoxy (EET) And Dihydroxyeicosatrienoic Acids (DHET) R-HSA-2142670"
  (P = 3.91e-04), "Synthesis Of (16-20)-Hydroxyeicosatetraenoic Acids (HETE) R-HSA-2142816",
  "Cytochrome P450 - Arranged By Substrate Type", "Xenobiotics", and
  "Phase I - Functionalization Of Compounds". All are driven by the same
  markers: AHRR (log2FC +2.13, padj 9.8e-26), CYP1A1 (log2FC +2.03, padj 1.1e-19), CYP1B1 (log2FC +1.12, padj 5.2e-06).
- Minimum Adjusted P-value across terms: 0.070.

## Primary mechanism call

**Induction of aryl hydrocarbon receptor (AhR)-associated xenobiotic and
eicosanoid metabolism.** The treatment combination elicits a coherent
xenobiotic-response transcriptional program (AHRR/CYP1A1/CYP1B1 induction)
that simultaneously enriches cytochrome P450, Phase-I xenobiotic metabolism,
and arachidonic-acid epoxygenase/hydroxylase (EET/DHET, HETE) pathways.
Secondary observation: CDKN1A induction and p53-regulated cell-death terms
indicate a concurrent cell-cycle arrest/DNA-damage component consistent with
cisplatin exposure.

## Enrichment vs causation

Pathway enrichment is an over-representation statistic: it quantifies whether
DE genes cluster in curated Reactome annotations more than expected under the
supplied background. It is evidence of *association* between the
transcriptional response and pathway annotations, not evidence of causation.
These data do not establish that AhR/CYP-mediated metabolism causes the
combination's cytotoxic effect, nor that it is required for it; causal claims
would require perturbation (e.g., CYP or AhR inhibition/knockdown) and
phenotypic rescue experiments. Enrichment also inherits annotation bias and
depends on the arbitrary significance cutoffs used to define the DE gene list.

## Limitations

Three replicates per group; no pathway passed BH-adjusted significance at
0.05 (min adjusted P = 0.070), so the call rests on convergent
raw-P enrichment plus individually highly significant marker genes. Results
are conditioned on the frozen `Reactome_2022` snapshot and its background.
