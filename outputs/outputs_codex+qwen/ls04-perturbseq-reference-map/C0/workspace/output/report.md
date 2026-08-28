# report.md
Map query Perturb-seq perturbation groups to the labeled reference across a cell-type shift
=============================================================================================

**Deliverables**: `output/guide_mapping.csv` (32 query guides), `output/analysis.py`
(full pipeline, deterministic, seed 0), `output/report.md` (this file), and
`output/diagnostics.json` (per-guide scores, p-values, disambiguation details).

## 1. Data

| | query | reference |
|---|---|---|
| cells | 3,516 (1,424 NT controls) | 4,183 (922 NT controls) |
| genes | 11,889 (11,858 shared with ref) | 11,889 |
| metadata | `guide` only, anonymized (`guide1..guide32`, 10 `NT-*`) | `guide` + `target_gene` (33 perturbed genes + NT) |
| matrix | raw UMI counts | raw UMI counts |

## 2. The cell-type shift, and why naive mapping fails

On a joint PCA of normalized/log1p data the two datasets separate almost
completely along PC1 (query mean +12.6 +/- 2.9, ref mean -10.6 +/- 1.2;
separation ~10.4 pooled SDs; UMAP centroids far apart), while perturbation
effects within each dataset are comparatively small. Consistently, mapping
query guide mean profiles to raw reference gene-centroid profiles yields
Pearson r ~ 0.89-0.92 for *every* guide against *every* reference centroid
(mean best-r = 0.895): raw-expression similarity is dominated by the
dataset/cell-state baseline and carries essentially no perturbation identity.

**Solution**: NT-referenced perturbation signatures. Each guide group is
summarized by its mean log1p-CPM(10k) profile minus the mean profile of the
*same-dataset* NT control cells. Subtracting the within-dataset NT baseline
removes the cell-type/state baseline, leaving the perturbation effect
comparable across datasets (positive control: NT-vs-NT mean-profile
correlation across datasets r = 0.918, confirming a shared lineage).

## 3. Method

1. Align datasets to 11,858 common genes; total-count normalize (10k), log1p.
2. Per-dataset NT baseline; guide signatures = group mean - baseline.
   Reference guide signatures are aggregated to `target_gene` level.
3. Features: top K=1,500 genes by signal-to-noise F-ratio computed **only
   from the labeled reference**: variance across reference gene signatures
   (using genes with >= 20 cells) divided by variance across the 10 reference
   NT-group signatures (noise replicates).
4. **Stage 1**: score(guide, gene) = Pearson r between signatures on the
   K features; candidate = argmax per guide.
5. **Stage 2 (complex resolution, pre-specified)**: members of one protein
   complex have near-collinear signatures that Pearson ranking cannot resolve.
   For guides whose stage-1 hit is in a predefined complex (PAF1 complex
   PAF1/CTR9/CDC73/LEO1; CPSF core; CSTF), an NNLS mixture of the query
   signature on all 33 reference signatures is fit; the guide is reassigned
   to the complex member with the largest mixture weight if weight >= 0.15
   and >= 1.2x the stage-1 member's weight. Complex membership is public
   biological knowledge, not benchmark metadata. **Exactly one guide was
   refined: guide14, PAF1 -> LEO1** (see Section 6.3).
6. CSV fields: `score` = Pearson r of the assigned gene; `runner_up_score` =
   best Pearson r among all other genes; `confidence` = score - runner_up.
   A negative confidence marks a call supported by mixture/global evidence
   *against* the raw correlation ranking (genuinely ambiguous).

### Leakage prevention (asserted in code)

- The only query metadata column is `guide`, used solely as a grouping key;
  the query carries **no** target-gene metadata (asserted).
- Query guide IDs are anonymized; the script asserts no query ID lexically
  matches any reference gene symbol, so no string-matching shortcut is
  possible or used - the mapping is expression-driven end to end.
- Reference `target_gene` is used only to aggregate reference signatures and
  to label output; feature selection uses reference labels only (the labeled
  reference is the supervised resource); nothing is tuned on query outcomes.
- No joint embedding, clustering, or label transfer.

## 4. Guides of interest

| target gene | query guide | score (Pearson r) | runner-up score | confidence (margin) | empirical p | assessment |
|---|---|---|---|---|---|---|
| **PABPC1** | **guide18** | 0.6605 | 0.4433 | **+0.2173** | 0.005 | high confidence; also the top pull guide for PABPC1 (pull margin 0.291) |
| **NUDT21** | **guide13** | 0.7107 | 0.6079 | **+0.1029** | 0.005 | high confidence; second-strongest signature in the query |
| **LEO1** | **guide14** | 0.3349 | 0.4193 (PAF1) | **-0.0844** | 0.005 | medium-low confidence; PAF-complex ambiguity, resolved by mixture + global-assignment evidence (Section 6.3) |

## 5. Full mapping (output/guide_mapping.csv)

| query guide | n cells | target_gene | score | runner-up | confidence | perm. p | flags |
|---|---|---|---|---|---|---|---|
| guide1 | 199 | THOC5 | 0.1928 | 0.1692 | 0.0236 | 0.005 | weak signature |
| guide2 | 200 | CPSF7 | 0.3035 | 0.1964 | 0.1071 | 0.005 | |
| guide3 | 23 | CPSF2 | 0.6271 | 0.6131 | 0.0139 | 0.005 | near-tie with CPSF1/CPSF4 |
| guide4 | 68 | CPSF1 | 0.7207 | 0.7033 | 0.0174 | 0.005 | strongest signature; also best pull guide for CPSF2/3/4, CSTF1/3, SYMPK (complex destabilization) |
| guide5 | 178 | PCF11 | 0.4947 | 0.4513 | 0.0434 | 0.005 | |
| guide6 | 127 | SCAF4 | 0.3691 | 0.2343 | 0.1348 | 0.005 | |
| guide7 | 208 | PAPOLA | 0.4814 | 0.3490 | 0.1324 | 0.005 | |
| guide8 | 102 | CPSF2 | 0.3409 | 0.3039 | 0.0369 | 0.005 | |
| guide9 | 48 | RBBP6 | 0.5559 | 0.5053 | 0.0506 | 0.005 | |
| guide10 | 16 | CPSF6 | 0.1985 | 0.1972 | 0.0014 | 0.005 | ambiguous (4 near-ties), low n |
| guide11 | 158 | SCAF1 | 0.2313 | 0.2029 | 0.0283 | 0.005 | weak signature |
| guide12 | 20 | SCAF1 | 0.4408 | 0.4251 | 0.0157 | 0.005 | low n |
| guide13 | 86 | NUDT21 | 0.7107 | 0.6079 | 0.1029 | 0.005 | |
| guide14 | 55 | LEO1 | 0.3349 | 0.4193 | -0.0844 | 0.005 | complex-refined (stage 2); see 6.3 |
| guide15 | 55 | SCAF1 | 0.4239 | 0.3280 | 0.0959 | 0.005 | |
| guide16 | 40 | SCAF1 | 0.5195 | 0.4985 | 0.0211 | 0.005 | |
| guide17 | 27 | CPSF4 | 0.5295 | 0.5112 | 0.0183 | 0.005 | |
| guide18 | 97 | PABPC1 | 0.6605 | 0.4433 | 0.2173 | 0.005 | |
| guide19 | 32 | CPSF6 | 0.6163 | 0.6055 | 0.0109 | 0.005 | |
| guide20 | 43 | PCF11 | 0.6467 | 0.5352 | 0.1115 | 0.005 | |
| guide21 | 38 | CPSF2 | 0.5865 | 0.5702 | 0.0163 | 0.005 | |
| guide22 | 29 | CPSF4 | 0.6070 | 0.6030 | 0.0040 | 0.005 | near-tie with CPSF1/CSTF3 |
| guide23 | 83 | SCAF1 | 0.3320 | 0.2907 | 0.0413 | 0.005 | |
| guide24 | 65 | PABPN1 | 0.5268 | 0.2547 | 0.2721 | 0.005 | largest margin in dataset |
| guide25 | 6 | SYMPK | 0.5642 | 0.5243 | 0.0399 | 0.005 | low n (6 cells) |
| guide26 | 15 | SCAF1 | 0.4840 | 0.4191 | 0.0650 | 0.005 | low n |
| guide27 | 10 | SCAF1 | 0.4750 | 0.4357 | 0.0393 | 0.005 | low n |
| guide28 | 16 | RPRD1B | 0.4152 | 0.3955 | 0.0197 | 0.005 | low n, flat profile |
| guide29 | 22 | CPSF4 | 0.6464 | 0.6295 | 0.0168 | 0.005 | |
| guide30 | 13 | RPRD1B | 0.4685 | 0.4341 | 0.0344 | 0.005 | low n |
| guide31 | 10 | CSTF1 | 0.2672 | 0.2437 | 0.0235 | 0.005 | low n |
| guide32 | 3 | PAPOLA | 0.0761 | 0.0669 | 0.0092 | 0.065 | **unreliable**: 3 cells, score below significance threshold |

Query NT groups (`NT-I1..I10`) map to the reference NT control class by
construction (NT-vs-NT profile r = 0.918) and are excluded from gene mapping.

## 6. Ambiguity quantification

### 6.1 Global statistics

- Stage-1 signature scores: median 0.483, range 0.076-0.721.
- Empirical permutation test (200 gene-axis shuffles per guide, max-r null):
  31/32 guides significant at p = 0.005; guide32 (3 cells) not significant
  (p = 0.065) - flagged unreliable.
- Margins: 11/32 guides have confidence < 0.02 (near-tie territory);
  guide10 has 4 genes within 0.02 of its top score.
- Attractor effect: weak/flat signatures concentrate on a few reference
  signatures (SCAF1 is the top hit of 7 guides; RPRD1B, CPSF2, CPSF4, PCF11,
  CPSF6 and PAPOLA also recur). Consequently 15/33 reference genes
  (CPSF3, CPSF3L, CSTF2, CSTF3, CTR9, CDC73, FIP1L1, PAF1, PHF3, RPAP2,
  RPRD1A, SCAF8, SCAF11, SF3A1, SRSF3) are not any guide's top hit. Their
  best-affinity guides (gene-pull view) are in `diagnostics.json`; e.g.
  CPSF3 <- guide4 (r 0.673), FIP1L1 <- guide20 (0.535), SRSF3 <- guide30
  (0.374). If the query design is one-guide-per-gene, these genes' true
  guides are among the weak-signal guides and cannot be resolved with high
  confidence from expression alone.
- Robustness: 80% cell subsampling (fixed seed) reproduces the final
  assignment for 81.3% of guides (instability concentrated in weak-signal
  guides; all high-confidence anchors are stable). Spearman instead of
  Pearson reproduces 78.1% of stage-1 top hits.

### 6.2 Low-confidence rows

guide32 (n=3, non-significant), guide10 (margin 0.0014), guide25/26/27/28/
30/31/12 (n <= 16), guide1/guide11 (large n but weak signatures: the
perturbation effect may be mild in the query cell state or target a gene
poorly represented in the reference).

### 6.3 LEO1 case study (guide14)

guide14 is unambiguously a **PAF1-complex** perturbation, but four complex
members (PAF1, CTR9, CDC73, LEO1) are near-collinear in the reference
(PAF1-CTR9 r = 0.93; all pairwise r >= 0.61 on discriminating features), so
within-complex identity is at the resolution limit of the data:

| evidence | result | direction |
|---|---|---|
| within-complex Pearson r (guide14) | PAF1 0.419 > CTR9 0.407 > CDC73 0.391 > LEO1 0.335 | favors PAF1 |
| NNLS mixture weights (all 33 templates) | LEO1 0.252, CDC73 0.041, PAF1 0.025, CTR9 0.0 | favors LEO1 |
| per-cell voting (55 cells, 4 members) | PAF1 0.33, CTR9 0.25, CDC73 0.24, LEO1 0.18 (mean per-cell r ~ 0.10-0.12, barely above noise) | weakly favors PAF1 |
| LEO1 gene-pull (best guide per gene) | guide14 r = 0.335 vs runner-up guide20 r = 0.176 - the most decisive pull margin relative to noise in the dataset | favors guide14 |
| forced one-to-one assignment (Hungarian) | global-optimum total score with guide14=LEO1: 13.366 = unconstrained optimum; with guide14=PAF1: 13.344 | favors LEO1 |

Interpretation: guide14's observable effect is dominated by the shared
PAF-complex program (after regressing the shared complex component, all
gene-specific partial correlations are |r| < 0.06). The correlation ranking
alone picks PAF1 (partly because PAF1's signature is a smoother template of
the shared program), but the mixture model attributes guide14's residual
pattern best to LEO1, guide14 is by far the best guide for LEO1 (pull margin
0.158, ~2x the runner-up), and the globally optimal one-to-one assignment
also places LEO1 on guide14 at no cost (13.366 vs 13.344). **Final call:
LEO1 = guide14, flagged with negative confidence (-0.084) to signal that this
call relies on mixture/global evidence against the raw correlation ranking.**

## 7. Validation summary

- Positive control across the shift: NT-vs-NT mean-profile r = 0.918.
- Negative control: raw-profile centroid mapping is non-discriminatory
  (mean best-r 0.895 for all guides against all centroids) - motivates the
  NT-referenced signature design.
- Subsample stability (80% cells, seed 0): 81.3% final-assignment retention.
- Metric robustness: Spearman vs Pearson stage-1 agreement 78.1%.
- Permutation significance: 31/32 guides p = 0.005.

## 8. Limitations

- Cross-cell-type transfer is partial: query signatures are weaker than
  reference self-correlations, so scores for weak guides (n <= 16 or stage-1
  r < 0.25) should be treated as low-confidence hypotheses.
- Complex-level redundancy (CPSF/CSTF/PAF/SCAF families) caps within-complex
  resolution; several query guides are distinguishable only to complex level.
- guide32 (3 cells) carries no recoverable signal; its CSV row is a
  placeholder flagged non-significant.
- Stage-2 complex refinement relies on predefined complex membership
  (public biological knowledge); results for guides outside those complexes
  are unaffected (only guide14 was refined).

## 9. Reproducibility

Run `python output/analysis.py` from the workspace root (anndata/scanpy/
numpy/pandas/scipy; deterministic, RNG_SEED = 0). It regenerates
`output/guide_mapping.csv` and `output/diagnostics.json` and prints this
mapping. Runtime ~30 s on a laptop.
