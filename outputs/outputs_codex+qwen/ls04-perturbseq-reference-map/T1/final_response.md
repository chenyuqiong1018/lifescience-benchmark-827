Done. All deliverables are in `output/` and validated.

## Headline answer
| target_gene | query_guide_id | score | runner_up_score | confidence |
|---|---|---|---|---|
| **PABPC1** | **guide18** | 0.532 | 0.339 | 0.035 |
| **NUDT21** | **guide13** | 0.636 | 0.513 | 0.038 |
| **LEO1** | **guide14** | 0.314 | 0.223 | 0.029 |

## Method (one primary path)
Query and reference are raw-count Perturb-seq datasets in shifted cellular contexts (NT-control profiles correlate only r = 0.92). So matching is done on **perturbation effects**, never absolute expression: each dataset is normalized independently (CPM 1e4 + log1p), each guide's signature is its mean minus its **own dataset's** NT-control mean over 11,858 shared genes, scored by Pearson correlation, and assigned one-to-one via Hungarian matching (reference has exactly one guide per gene).

## Leakage prevention
- Runtime assertion that query `obs` contains only the anonymized `guide` column (no target metadata exists to leak).
- No cross-dataset name matching and no parsing of guide-name strings for gene identity (only the explicit NT control designation is used).
- Within-dataset baselines only; reference guides aggregated to genes strictly via the `target_gene` label column.

## Validation & calibrated uncertainty
- **PABPC1→guide18, NUDT21→guide13**: mutually consistent in both directions, positive margins, stable when all signatures are re-derived on a random 50% subsample of cells — high confidence. NUDT21's documented confounder is CPSF6 (same pathway).
- **LEO1→guide14 is tentative**: guide14 scores CTR9 slightly higher from the query side (0.370 vs 0.314; both are PAF1-complex subunits with near-degenerate knockdown effects), and the LEO1 row flips under subsampling. Best read as a PAF1-complex-level assignment.
- Ambiguity quantified in `report.md`: 23/34 rows have negative margins (same-pathway signature degeneracy), 76.5% subsample stability (all 8 unstable rows in the PAF1/CPSF low-margin blocks), NT controls at noise floor (max |r| = 0.14), and a flagged unresolvable PHF3/SCAF8/guide32 block (query carries 32 perturbed guides vs 33 reference genes).

Deliverables: `output/guide_mapping.csv` (exactly 34 rows × 5 required columns, schema-asserted), `output/analysis.py` (reproducible, seeded), `output/report.md`, plus supplementary `output/supplementary_query_guide_mapping.csv` and `output/run_summary.json`.