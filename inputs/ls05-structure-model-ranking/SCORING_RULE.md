# Frozen local-extension rule

This fixture is a **benchmark-informed local extension**, not an upstream benchmark item.
The supplied values are already-computed comparison metrics against one frozen reference.
No coordinate-level claims may be made.

Rank every model by the following stable tuple:

1. complete chain mapping before incomplete mapping;
2. higher `tm_score`;
3. higher `lddt`;
4. lower `rmsd_a`;
5. more `aligned_residues`;
6. lower critical-region mean error; and
7. lexical `model_id` as the final tie break.

The critical region is residues 181--240. Its risk is the length-weighted mean of
`mean_error_a` over the overlapping rows in `residue_errors.csv`. In this fixture each
critical-region row covers the whole interval, so the reported risk equals that row's value.

For the required output columns use:

- `global_score = tm_score` (do not invent a weighted composite);
- `interface_score = 1` for complete chain mapping and `0` otherwise;
- `critical_residue_risk` in Å;
- `decision = preferred` for rank 1, `alternate` for other mapping-complete models,
  and `reject_incomplete_mapping` for an incomplete mapping.

Scientific basis: TM-score is a length-normalized global structure comparison measure
(Zhang & Skolnick, 2004, DOI 10.1002/prot.20264); lDDT is a superposition-free local
distance-difference measure (Mariani et al., 2013, DOI 10.1093/bioinformatics/btt473).
RMSD and aligned coverage are retained as secondary diagnostics rather than silently mixed
into a unitless score.

