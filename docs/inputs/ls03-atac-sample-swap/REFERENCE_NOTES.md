# Reference and analysis notes

`AmexT_v47-AmexG_v6.0-DD.gtf.gz` is the official axolotl transcript annotation downloaded from axolotl-omics. It is supplied so that gene/promoter accessibility can be reconstructed against the same AmexG v6 coordinate system as the count table.

For `sample_similarity.csv`, write one row for every unordered organ pair using columns `organ_a,organ_b,swap_score,rank,evidence_type`. Names must be case-sensitive and lexicographically ordered within a row. `rank=1` is the most strongly supported swap; scores must be finite and larger values mean stronger support. Evidence may combine promoter accessibility, organ marker coherence, and comparison with public source profiles, but the final decision must be supported by a unique top pair rather than total library size alone.
