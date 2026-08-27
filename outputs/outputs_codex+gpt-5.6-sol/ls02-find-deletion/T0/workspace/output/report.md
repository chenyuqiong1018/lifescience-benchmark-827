# chr22 large-deletion call (T0)

## Result

The independent T0 analysis calls an approximately **1,000,000 bp deletion** between 100-kb-rounded chr22 boundaries **20,000,000** and **21,000,000**.

## Evidence

- Depth segmentation found 10 consecutive callable 100-kb bins with no uniquely mapped read bases. The immediately flanking bins contain 256,104 and 261,150 mapped bases.
- 2 correctly oriented FR pairs bridge the entire interval. Subtracting the deleted reference span restores ordinary library spans; the normal FR median is 500 bp.
- 4 otherwise-unmapped reads match exactly across the inferred left/right reference join.

The selected chromosome-analysis skill motivated explicit assembly/chromosome validation; genome-annotation guidance motivated unambiguous coordinate reporting; code-execution guidance motivated the self-contained, deterministic standard-library implementation.

## Precision

The required values are rounded to 100 kb. The zero-depth segmentation cannot independently support finer precision, even though exact junction reads are compatible with the same displayed boundaries. Under a 1-based interval convention, the removed bases are approximately 20,000,001 through 21,000,000.
