# T0 ATAC label-swap analysis

## Call

**Cloaca and Stomach are swapped** (high confidence). The top reciprocal coherence gain is 1.2527; runner-up GallBladder–Stomach is 0.6748; margin 0.5780.

Whole-arm GTF TSS coordinates were converted to sequential count-table fragments using chrom sizes. Three-bin promoter counts were normalized by full library totals, scaled across organs per marker, and all 105 unordered swaps were ranked. Thus total library size cannot itself drive the call. Regulatory-region and tissue-specific skill guidance shaped the promoter/marker evidence; code-execution guidance shaped the reproducible full-table implementation.

Conserved markers and internal reciprocity do not replace matched axolotl reference epigenomes. The deterministic rule requires a positive top gain and margin >0.25; otherwise `swap_detected=false`.
