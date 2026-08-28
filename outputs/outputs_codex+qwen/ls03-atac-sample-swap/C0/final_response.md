Done. All four deliverables are written to `output/` and validated.

## Result: **swap detected — Cloaca ↔ Stomach** (confidence 0.9)

### How I reached the call
1. **Reconstructed coordinates.** The count table splits each GTF chromosome arm into ordered contigs (`chr1p_1`, `chr1p_2`…). Concatenating them in numeric order reproduces the AmexG v6 arm coordinates (verified on all 28 arms against the GTF's max gene positions).
2. **Promoter accessibility.** For every named gene I summed overlap-weighted ATAC signal over TSS ± 5 kb, converted to log2-CPM, and z-scored each gene across the 15 organ columns (variance-floored).
3. **Organ-marker coherence + pairwise swap scores** for all 105 pairs, scored by *mutual* cross-fit (`cross_min`) plus a directional win-rate — so a pair only wins if each tissue actually fits the other's column, not just by total library size or one-directional gain.

### Key evidence
- The **Cloaca column** expresses a gastric program: CLDN18.S z=+3.62, GKN1 +2.15, GKN2 +2.19, MUC5AC +2.63 (all argmax there).
- The **Stomach column** expresses a cloacal (hindgut + urothelium + genital-duct) program: CDX2 +2.39, MUC2 +3.34, FFAR4 +3.11, CALB1 +2.33, TCTE3/MSMB, UPK1B/UPK2.S.
- Cloaca↔Stomach is the **unique rank-1 pair**: cross-fit +0.81 (runner-up 0.32), win-rate 1.50/2 (75% of markers move in the swap direction), permutation p < 1e-4; robust across 2–10 kb promoter windows. All 13 other organs' marker sets peak in their own columns.

### Rejected alternatives (documented in report)
- **GallBladder↔Stomach** had the highest raw dsum but is strictly one-directional (cross_min 0.00) — an artifact of the weak GallBladder sample, exactly the library-size trap the notes warn about.
- **Gill↔Bladder** considered because the Gill column shows uroplakin/barrier genes, but the Gill column retains a genuine gill program (FOXI1, CCL28, DLX5) and the Bladder column shows no gill program while matching amphibian bladder (AQP-T2 z=+3.21 + detrusor smooth muscle).

### Deliverables
- `output/swap_call.json` — `swap_detected=true, organ_a=Cloaca, organ_b=Stomach, confidence=0.9, evidence=…`
- `output/sample_similarity.csv` — all 105 unordered pairs, lexicographically ordered names, finite scores, rank 1 = Cloaca/Stomach
- `output/analysis.py` — self-contained, regenerates both files (`python output/analysis.py`, ~1 min)
- `output/report.md` — full method, evidence tables, alternative-hypothesis analysis, and residual uncertainties (e.g., noisy ATP4B/TFF2 annotations, paralog averaging)