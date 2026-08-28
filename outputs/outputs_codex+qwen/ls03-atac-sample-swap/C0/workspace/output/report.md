# Axolotl bulk ATAC-seq organ-label swap detection — analysis report

**Date:** 2026-08-28  
**Call:** `swap_detected = true` — the labels **Cloaca** and **Stomach** are swapped (confidence 0.9).  
Deliverables: `output/swap_call.json`, `output/sample_similarity.csv`, `output/analysis.py`, this report.

---

## 1. Inputs

| File | Content |
|---|---|
| `inputs/sample.swap.atac.q1.tsv.gz` | 2,430,700 genomic bins (~10 kb tiling) × 15 organ columns of ATAC fragment counts |
| `inputs/sample.swap.atac.q1.chrom.sizes` | 63 contigs: each chromosome arm split into ordered parts (`chr1p_1`, `chr1p_2`, …) |
| `inputs/AmexT_v47-AmexG_v6.0-DD.gtf.gz` | Official axolotl annotation (AmexG v6); 99,088 genes, 22,695 with name symbols (`SYMBOL [nr]` / `SYMBOL [hs]`) |

Organs (column labels): Bladder, Brain, Cloaca, GallBladder, Gill, Heart, Intestine, Kidney, Limb, Liver, Lung, Pancreas, Prostate, Spleen, Stomach.

**Coordinate reconstruction.** The count table splits each GTF chromosome arm into numbered contigs. Concatenating contigs in numeric order reproduces the arm coordinate system: for all 28 arms, the contig-size sum matches the maximum GTF gene coordinate (max discrepancy ~1 Mb of a ~1.5 Gb arm). Bin positions were lifted to arm coordinates before intersecting with genes.

Library sizes (total bin counts) range from 3.98e8 (Liver) to 1.40e9 (Pancreas); all scores below are CPM-normalized, so library size alone cannot drive the call (per `REFERENCE_NOTES.md`).

## 2. Method (implemented in `output/analysis.py`)

1. **Promoter accessibility matrix.** For every named gene, accessibility = overlap-weighted ATAC counts over TSS ± 5 kb (strand-aware), converted to reads/bp, then CPM per organ (÷ library size × 1e6), then log2(x + 1e-3).
2. **Per-gene z-scores** across the 15 organ columns, with the standard deviation floored at the genome-wide median non-zero sd (0.531) so near-constant genes cannot dominate; genes with zero variance everywhere are treated as uninformative.
3. **Organ marker sets.** Curated sets of organ-specific genes using the human-ortholog symbols in `gene_name` (all annotated axolotl paralogs of a symbol are included and averaged — a conservative choice; see §6). Examples: stomach = GKN1, GKN2, CLDN18.S, PGA3, PGC, ANXA4, MUC5AC; cloaca (hindgut + urothelium + genital-duct receiving chamber) = CDX2, MUC2, TFF3.2, VIL1, UPK1B, UPK2.S, UPK3A, KRT13, FOXA2, PPARG, CALB1, SLC9A3, FFAR4; plus analogous sets for the other 13 organs.
4. **Pairwise swap scores** for all 105 unordered organ pairs (A, B), where MZ[O][c] = mean z of organ O's markers in column c:
   - directional gains `dA = MZ[A][col B] − MZ[A][col A]`, `dB` symmetric; `dsum = dA + dB`;
   - **`cross_min = min(MZ[A][col B], MZ[B][col A])`** — a true swap must place *both* tissues convincingly in the other's column; this rejects one-directional/library-size artifacts;
   - `win_rate` ∈ [0,2]: fraction of marker genes that move in the swap direction;
   - **`swap_score = cross_min + (win_rate − 1)`** (finite; larger = stronger support; rank 1 = best).
5. **Statistics.** Random-gene-set null for `dsum` (4,000 permutations, size-matched gene sets); robustness re-runs with promoter windows ±2 kb and ±10 kb.

## 3. Results

### 3.1 Organ self-consistency (marker-set argmax column)

Brain, Heart, Liver, Lung, Pancreas, Spleen, Intestine are unambiguously self-consistent (self z = +0.89…+2.28). Kidney column carries a full kidney program (proximal tubule SLC34A1 z=3.30, PAX8 2.87, SLC22A6, ACOT12 + distal ATP6V1B1 2.07). Limb, Gill, Bladder columns are also consistent with their labels (§5). The two clear outliers are reciprocal:

- **Stomach markers** (GKN1/GKN2/CLDN18.S/…) are *not* most accessible in the Stomach column (self z −0.35) — their argmax column is **Cloaca**.
- **Cloaca markers** (CDX2/MUC2/VIL1/UPK/PPARG/…) have argmax column **Stomach**.

### 3.2 Pairwise ranking (top of 105 pairs; full list in `sample_similarity.csv`)

| rank | pair | swap_score | cross_min | win_rate | dsum |
|---|---|---|---|---|---|
| **1** | **Cloaca ↔ Stomach** | **+1.311** | **+0.81** | **1.500** | **+1.44** |
| 2 | Bladder ↔ Stomach | +0.718 | +0.32 | 1.400 | +0.13 |
| 3 | Gill ↔ Lung | +0.305 | +0.37 | 0.939 | −0.77 |
| 4 | GallBladder ↔ Stomach | +0.303 | 0.00 | 1.300 | +1.54 |
| 5 | Prostate ↔ Stomach | +0.266 | −0.20 | 1.467 | +0.83 |

**Cloaca ↔ Stomach is the unique top pair on every decision-relevant metric:** largest bottleneck cross-fit (0.81 vs runner-up 0.32–0.37), largest win-rate (1.50, i.e. 75% of marker genes move in the swap direction), and near-top total improvement (dsum). It is the **only** pair combining strong positive cross-fit in *both* directions with win-rate ≥ 1.4. The result is invariant to the promoter window (rank 1 at ±2 kb, ±5 kb, ±10 kb with cross-fit 0.78–0.81).

Random-gene permutation of the directional gain: observed dsum 1.44 vs null mean 0.001, sd 0.198 → **z = 7.3, one-sided p < 1e-4**.

### 3.3 Decisive gene-level evidence

**The "Cloaca" column contains gastric tissue.** Stomach-specific genes peak in the Cloaca column (marker-level mean z; individual paralog maxima in parentheses):

| gene | function | Cloaca col | Stomach col |
|---|---|---|---|
| CLDN18.S | stomach-specific claudin | **+3.62** (argmax) | −0.50 |
| GKN1 | gastrokine (stomach-specific) | **+2.15** (paralogs 3.5) | −0.63 |
| GKN2 | gastrokine (stomach-specific) | **+2.19** | +0.11 |
| MUC5AC (primary locus) | gastric surface mucin | **+2.63** (argmax) | −0.04 |
| MUC1, ANXA4, CHGA, PGC | gastric epithelium | +1.80, +1.10, +0.71, +0.57 | +0.99, +1.45, +0.18, +0.50 |

plus gastric-wall smooth muscle/endoderm programs (ACTA2 +1.96, ACTG2 +2.15, ISL1 +2.09, GATA6 +1.74, SOX9 +1.71).

**The "Stomach" column contains cloacal tissue.** The cloaca receives the hindgut, ureters and genital ducts; exactly that composite program peaks in the Stomach column:

| gene | program | Stomach col | Cloaca col |
|---|---|---|---|
| CDX2 (CDX1 locus) | hindgut TF | **+2.39** (argmax) | +0.59 |
| MUC2 (primary locus) | goblet mucin | **+3.34** (argmax) | −0.52 |
| TCTE3, MSMB, TMPRSS3 | genital-duct/sperm duct | +3.25, +1.07, +1.84 | −0.00, −0.93, +0.23 |
| FFAR4, CALB1, SLC9A3 | hindgut absorption | +3.11, +2.33, +1.62 | +0.51, −0.10, −0.08 |
| UPK2.S, UPK1B | urothelium (cloacal/urinary) | +1.16, +1.11 | +0.06, −0.20 |
| PPARG, FOXA2, KRT18 | epithelial programs | +1.67, +0.58, +1.76 | +1.32, +2.37, +1.68 |
| MYOCD program (ACTA2/MYH11/CNN1) | cloacal smooth muscle | +1.73/+1.72/+1.19 | +1.96/+1.43/+1.41 |

Values are z-scores at the primary ortholog locus; where a marker symbol annotates multiple paralogs the scoring pipeline averages all of them (more conservative — the averaged scores are what enter the swap scores in §3.2).

Swapping the two labels restores both tissues to their correct columns (post-swap marker-set self-consistency: Stomach markers −0.35 → +0.81; Cloaca markers +0.72 → +1.00), while all 13 other organs are already self-consistent.

## 4. Alternative hypotheses considered and rejected

- **GallBladder ↔ Stomach** (highest raw dsum = +1.54). Rejected: strictly one-directional — stomach markers are flat in the GallBladder column (cross_min = 0.00), and the GallBladder column shows no gastric program (CLDN18.S −0.14, GKN1 +0.31). The GallBladder sample itself is the weakest in the panel (lowest marker self-consistency −0.55, lowest genome-wide correlations, ~5.7e8 reads): its generic epithelial markers (KRT8/18/19, EPCAM, MUC1, SOX9) "improve" in several epithelium-rich columns, inflating dsum without mutual fit. This is precisely the library-size/one-directional artifact that `REFERENCE_NOTES.md` warns against.
- **GallBladder ↔ Prostate.** One-directional (cross_min −0.27); prostate markers are negative in the GallBladder column. Rejected.
- **Gill ↔ Bladder.** The Gill column does show urothelial-barrier genes (UPK1A/1B/2/3A/3B z 1.8–3.0, KRT4 3.38, KRT76 3.34, GRHL3 2.28). However: (i) the Gill column simultaneously carries a genuine gill/branchial program (ionocyte TF FOXI1 +1.63, pharyngeal-arch patterning DLX5/MSX2, mucosal CCL28 +3.29, CNN2 +3.38); (ii) the Bladder column shows **no** gill program (globins z ≈ −0.04…−0.56, CA2 0.14; win-rate 0.65, cross_min −0.05); (iii) the Bladder column is consistent with genuine amphibian bladder: AQP-T2 (aquaporin water transport, the classical amphibian-bladder program) z = 3.21 argmax, detrusor smooth-muscle program (MYOCD-pathway ACTA2/ACTG2/MYH11), moderate uroplakin expression. Uroplakin/barrier-keratin expression in gill epidermis is a plausible amphibian biology (stratified barrier epithelium), whereas a Gill↔Bladder swap would leave the Stomach column's obvious hindgut/genital-duct program and the Cloaca column's gastric program unexplained. Rejected.
- **Prostate ↔ Stomach / Prostate ↔ Kidney.** The Prostate column's most prominent program (CA2 +2.92, SLC4A2 +3.00, RHCG +2.22, ATP6V1B1 +1.58, PAX2 +2.46, SLC12A3 +1.69) is an acid-base/ion-transport duct program consistent with Wolffian-duct-derived male reproductive tract (epididymal-type clear cells) of the urodele pelvic gland; stomach markers are low there (−0.30) and the Kidney column already contains a complete kidney program. No reciprocal fit. Rejected.
- **No swap.** Rejected: two columns (Stomach, Cloaca) are strongly inconsistent with their labels, and the inconsistency is exactly reciprocal and statistically significant (p < 1e-4).

## 5. Residual uncertainty (why confidence is 0.9, not 1.0)

- Some stomach markers are unreliable in this assembly/quantification: ATP4B and TFF2 have noisy paralog annotations (ATP4A absent; TFF2 maps to AB205-locus duplicates near TFF3.2), and PGA3 is nearly flat across organs; the gastric call therefore rests mainly on GKN1/GKN2/CLDN18.S/MUC5AC, which are, however, among the most stomach-specific genes known.
- TFF3.2 (intestinal trefoil) is high in the Cloaca column (+3.37); the axolotl TFF cluster on chr7p is ambiguous, but this slightly blurs the gastric interpretation.
- Symbol-level averaging sometimes dilutes true signals (e.g., the primary MUC2 locus is +3.34 in the Stomach column but its four MUC2-annotated paralogs pull the symbol average to −0.26; CDX2's true ortholog is CDX1 because CDX2 sits on an unplaced scaffold). The conservative averaging was kept for scoring; individual loci strengthen, not weaken, the call.
- The Gill column's UPK/barrier signal and the weak GallBladder sample are noted above; neither supports any alternative *unique* swap.

## 6. Reproducibility

Run from the workspace root:

```powershell
python output\analysis.py
```

It regenerates `output/sample_similarity.csv` (105 rows; columns `organ_a,organ_b,swap_score,rank,evidence_type`; names case-sensitive and lexicographically ordered within each row; rank 1 = best; larger score = stronger support) and `output/swap_call.json`. Runtime ≈ 1–2 min. Per `REFERENCE_NOTES.md`, evidence combines promoter accessibility and organ-marker coherence computed solely from the provided inputs; no external data were used.

## 7. Final call

```json
{
  "swap_detected": true,
  "organ_a": "Cloaca",
  "organ_b": "Stomach",
  "confidence": 0.9
}
```
