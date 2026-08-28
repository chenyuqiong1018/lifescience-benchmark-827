# Axolotl bulk ATAC-seq sample-swap analysis

## Question
Are two organ labels swapped in `sample.swap.atac.q1.tsv.gz` (15 organs x genome-wide 10 kb ATAC-seq bins, AmexG v6.0)?

## Answer

**swap_detected = true** - the labels of **Cloaca** and **Stomach** are swapped (confidence 0.969).

## Method

Swapping exactly two column labels leaves every label-free statistic of the matrix (e.g. the inter-column correlation matrix) unchanged, so detection requires external biological knowledge. The analysis uses organ marker genes:

1. Parse the AmexG v6.0 GTF; keep genes with a human ortholog symbol (`... [hs]` in `gene_name`; 8,996 genes) and take the TSS per strand.
2. Map base-chromosome coordinates (chr1p, chr2q, ...) onto the chunked assembly used by the count table (chr1p_1, chr1p_2, ...) via cumulative chrom sizes.
3. Promoter accessibility per gene = sum of the 10 kb bins overlapping TSS +/- 2 kb, normalized to counts-per-million (CPM) per organ. CPM is essential: library sizes range 398M (Liver) to 1,399M (Pancreas); without depth normalization the largest library dominates every score.
4. Standardize each gene across the 15 organs (z-score) and also compute a per-gene rank score; average both to reduce paralog/outlier effects.
5. Marker coherence matrix `T[organ o, column c]` = mean standardized accessibility of organ o's marker genes in column c. If labels are correct, `T[o,o]` is maximal for every o. If labels a<->b are swapped, `T[a,b] > T[a,a]` and `T[b,a] > T[b,b]`.
6. `swap_score(a,b) = (T[a,b]-T[a,a]) + (T[b,a]-T[b,b])`; all 105 unordered pairs are ranked (see `output/sample_similarity.csv`).

Marker sets are human ortholog symbols present in the axolotl annotation, chosen from prior knowledge of vertebrate organ identity. Decisive prior-knowledge facts: CLDN18/GKN1/CHIA/CCKAR are stomach-specific; CDX2/CDX4 are posterior (caudal) homeobox genes (the cloaca is the most caudal organ in the panel); REG4/ALPI mark hindgut/intestinal epithelium contained in the cloaca; CD8B/PRF1 mark lymphoid tissue, which the cloaca harbors (cloacal immune tissue).

## Results

### Marker coherence matrix T (rows = marker set, columns = data column)

```
                     Bladder  Brain  Cloaca  GallBladder  Gill  Heart  Intestine  Kidney  Limb  Liver  Lung  Pancreas  Prostate  Spleen  Stomach
Bladder-markers         0.80   0.17   -0.13         0.63  1.01  -0.44      -0.62   -0.61  0.18  -0.91 -0.40     -0.10      0.51   -0.95     0.88
Brain-markers          -0.03   1.16    0.16        -0.41  0.15   0.48      -0.13   -0.05  0.01  -0.53  0.16     -0.12     -0.30   -0.25    -0.30
Cloaca-markers         -0.18  -0.29    0.13         0.09 -0.03   0.33       0.63   -0.36 -0.28  -0.01  0.08     -0.11     -0.35   -0.42     0.78
GallBladder-markers     0.09  -0.48   -0.22         0.68  0.54  -0.75       0.35   -0.18 -0.01  -0.17  0.23      0.31      0.01   -0.69     0.28
Gill-markers            0.22  -0.44   -0.12         0.08  1.44  -0.61      -0.06   -0.37  0.05  -0.12  0.49     -0.35     -0.28   -0.84     0.92
Heart-markers           0.13   0.15    0.09        -0.36  0.17   1.03      -0.16   -0.14  0.14  -0.41 -0.00     -0.18     -0.26   -0.23     0.03
Intestine-markers       0.10  -0.31   -0.05        -0.06  0.01   0.13       1.35   -0.14  0.06   0.06 -0.05     -0.10     -0.41   -0.34    -0.26
Kidney-markers         -0.08  -0.14   -0.13        -0.20 -0.25   0.28      -0.03    1.21 -0.15   0.02 -0.08      0.05     -0.03   -0.41    -0.06
Limb-markers            0.23  -0.06    0.01        -0.29  0.65   0.27      -0.10    0.19  0.47  -0.57  0.04     -0.21      0.02   -0.43    -0.23
Liver-markers           0.11  -0.50    0.02        -0.05 -0.16  -0.30       0.04    0.03  0.04   1.06  0.12      0.06     -0.37   -0.19     0.08
Lung-markers            0.50  -0.14    0.22         0.03  0.37   0.76      -0.87    0.01  0.04  -0.40  0.94     -0.08     -0.56   -0.37    -0.45
Pancreas-markers       -0.11  -0.24    0.15        -0.06 -0.11   0.12       0.17   -0.18 -0.29  -0.19  0.00      1.71     -0.48   -0.37    -0.12
Prostate-markers       -0.11   0.13   -0.53        -0.51  0.20   0.52      -0.37    0.56  0.20  -0.68  0.12     -0.20      0.77   -0.30     0.20
Spleen-markers         -0.08  -0.25   -0.16        -0.59  0.21   0.40      -0.14   -0.33 -0.03   0.18 -0.09     -0.11     -0.30    1.44    -0.14
Stomach-markers        -0.23  -0.29    1.57         0.07  0.07  -0.33       0.12   -0.25  0.14  -0.21  0.10      0.30     -0.42   -0.44    -0.19
```

### Ranked organ pairs (top 10 of 105)

| rank | organ_a | organ_b | swap_score | d(a->b col) | d(b->a col) |
|---|---|---|---|---|---|
| 1 | Cloaca | Stomach | 2.427 | 0.658 | 1.769 |
| 2 | Bladder | Stomach | 0.043 | 0.078 | -0.035 |
| 3 | GallBladder | Stomach | -0.143 | -0.407 | 0.263 |
| 4 | Gill | Stomach | -0.259 | -0.521 | 0.261 |
| 5 | Limb | Stomach | -0.364 | -0.700 | 0.336 |
| 6 | Cloaca | Heart | -0.739 | 0.202 | -0.941 |
| 7 | Cloaca | Lung | -0.760 | -0.045 | -0.715 |
| 8 | Bladder | GallBladder | -0.761 | -0.171 | -0.590 |
| 9 | Prostate | Stomach | -0.795 | -0.569 | -0.226 |
| 10 | Bladder | Limb | -0.858 | -0.622 | -0.236 |

The top pair **Cloaca<->Stomach** scores **2.427**, i.e. 57x the second-best pair (Bladder<->Stomach, 0.043). Both directional deltas are positive and each organ's marker set peaks in the other organ's column (mutual-best match). No other pair has positive deltas in both directions with mutual-best columns.

### Key evidence genes

**Stomach marker genes peak in the `Cloaca` column** (they should peak in `Stomach` if labels were correct):

| gene | function | Bladder | Brain | Cloaca | GallBladder | Gill | Heart | Intestine | Kidney | Limb | Liver | Lung | Pancreas | Prostate | Spleen | Stomach |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CLDN18 | claudin-18 (stomach) | 1.2 | 0.8 | 8.5 | 1.1 | 1.2 | 1.2 | 1.2 | 0.9 | 1.3 | 0.9 | 0.9 | 1.1 | 1.0 | 1.0 | 0.9 |
| GKN1 | gastrokine-1 (stomach) | 0.5 | 0.3 | 7.1 | 0.4 | 0.5 | 0.4 | 1.1 | 0.4 | 0.5 | 0.5 | 0.7 | 0.6 | 0.4 | 0.3 | 0.7 |
| CHIA | acidic mammalian chitinase (stomach) | 0.6 | 0.6 | 6.9 | 0.8 | 0.7 | 0.0 | 0.8 | 0.5 | 0.8 | 1.0 | 0.7 | 1.0 | 0.4 | 0.5 | 0.7 |
| CHIA | acidic mammalian chitinase (stomach) | 0.8 | 0.8 | 7.0 | 1.1 | 1.1 | 0.0 | 1.0 | 1.0 | 1.2 | 1.4 | 1.3 | 1.2 | 1.0 | 1.3 | 1.2 |
| CHIA | acidic mammalian chitinase (stomach) | 0.8 | 0.8 | 7.0 | 1.1 | 1.1 | 0.0 | 1.0 | 1.0 | 1.2 | 1.4 | 1.3 | 1.2 | 1.0 | 1.3 | 1.2 |
| CCKAR | CCK-A receptor (stomach/pancreas) | 1.0 | 2.7 | 4.3 | 1.4 | 1.7 | 1.7 | 1.6 | 1.2 | 1.4 | 1.1 | 1.7 | 1.8 | 1.7 | 1.2 | 1.5 |

**Cloaca/caudal marker genes peak in the `Stomach` column:**

| gene | function | Bladder | Brain | Cloaca | GallBladder | Gill | Heart | Intestine | Kidney | Limb | Liver | Lung | Pancreas | Prostate | Spleen | Stomach |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CDX2 | caudal homeobox (posterior gut) | 0.4 | 0.6 | 0.7 | 0.5 | 0.6 | 0.6 | 3.5 | 0.4 | 0.5 | 0.4 | 0.5 | 0.5 | 0.5 | 0.4 | 4.0 |
| CDX4 | caudal homeobox (posterior gut) | 0.7 | 1.0 | 1.3 | 0.4 | 1.2 | 1.4 | 2.4 | 0.8 | 0.8 | 0.6 | 1.0 | 0.7 | 1.2 | 1.0 | 3.7 |
| CD8B | T-cell co-receptor (lymphoid) | 0.1 | 0.1 | 0.2 | 0.2 | 0.1 | 0.1 | 0.3 | 0.1 | 0.1 | 0.2 | 0.1 | 0.2 | 0.1 | 0.2 | 1.5 |
| REG4 | regenerating islet-derived 4 (hindgut) | 1.3 | 1.1 | 1.4 | 1.0 | 0.9 | 1.6 | 1.2 | 1.1 | 0.7 | 1.5 | 1.4 | 0.8 | 0.6 | 0.8 | 2.4 |
| ALPI | intestinal alkaline phosphatase | 1.2 | 0.9 | 3.7 | 1.3 | 1.8 | 1.4 | 2.5 | 1.1 | 1.4 | 1.5 | 2.1 | 1.3 | 1.4 | 1.2 | 4.5 |
| PRF1 | perforin (cytotoxic lymphoid) | 0.6 | 0.5 | 0.7 | 0.6 | 0.5 | 0.6 | 0.6 | 0.5 | 0.4 | 0.3 | 0.5 | 0.6 | 0.5 | 0.5 | 1.5 |

CLDN18 (a canonical stomach-specific tight-junction gene) is the single most discordant locus: its promoter is ~6-8x more accessible in the `Cloaca` column than in any other column, including `Stomach`. GKN1 (gastrokine-1), CHIA (gastric chitinase) and CCKAR show the same pattern. Conversely, the caudal homeobox genes CDX2/CDX4, the hindgut marker REG4 and lymphoid genes CD8B/PRF1 peak in the `Stomach` column - the cloaca is the most caudal organ and contains terminal hindgut and cloacal lymphoid tissue.

### Robustness checks

- Leave-one-locus-out: the Stomach-marker preference for the `Cloaca` column stays positive for every omitted locus (range 1.69..1.90).
- Leave-one-locus-out: the Cloaca-marker preference for the `Stomach` column stays positive for every omitted locus (range 0.54..0.81).
- Self-consistency: before correction 11/15 organs' marker sets peak in their own column; after relabeling Cloaca<->Stomach it is 13/15. Organs not involved in the swap with informative marker sets: 11/12 self-consistent.
- Library size is not the driver: scores are CPM-normalized and the decision uses reciprocal marker coherence, not total read counts (the Liver column has the smallest library yet scores highest for liver markers; Pancreas has the largest yet pancreas markers, not global signal, decide its column).

### Why no other pair is a credible swap

- A genuine swap requires BOTH directions: a's markers peak in b's column AND b's markers peak in a's column. No other pair satisfies this with positive deltas (see ranked table; e.g. Bladder->Gill is +1.01 but Gill->Bladder is -1.22).
- Limb markers partially match the Gill column because the limb set is skeletal-muscle/cartilage-centric (ACTA1, TNNC2, TNNI2, COL2A1, MYL1) and axolotl gills contain branchial muscle and cartilage; the reverse direction is ~0, so this is tissue composition, not a swap.
- The Bladder marker set has <3 informative loci in this annotation (UPK/KRT20 orthologs are absent or mis-assigned), so its row is near-noise and not evidence for or against any swap.

### Context: genome-wide profile correlations (label-free)

The correlation matrix is invariant to relabeling and therefore cannot by itself identify a swap; it is shown as biological context (Stomach and Cloaca profiles are moderately similar as endoderm-derived tissues, which makes the swap harder to notice by eye).

```
             Bladder  Brain  Cloaca  GallBladder  Gill  Heart  Intestine  Kidney  Limb  Liver  Lung  Pancreas  Prostate  Spleen  Stomach
Bladder         1.00   0.88    0.83         0.70  0.87   0.75       0.85    0.91  0.93   0.80  0.91      0.85      0.85    0.84     0.81
Brain           0.88   1.00    0.80         0.73  0.82   0.74       0.82    0.88  0.90   0.80  0.86      0.86      0.82    0.81     0.78
Cloaca          0.83   0.80    1.00         0.62  0.85   0.70       0.89    0.79  0.83   0.83  0.88      0.82      0.82    0.79     0.88
GallBladder     0.70   0.73    0.62         1.00  0.64   0.54       0.70    0.75  0.74   0.72  0.67      0.83      0.58    0.59     0.63
Gill            0.87   0.82    0.85         0.64  1.00   0.73       0.85    0.82  0.90   0.81  0.91      0.82      0.85    0.83     0.85
Heart           0.75   0.74    0.70         0.54  0.73   1.00       0.71    0.72  0.74   0.67  0.79      0.70      0.69    0.70     0.67
Intestine       0.85   0.82    0.89         0.70  0.85   0.71       1.00    0.83  0.86   0.87  0.87      0.87      0.82    0.80     0.90
Kidney          0.91   0.88    0.79         0.75  0.82   0.72       0.83    1.00  0.91   0.82  0.87      0.87      0.83    0.82     0.77
Limb            0.93   0.90    0.83         0.74  0.90   0.74       0.86    0.91  1.00   0.84  0.91      0.88      0.85    0.87     0.82
Liver           0.80   0.80    0.83         0.72  0.81   0.67       0.87    0.82  0.84   1.00  0.84      0.87      0.77    0.81     0.84
Lung            0.91   0.86    0.88         0.67  0.91   0.79       0.87    0.87  0.91   0.84  1.00      0.85      0.85    0.87     0.85
Pancreas        0.85   0.86    0.82         0.83  0.82   0.70       0.87    0.87  0.88   0.87  0.85      1.00      0.78    0.78     0.81
Prostate        0.85   0.82    0.82         0.58  0.85   0.69       0.82    0.83  0.85   0.77  0.85      0.78      1.00    0.82     0.81
Spleen          0.84   0.81    0.79         0.59  0.83   0.70       0.80    0.82  0.87   0.81  0.87      0.78      0.82    1.00     0.78
Stomach         0.81   0.78    0.88         0.63  0.85   0.67       0.90    0.77  0.82   0.84  0.85      0.81      0.81    0.78     1.00
```

## Conclusion

The `Cloaca` and `Stomach` labels are swapped: the column labeled `Cloaca` carries the stomach chromatin profile (GKN1/CLDN18/CHIA/CCKAR accessibility), and the column labeled `Stomach` carries the cloacal (caudal/hindgut + lymphoid) profile (CDX2/CDX4/REG4/CD8B). Relabeling Cloaca<->Stomach restores marker coherence for 13/15 organs. Confidence: 0.969.
