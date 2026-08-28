## Result: a label swap was detected — **Cloaca ↔ Stomach** (confidence 0.969)

### Deliverables (all in `output/`)
| File | Content |
|---|---|
| `swap_call.json` | `swap_detected=true, organ_a=Cloaca, organ_b=Stomach, confidence=0.969` + evidence summary |
| `sample_similarity.csv` | All 105 unordered organ pairs, `organ_a,organ_b,swap_score,rank,evidence_type`, lexicographically ordered, finite scores, rank 1 = strongest |
| `analysis.py` | Self-contained, re-runnable pipeline (reads `./inputs`, writes `./output`) |
| `report.md` | Full methods, coherence matrices, evidence tables, robustness checks, limitations |

### How it was done
1. **Promoter accessibility**: parsed 8,996 AmexG v6 genes with human ortholog symbols (`[hs]` tags) from the GTF, mapped base-chromosome TSSs onto the chunked assembly via cumulative chrom sizes, and summed 10 kb bins within TSS±2 kb.
2. **CPM normalization was essential**: library sizes span 398M (Liver) → 1,399M (Pancreas); without depth correction the scores just track library size (I verified and corrected this).
3. **Marker coherence matrix** T[organ, column] using per-gene z + rank scores; `swap_score(a,b) = (T[a,b]−T[a,a]) + (T[b,a]−T[b,b])`. Since a two-label swap leaves the label-free correlation structure invariant, marker knowledge is what makes detection possible.

### Key evidence
- **Stomach markers peak in the `Cloaca` column**: CLDN18 (8.5 vs ≤1.3 everywhere else — the single most discordant locus), GKN1 (7.1), CHIA chitinases (6.9–7.0), CCKAR (4.3).
- **Caudal/cloacal markers peak in the `Stomach` column**: CDX2 (4.0 vs 0.7), CDX4 (3.7 vs 1.3), REG4, ALPI, and lymphoid CD8B/PRF1 (cloaca harbors terminal hindgut and cloacal lymphoid tissue).
- **Unique top pair**: swap_score 2.427 vs 0.043 for the second pair (~57× margin); mutual-best match holds; no other pair has both directions positive; leave-one-locus-out keeps both directional preferences positive; relabeling restores self-consistency from 11/15 to 13/15 organs.

### Documented uncertainties
The residual non-self-consistent organs are explained, not hidden: the Bladder set has <3 informative loci in this annotation (UPK/KRT20 orthologs missing/mis-assigned), and Limb markers partially match Gill because gills genuinely contain skeletal muscle/cartilage (ACTA1, TNNC2, COL2A1) — its reverse direction is ~0, so neither forms a competing swap candidate.