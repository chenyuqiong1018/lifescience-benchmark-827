# ATAC-seq organ-label swap analysis

## Result

**Swap detected: yes**

**Pair: Cloaca ↔ Stomach**

**Confidence: high**

The top swap score is `1.252715`. The runner-up is
GallBladder ↔ Stomach at `0.674759`,
leaving a margin of `0.577955`. The leading pair is therefore distinct rather
than one of several nearly tied alternatives.

## Evidence

- Cloaca markers score 1.9756 in the Stomach-labelled sample versus 0.9943 in the Cloaca-labelled sample; Stomach markers score 1.1183 versus 0.8469.
- Leave-one-marker-locus-out re-ranking retained the same top pair in 99.3% of 150 trials; minimum top-versus-runner-up margin was 0.052243.
- Signals were normalized by complete-library CPM and expressed relative to each promoter locus's across-organ mean; total library size alone was not used as swap evidence.

## Method

The GTF gene TSS coordinates were mapped from whole chromosome arms onto the
split-contig coordinate system in the supplied chromosome-size file. For a
predeclared panel of organ identity markers, accessibility was summed over the TSS
bin and one adjacent 10-kb bin on each side. Each sample was normalized by its full
ATAC library count (CPM), then divided by that promoter locus's mean across the 15
organs so that each locus contributed a comparable relative-accessibility profile.
For every one of the 105 unordered organ pairs,
the `swap_score` is the gain in reciprocal marker-to-label coherence after exchanging
that pair. A deterministic leave-one-marker-locus-out analysis re-ranked all pairs
to assess robustness; it is an effect-size stability check, not a p-value.

The table contained 2,430,700 data rows with observed bin widths [1596, 2000, 2534, 3680, 4000, 4080, 4180, 4533, 4587, 4774, 4796, 4920, 5120, 5920, 5962, 6000, 6200, 6413, 6515, 6800, 6920, 7290, 7760, 7773, 8000, 8240, 8613, 8621, 8667, 9040, 10000];
410 requested promoter-region bins were present. Marker coverage:
Bladder=4 genes/5 loci, Brain=6 genes/11 loci, Cloaca=3 genes/3 loci, GallBladder=5 genes/12 loci, Gill=5 genes/7 loci, Heart=6 genes/17 loci, Intestine=5 genes/13 loci, Kidney=5 genes/29 loci, Limb=6 genes/10 loci, Liver=5 genes/7 loci, Lung=2 genes/2 loci, Pancreas=7 genes/11 loci, Prostate=3 genes/3 loci, Spleen=3 genes/3 loci, Stomach=3 genes/17 loci.

No figure was generated because the requested ranked CSV directly represents the
complete 105-pair comparison without loss of information.

## Reproduction

Run `python analysis.py` from any directory. The script locates the repository from
its own path and rewrites `swap_call.json`, `sample_similarity.csv`, and this report
deterministically.
