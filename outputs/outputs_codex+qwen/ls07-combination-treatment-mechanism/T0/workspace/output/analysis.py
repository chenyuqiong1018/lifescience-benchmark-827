#!/usr/bin/env python
"""LS07 combination-treatment mechanism analysis.

Contrast: Cisplatin_IC50_CBD_IC50 vs DMSO (frozen BixBench bix-43 inputs).
Pipeline: DESeq2 via PyDESeq2 0.5.0 -> DEG list -> GSEApy 1.1.4 enrichr ORA against
the evaluator-supplied frozen Reactome_2022 GMT, using the supplied background universe
(Reactome_2022.background.txt). No pathway library download/substitution and no
identifier mapping beyond the supplied ensg_to_gene_name.tsv.
"""
import hashlib
import json
import os
import platform
from collections import Counter

import numpy as np
import pandas as pd
import scipy
import gseapy as gp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

IN = "inputs"
OUT = "output"
os.makedirs(OUT, exist_ok=True)

COMPARISON = "Cisplatin_IC50_CBD_IC50 vs DMSO"
GROUP_TREAT = "Cisplatin_IC50_CBD_IC50"
GROUP_CTRL = "DMSO"
PADJ_CUT = 0.05
LFC_CUT = 1.0

# ---------------------------------------------------------------------------
# Compatibility shim: anndata >= 0.13 auto-reshapes 1-D varm arrays to (n, 1),
# which breaks PyDESeq2 0.5.0 (it indexes with 1-D boolean masks). Keep any
# value that entered validation as a 1-D ndarray as 1-D. Statistics unaffected.
# ---------------------------------------------------------------------------
import anndata._core.aligned_mapping as _am
_orig_validate = _am.AlignedMappingBase._validate_value
def _validate_keep_1d(self, val, key):
    was_1d = isinstance(val, np.ndarray) and val.ndim == 1
    out = _orig_validate(self, val, key)
    if was_1d and isinstance(out, np.ndarray) and out.ndim == 2 and out.shape[1] == 1:
        out = out.reshape(-1)
    return out
_am.AlignedMappingBase._validate_value = _validate_keep_1d

# ---------- integrity checks vs input manifest ----------
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

with open(os.path.join(IN, "Reactome_2022.manifest.json"), encoding="utf-8") as f:
    input_manifest = json.load(f)

gmt_path = os.path.join(IN, "Reactome_2022.gmt")
bg_path = os.path.join(IN, "Reactome_2022.background.txt")
map_path = os.path.join(IN, "ensg_to_gene_name.tsv")
gmt_sha, bg_sha, map_sha = sha256(gmt_path), sha256(bg_path), sha256(map_path)
assert gmt_sha == input_manifest["gmt_sha256"], "GMT sha256 mismatch vs frozen manifest"
assert bg_sha == input_manifest["background_sha256"], "background sha256 mismatch"
assert map_sha == input_manifest["mapping_sha256"], "mapping sha256 mismatch"
print("[integrity] GMT/background/mapping sha256 all match the frozen input manifest")

# ---------- load expression and layout ----------
counts = pd.read_csv(os.path.join(IN, "counts_raw_unfiltered.csv"), index_col=0)
layout = pd.read_csv(os.path.join(IN, "sample_layout.csv"))
sel = layout[layout["Group"].isin([GROUP_TREAT, GROUP_CTRL])].copy()
assert len(sel) == 6, "expected 3 vs 3 samples"
col_map = {sid: sid.replace("-", "_") for sid in sel["SampleID"]}
sel["column"] = sel["SampleID"].map(col_map)
assert sel["column"].isin(counts.columns).all(), "layout samples missing from counts"
print("[layout]", dict(zip(sel["SampleID"], sel["Group"])))

sub = counts[sel["column"].tolist()].copy()
sub.columns = sel["SampleID"].tolist()

# ---------- ENSG -> symbol, aggregate duplicates by sum ----------
mapping = pd.read_csv(map_path, sep="\t")
ensg2sym = dict(zip(mapping["ENSG"], mapping["gene_name"]))
sub["symbol"] = sub.index.map(lambda g: ensg2sym.get(g))
n_unmapped = int(sub["symbol"].isna().sum())
sub = sub.dropna(subset=["symbol"])
sub = sub.groupby("symbol").sum()
print(f"[mapping] {len(mapping)} ENSG rows; unmapped dropped: {n_unmapped}; symbols kept: {sub.shape[0]}")

# ---------- PyDESeq2 differential expression ----------
meta = pd.DataFrame({"condition": sel["Group"].tolist()}, index=sub.columns)
order = ["3-1", "3-2", "3-3", "9-1", "9-2", "9-3"]
meta = meta.loc[order]
sub = sub[order]
dds = DeseqDataSet(counts=sub.T, metadata=meta, design="~condition", refit_cooks=True)
dds.deseq2()
stat_res = DeseqStats(dds, contrast=["condition", GROUP_TREAT, GROUP_CTRL], quiet=True)
stat_res.summary()
res = stat_res.results_df.astype(float).copy()
res.to_csv(os.path.join(OUT, "deseq2_results.csv"))
print(f"[deseq2] tested genes: {res.shape[0]}; contrast: {GROUP_TREAT} vs {GROUP_CTRL}")

sig_full = res[(res["padj"] < PADJ_CUT) & (res["log2FoldChange"].abs() >= LFC_CUT)]
if len(sig_full) >= 10:
    degs = sig_full
    rule = f"padj<{PADJ_CUT} and |log2FC|>={LFC_CUT}"
else:
    degs = res[res["padj"] < PADJ_CUT]
    rule = f"padj<{PADJ_CUT} (|log2FC| filter relaxed: only {len(sig_full)} genes passed full cutoff)"
deg_list = degs.index.tolist()
n_up = int((degs["log2FoldChange"] > 0).sum())
n_down = int((degs["log2FoldChange"] < 0).sum())
print(f"[degs] {len(deg_list)} ({n_up} up / {n_down} down) by rule: {rule}")

# ---------- enrichment: GSEApy 1.1.4 enrichr on frozen local GMT + supplied universe ----------
background = [l.strip() for l in open(bg_path, encoding="utf-8") if l.strip()]
assert len(background) == input_manifest["background_gene_count"], "background size mismatch"
bg_set = set(background)
deg_bg = set(deg_list) & bg_set
print(f"[enrichr] DEGs in supplied universe: {len(deg_bg)}/{len(deg_list)}; universe size {len(bg_set)}")

enr = gp.enrichr(gene_list=deg_list, gene_sets=gmt_path, background=background,
                 outdir=None, no_plot=True, verbose=False)
tab = enr.res2d.copy()
tab["P-value"] = tab["P-value"].astype(float)
tab["Adjusted P-value"] = tab["Adjusted P-value"].astype(float)
tab = tab.sort_values(["Adjusted P-value", "P-value"], ascending=True).reset_index(drop=True)

# independent consistency checks vs the frozen GMT
gmt_dict = {l.rstrip("\n").split("\t")[0]: l.rstrip("\n").split("\t")[2:]
            for l in open(gmt_path, encoding="utf-8") if l.strip()}
assert len(gmt_dict) == input_manifest["term_count"], "GMT term count mismatch vs frozen manifest"
expected_overlap_terms = sum(1 for genes in gmt_dict.values() if set(genes) & deg_bg)
assert tab.shape[0] == expected_overlap_terms, (
    f"enrichment rows {tab.shape[0]} != independently counted overlapping terms {expected_overlap_terms}")
tab.to_csv(os.path.join(OUT, "pathway_enrichment.csv"), index=False)
sig_terms = tab[tab["Adjusted P-value"] < 0.05]
print(f"[enrichr] terms tested: {len(gmt_dict)}; terms with overlap (reported rows): {tab.shape[0]}; significant (BH<0.05): {len(sig_terms)}")
print(tab.head(15)[["Term", "Overlap", "P-value", "Adjusted P-value"]].to_string(index=False))

# ---------- theme classification for mechanism call ----------
THEMES = [
    ("Xenobiotic / drug metabolism (CYP450, Phase I/II)", ["cytochrome p450", "xenobiotic", "phase i", "phase ii", "eicosanoid", "hydroxyeicosatetraenoic", "epoxy", "dihydroxyeicosatrienoic", "arachidonic acid", "drug"]),
    ("DNA damage response / DNA repair", ["DNA damage", "DNA repair", "DNA replication", "homologous recombination", "nonhomologous end joining", "nucleotide excision", "base excision", "mismatch repair", "double-strand break", "checkpoint", "p53"]),
    ("Cell cycle / mitosis", ["cell cycle", "mitotic", "mitosis", "chromosome", "cytokinesis", "M phase", "G1/S", "G2/M", "segregation", "spindle", "centrosome"]),
    ("Apoptosis / programmed cell death", ["apoptosis", "apoptotic", "cell death", "caspase"]),
    ("Translation / ribosome", ["translation", "ribosom", "tRNA", "peptide", "nonsense-mediated", "cap-dependent", "IRES"]),
    ("Transcription / RNA processing / chromatin", ["transcription", "chromatin", "histone", "RNA polymerase", "splicing", "methylation", "epigenetic"]),
    ("Lipid / general metabolism", ["metabolism", "metabolic", "synthesis", "degradation", "catabolism", "biosynthesis", "fatty acid", "amino acid", "glycolysis", "TCA", "oxidative", "respiratory", "lipid", "cholesterol", "sterol", "glucose", "electron transport"]),
    ("Immune / signal transduction", ["interferon", "cytokine", "interleukin", "immune", "inflammasome", "NF-kB", "signaling", "signal transduction", "receptor", "kinase", "MAPK", "PI3K", "Wnt", "Notch", "Hedgehog", "PPAR"]),
    ("Protein homeostasis / trafficking / glycosylation", ["vesicle", "Golgi", "endosom", "lysosom", "autophagy", "trafficking", "secretion", "exocyt", "ubiquitin", "proteasom", "SUMOylation", "unfolded protein", "glycosylation", "glycan"]),
    ("Extracellular matrix / adhesion", ["extracellular matrix", "collagen", "integrin", "adhesion"]),
    ("Hemostasis / development", ["hemostasis", "coagulation", "platelet", "development", "differentiation"]),
]
def classify(term):
    t = term.lower()
    for theme, kws in THEMES:
        for kw in kws:
            if kw.lower() in t:
                return theme
    return "other"

TOP_N = 15
top_ranked = tab.head(TOP_N).copy()
top_ranked["theme"] = top_ranked["Term"].map(classify)
top_theme_counts = Counter(top_ranked["theme"])
top_theme, top_theme_n = top_theme_counts.most_common(1)[0]

sig_terms = sig_terms.copy()
sig_terms["theme"] = sig_terms["Term"].map(classify)
sig_theme_counts = Counter(sig_terms["theme"])

def term_record(row):
    return {"term": row["Term"], "adjusted_p_value": float(row["Adjusted P-value"]),
            "p_value": float(row["P-value"]), "overlap": row["Overlap"],
            "genes": row["Genes"]}

top10 = [term_record(r) for _, r in tab.head(10).iterrows()]
primary_term = tab.iloc[0]["Term"]

if len(sig_terms) > 0:
    evidence_basis = f"{len(sig_terms)} terms significant at BH<0.05"
else:
    evidence_basis = (
        f"no term reached BH<0.05; call based on top-ranked terms (best: adjusted p="
        f"{tab.iloc[0]['Adjusted P-value']:.3e})"
    )

mechanism = {
    "task": "ls07-combination-treatment-mechanism",
    "comparison": COMPARISON,
    "differential_expression": {
        "method": "DESeq2 via PyDESeq2 0.5.0 (Wald test, BH adjustment, Cooks outlier refit)",
        "design": "~condition",
        "contrast": ["condition", GROUP_TREAT, GROUP_CTRL],
        "samples": {"control": ["3-1", "3-2", "3-3"], "treatment": ["9-1", "9-2", "9-3"]},
        "gene_id_mapping": "provided ensg_to_gene_name.tsv (counts of duplicate symbols summed)",
        "tested_genes": int(res.shape[0]),
        "selection_rule": rule,
        "n_degs": int(len(deg_list)), "n_up": n_up, "n_down": n_down,
    },
    "enrichment": {
        "method": "Over-representation analysis via GSEApy 1.1.4 enrichr (local GMT, one-sided hypergeometric/Fisher, Benjamini-Hochberg)",
        "library": "frozen Reactome_2022 GMT (1818 terms); not downloaded or substituted",
        "background_universe": "provided Reactome_2022.background.txt (10489 genes)",
        "terms_tested": int(len(gmt_dict)),
        "terms_with_overlap_reported": int(tab.shape[0]),
        "terms_significant_BH_lt_0.05": int(len(sig_terms)),
    },
    "primary_mechanism": primary_term,
    "primary_mechanism_theme": top_theme,
    "evidence_basis": evidence_basis,
    "mechanism_description": (
        f"Best-supported primary cellular mechanism: '{primary_term}' (BH-adjusted p="
        f"{tab.iloc[0]['Adjusted P-value']:.3e}; {evidence_basis}). Dominant theme among the top "
        f"{TOP_N} ranked terms: '{top_theme}' ({top_theme_n}/{TOP_N})."
    ),
    "top_pathways": top10,
    "top_ranked_theme_counts": dict(top_theme_counts.most_common()),
    "significant_theme_counts": dict(sig_theme_counts.most_common()),
    "caveats": [
        "Pathway enrichment identifies annotated gene sets over-represented among DE genes; it is statistical association, not demonstrated causation.",
        "No Reactome term passed BH<0.05 in this contrast; the mechanism call reflects the strongest ranked enrichment evidence and is exploratory.",
        "The contrast measures the combined treatment's transcriptomic response vs DMSO; it does not separate the causal contribution of each drug.",
    ],
}
with open(os.path.join(OUT, "mechanism_call.json"), "w", encoding="utf-8") as f:
    json.dump(mechanism, f, indent=2, ensure_ascii=False)

# ---------- resource manifest ----------
manifest = {
    "resource_name": input_manifest["resource_name"],
    "source_manifest": "inputs/Reactome_2022.manifest.json",
    "gmt_file": "inputs/Reactome_2022.gmt",
    "gmt_sha256": gmt_sha,
    "gmt_sha256_matches_input_manifest": gmt_sha == input_manifest["gmt_sha256"],
    "background_file": "inputs/Reactome_2022.background.txt",
    "background_sha256": bg_sha,
    "background_sha256_matches_input_manifest": bg_sha == input_manifest["background_sha256"],
    "mapping_file": "inputs/ensg_to_gene_name.tsv",
    "mapping_sha256": map_sha,
    "mapping_sha256_matches_input_manifest": map_sha == input_manifest["mapping_sha256"],
    "term_count_in_gmt": int(len(gmt_dict)),
    "background_gene_count_observed": len(bg_set),
    "background_used_as_tested_universe": True,
    "tool": {"name": "gseapy", "version": gp.__version__, "entry_point": "gseapy.enrichr",
             "statistics": "one-sided hypergeometric (Fisher) test with Benjamini-Hochberg adjustment"},
    "environment": {"python": platform.python_version(), "numpy": np.__version__,
                    "pandas": pd.__version__, "scipy": scipy.__version__},
    "compatibility_note": "anndata>=0.13 reshapes 1-D varm arrays to (n,1); analysis.py applies a documented shim preserving 1-D masks for PyDESeq2 0.5.0 (no effect on statistics).",
    "policy_compliance": "No pathway library downloaded or substituted; only the frozen input GMT/background/mapping were used, per input manifest run_policy.",
    "run_date": "2026-08-30",
}
with open(os.path.join(OUT, "resource_manifest.json"), "w", encoding="utf-8") as f:
    json.dump(manifest, f, indent=2, ensure_ascii=False)

print("\n[done] artifacts written under output/: pathway_enrichment.csv, mechanism_call.json, resource_manifest.json, deseq2_results.csv (supplementary)")
print("[top term]", primary_term, "| adj.p =", f"{tab.iloc[0]['Adjusted P-value']:.3e}")
