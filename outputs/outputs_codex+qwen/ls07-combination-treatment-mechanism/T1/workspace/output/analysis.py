#!/usr/bin/env python
"""ls07-combination-treatment-mechanism / T1: approved DE + enrichment pipeline.

Comparison : Cisplatin_IC50_CBD_IC50 (combination treatment) vs DMSO (vehicle)
DE method  : PyDESeq2 on raw counts (6 samples: 3 per group)
Enrichment : GSEApy 1.1.4 `enrichr` in LOCAL mode against the evaluator-supplied
             frozen `Reactome_2022` GMT, using the supplied background universe.
             One-sided hypergeometric (Fisher) test, Enrichr odds-ratio and
             combined score, Benjamini-Hochberg FDR correction.

Fully offline: no pathway library is downloaded or substituted, and no
identifier mapping other than the supplied ensg_to_gene_name.tsv is used.
GSEApy reports only terms sharing >=1 gene with the query list (native
behaviour of calc_pvalues); those are the terms tested.
"""
import hashlib
import json
import os
import platform
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import gseapy as gp
import pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

BASE = os.path.dirname(os.path.abspath(__file__))          # output/
ROOT = os.path.dirname(BASE)                               # workspace
INP = os.path.join(ROOT, "inputs")

COUNTS_CSV = os.path.join(INP, "counts_raw_unfiltered.csv")
LAYOUT_CSV = os.path.join(INP, "sample_layout.csv")
MAP_TSV = os.path.join(INP, "ensg_to_gene_name.tsv")
GMT = os.path.join(INP, "Reactome_2022.gmt")
BG = os.path.join(INP, "Reactome_2022.background.txt")
IN_MANIFEST = os.path.join(INP, "Reactome_2022.manifest.json")

CONTROL, TREAT = "DMSO", "Cisplatin_IC50_CBD_IC50"
MIN_TOTAL_COUNT = 10        # minimal pre-filter across the 6 samples
PADJ_CUTOFF = 0.05
L2FC_CUTOFF = 1.0

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def parse_gmt(path):
    terms = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            terms[parts[0]] = [g for g in parts[2:] if g]
    return terms

# ---------------------------------------------------------------- inputs
counts = pd.read_csv(COUNTS_CSV, index_col=0)
layout = pd.read_csv(LAYOUT_CSV)
mp = pd.read_csv(MAP_TSV, sep="\t")
in_manifest = json.load(open(IN_MANIFEST, encoding="utf-8"))

sel = layout[layout["Group"].isin([CONTROL, TREAT])].copy()
assert len(sel) == 6, "expected 3 vs 3 samples"
sel["col"] = sel["SampleID"].str.replace("-", "_")
missing = [c for c in sel["col"] if c not in counts.columns]
assert not missing, f"missing count columns: {missing}"

sub = counts[sel["col"].tolist()].copy()
meta = pd.DataFrame({"condition": [TREAT if g == TREAT else CONTROL
                                   for g in sel["Group"]]}, index=sel["col"])

keep = sub.sum(axis=1) >= MIN_TOTAL_COUNT
sub = sub.loc[keep]
print(f"[de] genes after count filter: {sub.shape[0]}", flush=True)

# ---------------------------------------------------------------- DE
dds = DeseqDataSet(counts=sub.T.astype(int), metadata=meta, design="~condition",
                   quiet=True)
dds.deseq2()
ds = DeseqStats(dds, contrast=["condition", TREAT, CONTROL], quiet=True)
ds.summary()
res = ds.results_df.copy()
res.index.name = "ensembl_id"

ensg2sym = dict(zip(mp["ENSG"].astype(str), mp["gene_name"].astype(str)))
res["gene_symbol"] = [ensg2sym.get(str(i), np.nan) for i in res.index]
res = res.dropna(subset=["gene_symbol"])
res = res[res["gene_symbol"] != "nan"]
res = res.sort_values(["padj", "pvalue", "baseMean"],
                      ascending=[True, True, False], na_position="last")
res = res.drop_duplicates(subset="gene_symbol", keep="first")
res = res.sort_values(["padj", "pvalue"], ascending=[True, True],
                      na_position="last").reset_index()

de = res[(res["padj"] < PADJ_CUTOFF) &
         (res["log2FoldChange"].abs() >= L2FC_CUTOFF)].copy()
up = de.loc[de["log2FoldChange"] > 0, "gene_symbol"].tolist()
down = de.loc[de["log2FoldChange"] < 0, "gene_symbol"].tolist()
de_genes = de["gene_symbol"].tolist()
print(f"[de] DE genes: {len(de_genes)} (up {len(up)} / down {len(down)})", flush=True)

res.to_csv(os.path.join(BASE, "de_results.csv"), index=False)

# ------------------------------------------------- resource verification
terms = parse_gmt(GMT)
union = sorted({g for genes in terms.values() for g in genes})
background = [l.strip() for l in open(BG, encoding="utf-8") if l.strip()]
assert sorted(set(background)) == background, "background not sorted unique"
assert union == sorted(background), "GMT union != supplied background universe"

sha = {
    "gmt": sha256_file(GMT),
    "background": sha256_file(BG),
    "mapping": sha256_file(MAP_TSV),
}
for key, field in [("gmt", "gmt_sha256"), ("background", "background_sha256"),
                   ("mapping", "mapping_sha256")]:
    assert sha[key] == in_manifest[field], f"sha256 mismatch for {key}"
assert len(terms) == in_manifest["term_count"]
assert len(background) == in_manifest["background_gene_count"]
print(f"[res] verified: {len(terms)} terms, {len(background)} background genes", flush=True)

bg_set = set(background)
query_in_bg = [g for g in de_genes if g in bg_set]
print(f"[enr] query genes in background: {len(query_in_bg)}/{len(de_genes)}", flush=True)

# ------------------------------------------------------------- enrichment
enr = gp.enrichr(gene_list=de_genes, gene_sets=GMT, background=BG,
                 organism="human", outdir=None, no_plot=True)
tab = enr.res2d.copy()
if "Gene_set" in tab.columns:
    tab = tab.drop(columns=["Gene_set"])
tab["P-value"] = tab["P-value"].astype(float)
tab["Adjusted P-value"] = tab["Adjusted P-value"].astype(float)
tab = tab.sort_values("P-value", ascending=True).reset_index(drop=True)
tab.to_csv(os.path.join(BASE, "pathway_enrichment.csv"), index=False)
print(f"[enr] terms tested (>=1 query gene): {len(tab)}", flush=True)

sig = tab[tab["Adjusted P-value"] < PADJ_CUTOFF]
print(f"[enr] significant (BH < {PADJ_CUTOFF}): {len(sig)}", flush=True)

# ------------------------------------------------------- mechanism call
# (family label, meta-mechanism, keywords) - first match wins
FAMILIES = [
    ("Xenobiotic/drug metabolism (cytochrome P450)",
     "Cytochrome P450-mediated xenobiotic/drug and lipid metabolism",
     ["cytochrome p450", "xenobiotic", "phase i", "phase ii",
      "conjugation", "cyp"]),
    ("Arachidonic acid / eicosanoid metabolism",
     "Cytochrome P450-mediated xenobiotic/drug and lipid metabolism",
     ["arachidonic", "eicosanoid", "epoxy", "dihydroxyeicosatrienoic",
      "hete", "prostaglandin", "leukotriene", "lipoxin"]),
    ("Lipid metabolism / PPAR signaling",
     "Cytochrome P450-mediated xenobiotic/drug and lipid metabolism",
     ["lipid", "ppar", "fatty acid", "sterol", "cholesterol",
      "triglyceride", "ketone"]),
    ("DNA damage response / DNA repair",
     "DNA damage response and genome maintenance",
     ["dna repair", "dna damage", "dna replication", "mismatch repair",
      "homologous recombination", "nucleotide excision", "base excision",
      "double-strand break", "checkpoint"]),
    ("Cell cycle / mitotic division",
     "Cell cycle and mitosis",
     ["cell cycle", "mitotic", "chromosome", "cohesin", "condensin",
      "cytokinesis", "spindle", "centromere"]),
    ("Apoptosis / p53-mediated cell death",
     "Apoptosis and p53-mediated cell death",
     ["apoptosis", "p53", "programmed cell death", "cell death", "caspase"]),
    ("Interferon / innate immune signaling",
     "Immune and inflammatory signaling",
     ["interferon", "immune", "cytokine", "inflammasome", "antigen",
      "nf-kb", "innate", "complement"]),
    ("Signal transduction / growth factor signaling",
     "Signal transduction",
     ["signaling", "receptor", "kinase", "mapk", "pi3k", "wnt", "notch",
      "hedgehog", "egfr", "vegf", "gpcr"]),
    ("Transcription / translation (gene expression)",
     "Gene expression regulation",
     ["transcription", "translation", "rna polymerase", "mrna", "ribosome",
      "splicing", "chromatin", "histone"]),
    ("Extracellular matrix / adhesion / migration",
     "Extracellular matrix organization and cell adhesion",
     ["extracellular matrix", "collagen", "integrin", "adhesion",
      "migration", "motility"]),
    ("Other metabolism",
     "Cytochrome P450-mediated xenobiotic/drug and lipid metabolism",
     ["metabol", "glycolysis", "oxidative", "respiratory", "tca cycle"]),
]

def classify(term):
    t = term.lower()
    for label, meta, kws in FAMILIES:
        if any(k in t for k in kws):
            return label, meta
    return "Other", "Other"

top = tab.head(15).copy()
top["family"] = [classify(t)[0] for t in top["Term"]]
top["meta"] = [classify(t)[1] for t in top["Term"]]

primary_label, primary_meta = classify(tab["Term"].iloc[0])
# strengthen the call if the top-5 ranking is dominated by one meta-mechanism
top5_meta = [classify(t)[1] for t in tab.head(5)["Term"]]
meta_counts = pd.Series(top5_meta).value_counts()
if meta_counts.max() >= 3:
    primary_meta = meta_counts.idxmax()

sigf = sig.copy()
if len(sigf):
    sigf["meta"] = [classify(t)[1] for t in sigf["Term"]]
    fam_counts_sig = sigf["meta"].value_counts().to_dict()
else:
    fam_counts_sig = {}

top_list = []
for _, r in tab.head(10).iterrows():
    lab, met = classify(r["Term"])
    top_list.append({
        "term": r["Term"],
        "overlap": r["Overlap"],
        "p_value": float(r["P-value"]),
        "adjusted_p_value": float(r["Adjusted P-value"]),
        "odds_ratio": float(r["Odds Ratio"]),
        "combined_score": float(r["Combined Score"]),
        "genes": r["Genes"].split(";") if isinstance(r["Genes"], str) else [],
        "family": lab,
        "meta_mechanism": met,
    })

n_fdr_sig = int(len(sig))
mechanism = {
    "task": "ls07-combination-treatment-mechanism",
    "comparison": f"{TREAT} vs {CONTROL}",
    "primary_mechanism": primary_meta,
    "primary_mechanism_detail": primary_label,
    "mechanism_evidence": {
        "n_de_genes": int(len(de_genes)),
        "n_up": int(len(up)),
        "n_down": int(len(down)),
        "de_criteria": {"padj_cutoff": PADJ_CUTOFF,
                        "abs_log2fc_cutoff": L2FC_CUTOFF},
        "query_genes_in_background": int(len(query_in_bg)),
        "background_gene_count": int(len(background)),
        "library_term_count": int(len(terms)),
        "terms_tested_with_overlap": int(len(tab)),
        "n_pathways_fdr_lt_0.05": n_fdr_sig,
        "min_raw_p_value": float(tab["P-value"].min()),
        "min_adjusted_p_value": float(tab["Adjusted P-value"].min()),
        "significant_meta_family_counts": fam_counts_sig,
    },
    "top_pathways": top_list,
    "method": {
        "de": (f"pydeseq2 {pydeseq2.__version__} (DESeq2-style negative-"
               f"binomial GLM), contrast condition [{TREAT} vs {CONTROL}], "
               "genes with total counts >= 10 across the 6 samples tested, "
               "padj < 0.05 and |log2FC| >= 1 called DE"),
        "enrichment": ("gseapy 1.1.4 enrichr, local offline mode; one-sided "
                       "hypergeometric/Fisher exact test against the supplied "
                       "background universe; Enrichr odds ratio and combined "
                       "score; Benjamini-Hochberg FDR"),
        "gene_set_library": in_manifest["resource_name"] + " (frozen GMT)",
        "background_universe": ("Reactome_2022.background.txt "
                                f"({len(background)} genes)"),
    },
    "causation_note": ("Pathway enrichment identifies over-represented "
                       "annotated gene sets among differentially expressed "
                       "genes. It is associative statistical evidence about "
                       "transcriptional state and does not by itself "
                       "demonstrate a causal mechanism."),
    "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
}
json.dump(mechanism, open(os.path.join(BASE, "mechanism_call.json"), "w",
                          encoding="utf-8"), indent=2, ensure_ascii=False)

# ------------------------------------------------------ resource manifest
out_manifest = dict(in_manifest)
out_manifest["verification"] = {
    "verified_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "gmt_sha256_verified": True,
    "background_sha256_verified": True,
    "mapping_sha256_verified": True,
    "term_count_verified": True,
    "background_gene_count_verified": True,
    "gmt_gene_union_equals_background": True,
    "recomputed_sha256": sha,
    "software": {
        "python": platform.python_version(),
        "gseapy": gp.__version__,
        "pydeseq2": pydeseq2.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
    },
    "network_used": False,
    "policy_compliance": ("No pathway library or identifier mapping was "
                          "downloaded or substituted; enrichment ran offline "
                          "against the frozen local GMT and supplied "
                          "background universe."),
}
json.dump(out_manifest, open(os.path.join(BASE, "resource_manifest.json"), "w",
                             encoding="utf-8"), indent=2, ensure_ascii=False)

print("[done] wrote pathway_enrichment.csv, mechanism_call.json, "
      "resource_manifest.json, de_results.csv", flush=True)
print(f"[call] primary mechanism: {primary_meta}", flush=True)
print("[top terms]")
print(tab.head(15)[["Term", "Overlap", "P-value", "Adjusted P-value"]].to_string())
