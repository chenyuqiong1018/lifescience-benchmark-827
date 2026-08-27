#!/usr/bin/env python3
import gzip,json,hashlib
from pathlib import Path

TASK="ls02-infer-genome-build"
for root in Path(__file__).resolve().parents:
    if (root/"inputs"/TASK).is_dir(): break
else: raise RuntimeError("repo not found")
inp=root/"inputs"/TASK; out=Path(__file__).parent; refs=inp/"references"
vcf=inp/"vcf.infer.build.q1.vcf.gz"
paths={b:refs/f"{b}_chr20.fa.gz" for b in ("hg18","hg19","hg38")}
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for x in iter(lambda:f.read(1048576),b""): h.update(x)
    return h.hexdigest()
def fasta(p):
    with gzip.open(p,"rt") as f: return f.readline().strip(),"".join(x.strip().upper() for x in f)
variants=[]; chroms=set()
with gzip.open(vcf,"rt") as f:
    for line in f:
        if line.startswith("#"): continue
        x=line.rstrip().split("\t"); chroms.add(x[0]); variants.append((int(x[1]),x[3].upper()))
if chroms not in ({"20"},{"chr20"}): raise ValueError(chroms)
checks={}
for build,p in paths.items():
    header,seq=fasta(p); match=oor=0; examples=[]
    for pos,ref in variants:
        observed=seq[pos-1:pos-1+len(ref)] if 0<pos<=len(seq) else ""
        if observed==ref: match+=1
        else:
            if len(observed)!=len(ref): oor+=1
            if len(examples)<5: examples.append({"pos":pos,"ref":ref,"observed":observed})
    checks[build]={"matches":match,"mismatches":len(variants)-match,"out_of_range":oor,"fraction":match/len(variants),"length_bp":len(seq),"header":header,"sha256":sha(p),"mismatch_examples":examples}
rank=sorted(checks,key=lambda b:checks[b]["matches"],reverse=True); winner=rank[0]; best=checks[winner]
confidence="high" if best["matches"]==len(variants) and checks[rank[1]]["matches"]<len(variants) else "medium" if best["matches"]>checks[rank[1]]["matches"] else "low"
evidence={"method":"exhaustive FASTA[POS-1:POS-1+len(REF)] comparison","chromosome_labels":sorted(chroms),"candidate_checks":checks,"winner_margin_matches":best["matches"]-checks[rank[1]]["matches"],"t2t":{"status":"not_evaluated_reference_absent","reason":"No T2T/hs1 reference was supplied; absence is not a mismatch."},"vcf_sha256":sha(vcf)}
call={"build":winner,"confidence":confidence,"n_variants_checked":len(variants),"n_ref_matches":best["matches"],"n_ref_mismatches":best["mismatches"],"evidence":evidence}
(out/"build_call.json").write_text(json.dumps(call,indent=2)+"\n",encoding="utf-8",newline="\n")
table="\n".join(f"| {b} | {d['matches']:,} | {d['mismatches']:,} | {d['fraction']:.6f} |" for b,d in checks.items())
report=f"""# T0 genome-build inference

## Call

The VCF uses **{winner}** coordinates with **{confidence} confidence**: {best['matches']:,}/{len(variants):,} REF alleles match and {best['mismatches']:,} mismatch.

| Candidate | REF matches | REF mismatches | Fraction |
|---|---:|---:|---:|
{table}

Every VCF REF was checked at its declared 1-based coordinate; chromosome naming was QC only. T2T is explicitly untested because no T2T reference was supplied. The selected chromosome/UCSC skills informed assembly-aware coordinate checks, while code-execution guidance informed the auditable exhaustive script.
"""
(out/"report.md").write_text(report,encoding="utf-8",newline="\n")
print(json.dumps({"build":winner,"confidence":confidence,"checked":len(variants),"matches":best["matches"],"mismatches":best["mismatches"]},sort_keys=True))
