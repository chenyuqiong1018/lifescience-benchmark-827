#!/usr/bin/env python3
"""Detect a unique organ-label swap from normalized promoter-marker coherence."""
import csv,gzip,hashlib,json,re
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path

TASK="ls03-atac-sample-swap"
for repo in Path(__file__).resolve().parents:
    if (repo/"inputs"/TASK).is_dir(): break
else: raise RuntimeError("repository not found")
inp=repo/"inputs"/TASK; out=Path(__file__).resolve().parent
gtf=inp/"AmexT_v47-AmexG_v6.0-DD.gtf.gz"; table=inp/"sample.swap.atac.q1.tsv.gz"; sizes=inp/"sample.swap.atac.q1.chrom.sizes"
MARKERS={
"Bladder":"UPK1A UPK1B UPK2 UPK3A KRT20 KRT13","Brain":"SOX2 GFAP MAP2 RBFOX3 NEUROD1 SLC17A7 SNAP25","Cloaca":"HOXA13 HOXD13 KRT13 KRT14","GallBladder":"KRT19 MUC1 EPCAM SOX17 KRT8","Gill":"FOXI1 KRT4 ATP6V1B1 CA2 GATA3","Heart":"MYH6 MYH7 TNNT2 NKX2-5 ACTC1 PLN","Intestine":"CDX2 VIL1 LGR5 SI MUC2 FABP2","Kidney":"PAX2 PAX8 SLC12A1 UMOD NPHS1 AQP2","Limb":"HOXA13 HOXD13 TBX5 TBX4 FGF8 PRRX1","Liver":"ALB AFP APOA1 TTR HNF4A CYP3A4 FGA","Lung":"SFTPA1 SFTPB SFTPC NKX2-1 SCGB1A1","Pancreas":"INS GCG PDX1 PRSS1 AMY2A CPA1 CEL","Prostate":"AR NKX3-1 KLK3 HOXB13","Spleen":"SPI1 PTPRC CD3D CD79A MS4A1","Stomach":"GAST GIF ATP4A ATP4B MUC5AC PGA3"}
MARKERS={k:set(v.split()) for k,v in MARKERS.items()}; wanted=set().union(*MARKERS.values())
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for x in iter(lambda:f.read(1<<20),b""):h.update(x)
 return h.hexdigest()

parts=defaultdict(list)
with sizes.open() as f:
 for line in f:
  name,n=line.rstrip().split("\t");parts[re.sub(r'_\d+$','',name)].append((name,int(n)))
def remap(chrom,pos1):
 offset=0
 for name,n in parts.get(chrom,[]):
  if pos1<=offset+n:return name,pos1-offset
  offset+=n
 return None

genes=defaultdict(list); gene_records=0
with gzip.open(gtf,"rt") as f:
 for line in f:
  if line.startswith("#"):continue
  x=line.rstrip().split("\t")
  if len(x)<9 or x[2]!="gene":continue
  gene_records+=1;m=re.search(r'gene_name "([^"]+)"',x[8])
  if not m:continue
  raw=m.group(1); symbols={s.upper() for s in re.findall(r'([A-Za-z0-9-]+) \[[a-z]{2}\]',raw)}|{raw.split()[0].upper()}
  for symbol in symbols&wanted:
   pos=int(x[3]) if x[6]=="+" else int(x[4]);mapped=remap(x[0],pos)
   if mapped:genes[symbol].append(mapped)

targets=set()
for locations in genes.values():
 for chrom,pos in locations:
  b=(pos-1)//10000
  for delta in (-1,0,1):targets.add((chrom,b+delta))
values={};totals=None;row_count=0
with gzip.open(table,"rt") as f:
 header=f.readline().rstrip().split("\t");organs=header[3:]
 for line in f:
  x=line.rstrip().split("\t");vec=list(map(float,x[3:]));row_count+=1
  if totals is None:totals=[0.0]*len(vec)
  totals=[a+b for a,b in zip(totals,vec)];key=(x[0],int(x[1])//10000)
  if key in targets:values[key]=vec
if set(organs)!=set(MARKERS):raise ValueError("organ columns differ from marker panel")

gene_vectors=defaultdict(list)
for organ,markers in MARKERS.items():
 for symbol in sorted(markers):
  for chrom,pos in genes.get(symbol,[]):
   b=(pos-1)//10000;rows=[values[(chrom,b+d)] for d in (-1,0,1) if (chrom,b+d) in values]
   if not rows:continue
   cpm=[sum(r[i] for r in rows)/totals[i]*1e6 for i in range(len(organs))];mean=sum(cpm)/len(cpm)
   if mean>0:gene_vectors[organ].append([v/mean for v in cpm])
scores={organ:[sum(v[i] for v in gene_vectors[organ])/len(gene_vectors[organ]) for i in range(len(organs))] for organ in organs}
index={o:i for i,o in enumerate(organs)};pairs=[]
for a,b in combinations(sorted(organs),2):
 gain=scores[a][index[b]]+scores[b][index[a]]-scores[a][index[a]]-scores[b][index[b]]
 pairs.append({"organ_a":a,"organ_b":b,"swap_score":gain,"evidence_type":"library-normalized promoter-marker coherence gain"})
pairs.sort(key=lambda x:(-x["swap_score"],x["organ_a"],x["organ_b"]))
for rank,row in enumerate(pairs,1):row["rank"]=rank
with (out/"sample_similarity.csv").open("w",newline="",encoding="utf-8") as f:
 w=csv.DictWriter(f,fieldnames=["organ_a","organ_b","swap_score","rank","evidence_type"]);w.writeheader()
 for r in pairs:w.writerow({**r,"swap_score":f"{r['swap_score']:.8f}"})
top,second=pairs[:2]; detected=top["swap_score"]>0 and top["swap_score"]-second["swap_score"]>0.25
evidence={"method":"three-bin TSS promoter accessibility, CPM library normalization, per-marker across-organ scaling, reciprocal diagonal-coherence gain","coordinate_remapping":"whole chromosome-arm GTF coordinates converted to sequential chrom-size fragments","table_rows":row_count,"gtf_gene_records":gene_records,"marker_instances_by_organ":{o:len(gene_vectors[o]) for o in organs},"library_totals":dict(zip(organs,map(int,totals))),"current_diagonal_scores":{o:scores[o][index[o]] for o in organs},"top_pair_score":top["swap_score"],"second_pair":{"organ_a":second["organ_a"],"organ_b":second["organ_b"],"score":second["swap_score"]},"top_margin":top["swap_score"]-second["swap_score"],"reciprocal_components":{"organ_a_self":scores[top["organ_a"]][index[top["organ_a"]]],"organ_a_in_b":scores[top["organ_a"]][index[top["organ_b"]]],"organ_b_self":scores[top["organ_b"]][index[top["organ_b"]]],"organ_b_in_a":scores[top["organ_b"]][index[top["organ_a"]]]},"input_sha256":{"gtf":sha(gtf),"count_table":sha(table),"chrom_sizes":sha(sizes)}}
call={"swap_detected":detected,"organ_a":top["organ_a"] if detected else None,"organ_b":top["organ_b"] if detected else None,"confidence":"high" if detected and evidence["top_margin"]>0.5 else "moderate" if detected else "low","evidence":evidence}
(out/"swap_call.json").write_text(json.dumps(call,indent=2)+"\n",encoding="utf-8",newline="\n")
report=f"""# ATAC-seq organ-label swap analysis

## Call

**Swap detected: {str(detected).lower()} — {top['organ_a']} and {top['organ_b']}** ({call['confidence']} confidence).

The unique top reciprocal marker-coherence gain is {top['swap_score']:.4f}; the runner-up is {second['organ_a']}–{second['organ_b']} at {second['swap_score']:.4f}, a margin of {evidence['top_margin']:.4f}. {top['organ_a']} markers rise from {evidence['reciprocal_components']['organ_a_self']:.3f} in their labeled column to {evidence['reciprocal_components']['organ_a_in_b']:.3f} in {top['organ_b']}; {top['organ_b']} markers rise from {evidence['reciprocal_components']['organ_b_self']:.3f} to {evidence['reciprocal_components']['organ_b_in_a']:.3f} in {top['organ_a']}.

## Method and safeguards

Official GTF gene TSS coordinates were converted from whole chromosome arms to sequential `_1/_2/...` ATAC fragments using the supplied chromosome sizes. Counts in the TSS bin and its two neighbors were divided by full-table library totals, so the decision is not based on library size. Conserved organ-marker promoter profiles were scaled across samples, and every one of the {len(pairs)} unordered pairs was scored by the improvement in reciprocal diagonal coherence. `sample_similarity.csv` contains all pairs, ordered by decreasing finite score.

## Limits

Marker sets are curated conserved vertebrate markers and are not a substitute for matched axolotl reference epigenomes. Confidence therefore reflects the unique internal reciprocal signal and its margin, not independent sample provenance. The decision rule requires a positive top gain and margin >0.25; otherwise the script returns `swap_detected=false`.
"""
(out/"report.md").write_text(report,encoding="utf-8",newline="\n")
print(json.dumps({"swap_detected":detected,"organ_a":top["organ_a"],"organ_b":top["organ_b"],"confidence":call["confidence"],"top_score":top["swap_score"],"margin":evidence["top_margin"]},sort_keys=True))
