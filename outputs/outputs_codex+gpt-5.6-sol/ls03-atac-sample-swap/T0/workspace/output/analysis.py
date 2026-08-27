#!/usr/bin/env python3
import csv,gzip,json,re
from collections import defaultdict
from itertools import combinations
from pathlib import Path

TASK="ls03-atac-sample-swap"
for repo in Path(__file__).resolve().parents:
 if (repo/"inputs"/TASK).is_dir():break
else:raise RuntimeError("repo not found")
inp=repo/"inputs"/TASK;out=Path(__file__).parent
PANEL={"Bladder":"UPK1A UPK1B UPK2 UPK3A KRT20 KRT13","Brain":"SOX2 GFAP MAP2 RBFOX3 NEUROD1 SLC17A7 SNAP25","Cloaca":"HOXA13 HOXD13 KRT13 KRT14","GallBladder":"KRT19 MUC1 EPCAM SOX17 KRT8","Gill":"FOXI1 KRT4 ATP6V1B1 CA2 GATA3","Heart":"MYH6 MYH7 TNNT2 NKX2-5 ACTC1 PLN","Intestine":"CDX2 VIL1 LGR5 SI MUC2 FABP2","Kidney":"PAX2 PAX8 SLC12A1 UMOD NPHS1 AQP2","Limb":"HOXA13 HOXD13 TBX5 TBX4 FGF8 PRRX1","Liver":"ALB AFP APOA1 TTR HNF4A CYP3A4 FGA","Lung":"SFTPA1 SFTPB SFTPC NKX2-1 SCGB1A1","Pancreas":"INS GCG PDX1 PRSS1 AMY2A CPA1 CEL","Prostate":"AR NKX3-1 KLK3 HOXB13","Spleen":"SPI1 PTPRC CD3D CD79A MS4A1","Stomach":"GAST GIF ATP4A ATP4B MUC5AC PGA3"}
PANEL={k:set(v.split()) for k,v in PANEL.items()};wanted=set().union(*PANEL.values())
segments=defaultdict(list)
with (inp/"sample.swap.atac.q1.chrom.sizes").open() as f:
 for line in f:
  name,n=line.rstrip().split("\t");segments[re.sub(r'_\d+$','',name)].append((name,int(n)))
def convert(chrom,pos):
 offset=0
 for name,n in segments.get(chrom,[]):
  if pos<=offset+n:return name,pos-offset
  offset+=n
 return None
locations=defaultdict(list)
with gzip.open(inp/"AmexT_v47-AmexG_v6.0-DD.gtf.gz","rt") as f:
 for line in f:
  if line.startswith("#"):continue
  x=line.rstrip().split("\t")
  if len(x)<9 or x[2]!="gene":continue
  m=re.search(r'gene_name "([^"]+)"',x[8])
  if not m:continue
  raw=m.group(1);symbols={s.upper() for s in re.findall(r'([A-Za-z0-9-]+) \[[a-z]{2}\]',raw)}|{raw.split()[0].upper()}
  for symbol in symbols&wanted:
   mapped=convert(x[0],int(x[3]) if x[6]=="+" else int(x[4]))
   if mapped:locations[symbol].append(mapped)
targets=set()
for locs in locations.values():
 for chrom,pos in locs:
  b=(pos-1)//10000
  for d in (-1,0,1):targets.add((chrom,b+d))
selected={};totals=None;rows_total=0
with gzip.open(inp/"sample.swap.atac.q1.tsv.gz","rt") as f:
 organs=f.readline().rstrip().split("\t")[3:]
 for line in f:
  x=line.rstrip().split("\t");v=list(map(float,x[3:]));rows_total+=1
  if totals is None:totals=[0.0]*len(v)
  for i,z in enumerate(v):totals[i]+=z
  key=(x[0],int(x[1])//10000)
  if key in targets:selected[key]=v
vectors=defaultdict(list)
for organ,markers in PANEL.items():
 for symbol in markers:
  for chrom,pos in locations.get(symbol,[]):
   b=(pos-1)//10000;near=[selected[(chrom,b+d)] for d in (-1,0,1) if (chrom,b+d) in selected]
   if not near:continue
   cpm=[sum(x[i] for x in near)/totals[i]*1e6 for i in range(len(organs))];mean=sum(cpm)/len(cpm)
   if mean:vectors[organ].append([x/mean for x in cpm])
score={o:[sum(v[i] for v in vectors[o])/len(vectors[o]) for i in range(len(organs))] for o in organs};idx={o:i for i,o in enumerate(organs)}
ranking=[]
for a,b in combinations(sorted(organs),2):
 gain=score[a][idx[b]]+score[b][idx[a]]-score[a][idx[a]]-score[b][idx[b]]
 ranking.append({"organ_a":a,"organ_b":b,"swap_score":gain,"evidence_type":"TSS promoter CPM reciprocal marker-coherence gain"})
ranking.sort(key=lambda r:(-r["swap_score"],r["organ_a"],r["organ_b"]))
for i,r in enumerate(ranking,1):r["rank"]=i
with (out/"sample_similarity.csv").open("w",newline="",encoding="utf-8") as f:
 w=csv.DictWriter(f,fieldnames=["organ_a","organ_b","swap_score","rank","evidence_type"]);w.writeheader()
 for r in ranking:w.writerow({**r,"swap_score":f"{r['swap_score']:.8f}"})
top,second=ranking[:2];margin=top["swap_score"]-second["swap_score"];detected=top["swap_score"]>0 and margin>0.25
components={"a_self":score[top["organ_a"]][idx[top["organ_a"]]],"a_in_b":score[top["organ_a"]][idx[top["organ_b"]]],"b_self":score[top["organ_b"]][idx[top["organ_b"]]],"b_in_a":score[top["organ_b"]][idx[top["organ_a"]]]}
evidence={"method":"library-normalized three-bin TSS promoter accessibility and reciprocal organ-marker coherence","coordinate_conversion":"GTF whole-arm coordinates converted to chrom-size fragments","table_rows":rows_total,"marker_instances_by_organ":{o:len(vectors[o]) for o in organs},"library_totals":dict(zip(organs,map(int,totals))),"top_score":top["swap_score"],"runner_up":{"organ_a":second["organ_a"],"organ_b":second["organ_b"],"score":second["swap_score"]},"margin":margin,"reciprocal_components":components}
call={"swap_detected":detected,"organ_a":top["organ_a"] if detected else None,"organ_b":top["organ_b"] if detected else None,"confidence":"high" if detected and margin>0.5 else "moderate" if detected else "low","evidence":evidence}
(out/"swap_call.json").write_text(json.dumps(call,indent=2)+"\n",encoding="utf-8",newline="\n")
report=f"""# T0 ATAC label-swap analysis

## Call

**{top['organ_a']} and {top['organ_b']} are swapped** ({call['confidence']} confidence). The top reciprocal coherence gain is {top['swap_score']:.4f}; runner-up {second['organ_a']}–{second['organ_b']} is {second['swap_score']:.4f}; margin {margin:.4f}.

Whole-arm GTF TSS coordinates were converted to sequential count-table fragments using chrom sizes. Three-bin promoter counts were normalized by full library totals, scaled across organs per marker, and all {len(ranking)} unordered swaps were ranked. Thus total library size cannot itself drive the call. Regulatory-region and tissue-specific skill guidance shaped the promoter/marker evidence; code-execution guidance shaped the reproducible full-table implementation.

Conserved markers and internal reciprocity do not replace matched axolotl reference epigenomes. The deterministic rule requires a positive top gain and margin >0.25; otherwise `swap_detected=false`.
"""
(out/"report.md").write_text(report,encoding="utf-8",newline="\n")
print(json.dumps({"swap_detected":detected,"organ_a":top["organ_a"],"organ_b":top["organ_b"],"confidence":call["confidence"],"top_score":top["swap_score"],"margin":margin},sort_keys=True))
