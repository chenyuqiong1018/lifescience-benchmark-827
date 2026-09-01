from pathlib import Path
import shutil
import sys

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


def main(input_dir: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    ann = pd.read_csv(input_dir / "ensembl112_gene_coordinates.tsv", sep="\t")
    rna = pd.read_csv(input_dir / "multiome.match.atac.rna.q1.rna.tsv.gz", sep="\t")
    atac = pd.read_csv(input_dir / "multiome.match.atac.rna.q1.atac.tsv.gz", sep="\t")

    gene_col = "Gene name"
    chrom_col = "Chromosome/scaffold name"
    start_col = "Gene start (bp)"
    end_col = "Gene end (bp)"
    strand_col = "Strand"
    valid_chr = {str(i) for i in range(1, 23)} | {"X", "Y"}
    ann = ann[ann[chrom_col].astype(str).isin(valid_chr) & ann[gene_col].notna()].copy()
    counts = ann[gene_col].value_counts()
    ann = ann[ann[gene_col].map(counts).eq(1) & ann[gene_col].isin(rna.iloc[:, 0])].copy()
    tss = np.where(ann[strand_col].astype(float) >= 0, ann[start_col], ann[end_col]).astype(np.int64)
    bin_start = ((tss - 1) // 10000) * 10000
    ann["peak"] = "chr" + ann[chrom_col].astype(str) + "_" + bin_start.astype(str) + "_" + (bin_start + 10000).astype(str)

    atac = atac.set_index(atac.columns[0])
    rna = rna.set_index(rna.columns[0])
    ann = ann[ann["peak"].isin(atac.index)].drop_duplicates(gene_col)
    genes = ann[gene_col].tolist()
    peaks = ann.set_index(gene_col).loc[genes, "peak"]
    r = np.log1p(rna.loc[genes].astype(float))
    a = np.log1p(atac.loc[peaks].astype(float))
    a.index = genes
    keep = r.var(axis=1, ddof=1).nlargest(min(2000, len(r))).index
    r = r.loc[keep]
    a = a.loc[keep]

    scores = np.empty((r.shape[1], a.shape[1]))
    for i, rc in enumerate(r.columns):
        for j, ac in enumerate(a.columns):
            x, y = r[rc].to_numpy(), a[ac].to_numpy()
            mask = np.isfinite(x) & np.isfinite(y)
            scores[i, j] = np.corrcoef(x[mask], y[mask])[0, 1]
    score_df = pd.DataFrame(scores, index=r.columns, columns=a.columns)
    score_df.index.name = "rna_population"
    score_df.to_csv(output_dir / "score_matrix.csv")

    rows, cols = linear_sum_assignment(-np.nan_to_num(scores, nan=-1e9))
    mapping = []
    for i, j in zip(rows, cols):
        row_scores = scores[i]
        runner = np.sort(row_scores[np.isfinite(row_scores)])[-2]
        mapping.append((str(r.columns[i]), str(a.columns[j]), scores[i, j], runner))
    mapping_df = pd.DataFrame(mapping, columns=["rna_population", "atac_column", "match_score", "runner_up_score"])
    mapping_df.to_csv(output_dir / "column_mapping.csv", index=False)

    shutil.copy2(Path(__file__), output_dir / "analysis.py")
    total = mapping_df["match_score"].sum()
    report = f"""# ATAC–RNA population column matching

The analysis followed the frozen rule: unique Ensembl release 112 gene symbols on chromosomes 1–22/X/Y were mapped by strand-aware TSS to 10 kb ATAC bins. RNA TPM and ATAC values were transformed with `log1p`; the 2,000 mapped genes with highest RNA variance were used for Pearson correlation. A Hungarian assignment maximized the total correlation while enforcing a bijection.

- Unique mapped genes before variance selection: {len(genes)}
- Genes used for scoring: {len(keep)}
- Total assigned correlation: {total:.6f}

The shared biological signal is coordinated population-specific gene expression and chromatin accessibility at the transcription start sites of the same variable genes. `runner_up_score` is the second-highest correlation within each RNA-population row and indicates assignment ambiguity.

## Final mapping

{mapping_df.to_markdown(index=False)}
"""
    (output_dir / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
