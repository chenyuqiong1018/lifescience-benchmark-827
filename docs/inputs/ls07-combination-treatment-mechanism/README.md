# Input contract: ls07-combination-treatment-mechanism

## Files and roles

- `counts_raw_unfiltered.csv`, `sample_layout.csv`, `ensg_to_gene_name.tsv`: the same frozen BixBench `bix-43` expression inputs used by the DE task.
- `Reactome_2022.gmt`: byte-frozen Enrichr library named `Reactome_2022`.
- `Reactome_2022.background.txt`: sorted unique union of all gene symbols in the GMT; this is the explicit tested background for the local deterministic extension.
- `Reactome_2022.manifest.json`: retrieval, license, format, count, and integrity metadata.

## Provenance and redistribution boundary

The expression files came from the official BixBench v1.5 capsule retrieved on 2026-08-14. The named gene-set library was downloaded from the official Enrichr library endpoint on 2026-08-16. Reactome annotation files are distributed under CC0; Enrichr must be acknowledged using the citations recorded in the manifest. The repository does not assert a new license over the BixBench experimental inputs.

## Schema, units, and missing values

The GMT uses pathway label plus stable Reactome ID in column 1, an empty description in column 2, and HGNC-style gene symbols in columns 3 onward. The background contains one symbol per line. Enrichment probabilities are unitless; missing statistics must be empty/null.

## Integrity and readiness boundary

All agent-visible files are covered by `../SHA256SUMS.tsv`. Packaging this resource resolves the missing-library and missing-universe input defect, but it does **not** by itself accept the LS07-2 oracle. A full frozen DE result, enrichment reference table, and 3/3 positive/negative acceptance suite are still required before formal scoring.
