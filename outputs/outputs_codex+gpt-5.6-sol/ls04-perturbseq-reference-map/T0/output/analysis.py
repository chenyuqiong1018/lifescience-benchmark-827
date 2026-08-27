"""Skill-guided latent reference mapping for perturb-seq guide groups.

scvi-tools itself is unavailable in the execution environment.  Following the
installed skill's guidance, this fallback preserves raw integer counts, removes
the dataset shift using within-dataset non-targeting controls, and transfers
labels in a low-dimensional latent representation.  It never reads reference
target metadata until all anonymous reference-by-query scores are computed.
"""

from __future__ import annotations

import csv
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse


N_FEATURES = 2_000
PRIMARY_LATENT_DIM = 30
LATENT_DIMENSIONS = (10, 20, 30, 40, 50, 60)
TARGETS = ("PABPC1", "NUDT21", "LEO1")


def decode(values):
    return np.asarray(
        [value.decode() if isinstance(value, bytes) else str(value) for value in values]
    )


def blind_guide_signatures(path: Path):
    """Create within-dataset perturbation signatures without target labels."""
    with h5py.File(path, "r") as handle:
        shape = tuple(int(value) for value in handle["X"].attrs["shape"])
        raw_counts = sparse.csc_matrix(
            (
                handle["X/data"][:],
                handle["X/indices"][:],
                handle["X/indptr"][:],
            ),
            shape=shape,
        ).tocsr()
        genes = decode(handle["var/_index"][:])
        guide_names = decode(handle["obs/guide/categories"][:])
        guide_codes = handle["obs/guide/codes"][:].astype(np.int64)

    if raw_counts.data.dtype.kind not in "iu" or np.any(raw_counts.data < 0):
        raise ValueError(f"Expected raw non-negative integer counts in {path.name}")
    if len(np.unique(genes)) != len(genes):
        raise ValueError(f"Gene IDs are not unique in {path.name}")

    # Raw counts are retained above. Log-CP10K is used only for this documented
    # CPU latent-space fallback, not passed to a negative-binomial scVI model.
    matrix = raw_counts.astype(np.float32)
    library_size = np.asarray(matrix.sum(axis=1)).ravel()
    scaling = np.divide(
        10_000.0,
        library_size,
        out=np.zeros_like(library_size, dtype=np.float32),
        where=library_size > 0,
    )
    matrix = matrix.multiply(scaling[:, None]).tocsr()
    matrix.data = np.log1p(matrix.data)

    counts_per_group = np.bincount(guide_codes, minlength=len(guide_names))
    membership = sparse.csr_matrix(
        (
            np.ones(len(guide_codes), dtype=np.float32),
            (guide_codes, np.arange(len(guide_codes))),
        ),
        shape=(len(guide_names), len(guide_codes)),
    )
    means = (membership @ matrix).toarray() / counts_per_group[:, None]
    is_nt_group = np.asarray([name.startswith("NT-") for name in guide_names])
    is_nt_cell = is_nt_group[guide_codes]
    nt_mean = np.asarray(matrix[is_nt_cell].mean(axis=0)).ravel()
    perturbation_groups = np.flatnonzero(~is_nt_group)
    signatures = means[perturbation_groups] - nt_mean

    return {
        "genes": genes,
        "signatures": signatures,
        "group_indices": perturbation_groups,
        "guide_names": guide_names,
        "cell_counts": counts_per_group[perturbation_groups],
    }


def align_shared_genes(query, reference):
    query_index = {gene: i for i, gene in enumerate(query["genes"])}
    reference_index = {gene: i for i, gene in enumerate(reference["genes"])}
    common = sorted(query_index.keys() & reference_index.keys())
    if not common:
        raise ValueError("No shared genes")
    q_columns = [query_index[gene] for gene in common]
    r_columns = [reference_index[gene] for gene in common]
    return (
        np.asarray(common),
        query["signatures"][:, q_columns],
        reference["signatures"][:, r_columns],
    )


def normalize_rows(values):
    return values / np.maximum(np.linalg.norm(values, axis=1, keepdims=True), 1e-12)


def latent_score_matrices(query_signatures, reference_signatures):
    """Compute anonymous score matrices across fixed latent dimensions."""
    pooled = np.vstack((query_signatures, reference_signatures))
    n_features = min(N_FEATURES, pooled.shape[1])
    selected = np.argsort(np.var(pooled, axis=0), kind="stable")[-n_features:]
    selected_signatures = pooled[:, selected]

    # For n_guides << n_genes, eigendecomposition of the guide Gram matrix is
    # the exact compact route to the left singular-vector latent coordinates.
    gram = selected_signatures @ selected_signatures.T
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues, kind="stable")[::-1]
    eigenvalues = np.maximum(eigenvalues[order], 0.0)
    eigenvectors = eigenvectors[:, order]

    n_query = len(query_signatures)
    score_matrices = {}
    for dimension in LATENT_DIMENSIONS:
        dimension = min(dimension, len(eigenvalues))
        latent = eigenvectors[:, :dimension] * np.sqrt(eigenvalues[:dimension])[None, :]
        query_latent = normalize_rows(latent[:n_query])
        reference_latent = normalize_rows(latent[n_query:])
        score_matrices[dimension] = reference_latent @ query_latent.T
    return score_matrices, selected


def delayed_reference_targets(path: Path, group_indices):
    """Read target labels only after every latent score matrix is complete."""
    with h5py.File(path, "r") as handle:
        target_categories = decode(handle["obs/target_gene/categories"][:])
        target_codes = handle["obs/target_gene/codes"][:].astype(np.int64)
        guide_codes = handle["obs/guide/codes"][:].astype(np.int64)
    targets = []
    for group_index in group_indices:
        codes = np.unique(target_codes[guide_codes == group_index])
        if len(codes) != 1 or codes[0] < 0:
            raise ValueError(f"Reference group {group_index} has ambiguous target metadata")
        targets.append(target_categories[codes[0]])
    return np.asarray(targets)


def main():
    repo = Path(__file__).resolve().parents[5]
    input_dir = repo / "inputs" / "ls04-perturbseq-reference-map"
    query_path = input_dir / "perturb.seq.align.q1.query.h5ad"
    reference_path = input_dir / "perturb.seq.align.q1.ref.h5ad"

    query = blind_guide_signatures(query_path)
    reference = blind_guide_signatures(reference_path)
    common_genes, query_signatures, reference_signatures = align_shared_genes(
        query, reference
    )
    score_matrices, selected = latent_score_matrices(
        query_signatures, reference_signatures
    )

    # Explicit leakage firewall: this is the first target-metadata access.
    reference_targets = delayed_reference_targets(
        reference_path, reference["group_indices"]
    )
    query_ids = query["guide_names"][query["group_indices"]]
    primary_scores = score_matrices[PRIMARY_LATENT_DIM]

    rows = []
    for target in TARGETS:
        matches = np.flatnonzero(reference_targets == target)
        if len(matches) != 1:
            raise ValueError(f"Expected one reference guide for {target}")
        reference_row = int(matches[0])
        order = np.argsort(primary_scores[reference_row], kind="stable")[::-1]
        winner, runner_up = int(order[0]), int(order[1])
        stable = sum(
            int(np.argmax(scores[reference_row]) == winner)
            for scores in score_matrices.values()
        )
        rows.append(
            {
                "target_gene": target,
                "query_guide_id": query_ids[winner],
                "score": f"{primary_scores[reference_row, winner]:.6f}",
                "runner_up_score": f"{primary_scores[reference_row, runner_up]:.6f}",
                "confidence": f"{stable / len(score_matrices):.3f}",
            }
        )

    output_path = Path(__file__).with_name("guide_mapping.csv")
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "target_gene",
                "query_guide_id",
                "score",
                "runner_up_score",
                "confidence",
            ),
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"shared_genes={len(common_genes)} selected_features={len(selected)}")
    print(f"primary_latent_dim={PRIMARY_LATENT_DIM} dimensions={LATENT_DIMENSIONS}")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
