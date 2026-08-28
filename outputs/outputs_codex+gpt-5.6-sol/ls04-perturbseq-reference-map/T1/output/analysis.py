"""Multi-view, leakage-safe perturb-seq reference mapping.

The controlled skills guide three complementary anonymous signature views and
a prespecified bootstrap sensitivity analysis.  The unavailable scVI/scGPT
model runtimes are represented by transparent local proxies; no trained-model
execution is claimed.
"""

from __future__ import annotations

import csv
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse
from scipy.stats import rankdata


N_FEATURES = 2_000
LATENT_DIM = 30
N_BOOTSTRAPS = 200
RANDOM_SEED = 82711
TARGETS = ("PABPC1", "NUDT21", "LEO1")


def decode(values):
    return np.asarray(
        [value.decode() if isinstance(value, bytes) else str(value) for value in values]
    )


def load_anonymous_signatures(path: Path):
    """Load expression/guide groups while excluding target metadata."""
    with h5py.File(path, "r") as handle:
        shape = tuple(int(value) for value in handle["X"].attrs["shape"])
        raw = sparse.csc_matrix(
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

    if raw.data.dtype.kind not in "iu" or np.any(raw.data < 0):
        raise ValueError(f"{path.name} is not a raw non-negative integer matrix")
    if len(np.unique(genes)) != len(genes):
        raise ValueError(f"Non-unique gene IDs in {path.name}")
    if np.any(guide_codes < 0):
        raise ValueError(f"Missing guide assignments in {path.name}")

    matrix = raw.astype(np.float32)
    library_sizes = np.asarray(matrix.sum(axis=1)).ravel()
    scaling = np.divide(
        10_000.0,
        library_sizes,
        out=np.zeros_like(library_sizes, dtype=np.float32),
        where=library_sizes > 0,
    )
    matrix = matrix.multiply(scaling[:, None]).tocsr()
    matrix.data = np.log1p(matrix.data)

    group_counts = np.bincount(guide_codes, minlength=len(guide_names))
    group_matrix = sparse.csr_matrix(
        (
            np.ones(len(guide_codes), dtype=np.float32),
            (guide_codes, np.arange(len(guide_codes))),
        ),
        shape=(len(guide_names), len(guide_codes)),
    )
    guide_means = (group_matrix @ matrix).toarray() / group_counts[:, None]
    is_nt_group = np.asarray([name.startswith("NT-") for name in guide_names])
    is_nt_cell = is_nt_group[guide_codes]
    nt_mean = np.asarray(matrix[is_nt_cell].mean(axis=0)).ravel()
    perturbation_groups = np.flatnonzero(~is_nt_group)

    return {
        "genes": genes,
        "signatures": guide_means[perturbation_groups] - nt_mean,
        "group_indices": perturbation_groups,
        "guide_names": guide_names,
        "group_counts": group_counts[perturbation_groups],
    }


def align_genes(query, reference):
    query_lookup = {gene: index for index, gene in enumerate(query["genes"])}
    reference_lookup = {
        gene: index for index, gene in enumerate(reference["genes"])
    }
    common = sorted(query_lookup.keys() & reference_lookup.keys())
    if not common:
        raise ValueError("No shared genes")
    query_columns = [query_lookup[gene] for gene in common]
    reference_columns = [reference_lookup[gene] for gene in common]
    return (
        np.asarray(common),
        query["signatures"][:, query_columns],
        reference["signatures"][:, reference_columns],
    )


def normalize_rows(values):
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    return values / np.maximum(norms, 1e-12)


def component_scores(query, reference):
    """Return direct, rank-based, latent, and equal-weight ensemble scores."""
    direct = normalize_rows(reference) @ normalize_rows(query).T

    # scGPT-guided transparent proxy: a gene-order representation is invariant
    # to monotone scale shifts and emphasizes the relative perturbation program.
    query_ranks = rankdata(query, axis=1, method="average")
    reference_ranks = rankdata(reference, axis=1, method="average")
    query_ranks -= query_ranks.mean(axis=1, keepdims=True)
    reference_ranks -= reference_ranks.mean(axis=1, keepdims=True)
    ranked = normalize_rows(reference_ranks) @ normalize_rows(query_ranks).T

    # scVI-guided transparent proxy: low-dimensional coordinates of the pooled
    # signatures after each dataset was independently NT-centered.
    pooled = np.vstack((query, reference))
    gram = pooled @ pooled.T
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues, kind="stable")[::-1]
    eigenvalues = np.maximum(eigenvalues[order], 0.0)
    eigenvectors = eigenvectors[:, order]
    dimension = min(LATENT_DIM, len(eigenvalues))
    latent = eigenvectors[:, :dimension] * np.sqrt(eigenvalues[:dimension])[None, :]
    query_latent = normalize_rows(latent[: len(query)])
    reference_latent = normalize_rows(latent[len(query) :])
    latent_scores = reference_latent @ query_latent.T

    ensemble = (direct + ranked + latent_scores) / 3.0
    return direct, ranked, latent_scores, ensemble


def select_features(query, reference):
    pooled = np.vstack((query, reference))
    n_features = min(N_FEATURES, pooled.shape[1])
    selected = np.argsort(np.var(pooled, axis=0), kind="stable")[-n_features:]
    return selected, query[:, selected], reference[:, selected]


def bootstrap_all_reference_rows(query, reference):
    """Feature bootstrap before labels; returns a winner for every ref row."""
    rng = np.random.default_rng(RANDOM_SEED)
    n_features = query.shape[1]
    winners = np.empty((N_BOOTSTRAPS, len(reference)), dtype=np.int16)
    for iteration in range(N_BOOTSTRAPS):
        columns = rng.integers(0, n_features, size=n_features)
        *_, scores = component_scores(query[:, columns], reference[:, columns])
        winners[iteration] = np.argmax(scores, axis=1)
    return winners


def delayed_targets(path: Path, reference_group_indices):
    """First and only access to target metadata, after all blind analyses."""
    with h5py.File(path, "r") as handle:
        categories = decode(handle["obs/target_gene/categories"][:])
        target_codes = handle["obs/target_gene/codes"][:].astype(np.int64)
        guide_codes = handle["obs/guide/codes"][:].astype(np.int64)
    labels = []
    for group_index in reference_group_indices:
        codes = np.unique(target_codes[guide_codes == group_index])
        if len(codes) != 1 or codes[0] < 0:
            raise ValueError(f"Reference group {group_index} lacks one target")
        labels.append(categories[codes[0]])
    return np.asarray(labels)


def main():
    repo = Path(__file__).resolve().parents[5]
    input_dir = repo / "inputs" / "ls04-perturbseq-reference-map"
    query_path = input_dir / "perturb.seq.align.q1.query.h5ad"
    reference_path = input_dir / "perturb.seq.align.q1.ref.h5ad"

    query = load_anonymous_signatures(query_path)
    reference = load_anonymous_signatures(reference_path)
    common_genes, query_signatures, reference_signatures = align_genes(
        query, reference
    )
    selected, query_selected, reference_selected = select_features(
        query_signatures, reference_signatures
    )
    direct, ranked, latent, ensemble = component_scores(
        query_selected, reference_selected
    )
    bootstrap_winners = bootstrap_all_reference_rows(
        query_selected, reference_selected
    )

    # Leakage firewall: all base scores and all bootstrap winners already exist.
    reference_targets = delayed_targets(reference_path, reference["group_indices"])
    query_ids = query["guide_names"][query["group_indices"]]

    rows = []
    for target in TARGETS:
        matches = np.flatnonzero(reference_targets == target)
        if len(matches) != 1:
            raise ValueError(f"Expected one reference guide for {target}")
        reference_row = int(matches[0])
        order = np.argsort(ensemble[reference_row], kind="stable")[::-1]
        winner, runner_up = int(order[0]), int(order[1])
        confidence = np.mean(bootstrap_winners[:, reference_row] == winner)
        rows.append(
            {
                "target_gene": target,
                "query_guide_id": query_ids[winner],
                "score": f"{ensemble[reference_row, winner]:.6f}",
                "runner_up_score": f"{ensemble[reference_row, runner_up]:.6f}",
                "confidence": f"{confidence:.3f}",
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
    print(f"latent_dim={LATENT_DIM} feature_bootstraps={N_BOOTSTRAPS}")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
