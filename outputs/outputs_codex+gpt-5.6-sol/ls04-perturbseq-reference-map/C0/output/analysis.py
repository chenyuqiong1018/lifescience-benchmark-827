"""Blind perturbation-reference mapping for task ls04-perturbseq-reference-map.

The score matrix is computed without reading ``obs/target_gene``.  Target labels
are opened only after feature selection and every reference-by-query score has
been frozen in memory.
"""

from __future__ import annotations

import csv
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse


N_FEATURES = 2_000
N_BOOTSTRAPS = 500
RANDOM_SEED = 827
TARGETS = ("PABPC1", "NUDT21", "LEO1")


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [value.decode() if isinstance(value, bytes) else str(value) for value in values]
    )


def load_blind_signatures(path: Path):
    """Load counts and guide groups, deliberately excluding target metadata."""
    with h5py.File(path, "r") as handle:
        x_group = handle["X"]
        shape = tuple(int(value) for value in x_group.attrs["shape"])
        matrix = sparse.csc_matrix(
            (
                handle["X/data"][:],
                handle["X/indices"][:],
                handle["X/indptr"][:],
            ),
            shape=shape,
        ).tocsr().astype(np.float32)
        genes = decode(handle["var/_index"][:])
        guide_categories = decode(handle["obs/guide/categories"][:])
        guide_codes = handle["obs/guide/codes"][:].astype(np.int64)

    if len(np.unique(genes)) != len(genes):
        raise ValueError(f"Non-unique gene IDs in {path.name}")
    if np.any(guide_codes < 0):
        raise ValueError(f"Missing guide codes in {path.name}")

    totals = np.asarray(matrix.sum(axis=1)).ravel()
    scales = np.divide(
        10_000.0,
        totals,
        out=np.zeros_like(totals, dtype=np.float32),
        where=totals > 0,
    )
    matrix = matrix.multiply(scales[:, None]).tocsr()
    matrix.data = np.log1p(matrix.data)

    counts = np.bincount(guide_codes, minlength=len(guide_categories))
    membership = sparse.csr_matrix(
        (
            np.ones(len(guide_codes), dtype=np.float32),
            (guide_codes, np.arange(len(guide_codes))),
        ),
        shape=(len(guide_categories), len(guide_codes)),
    )
    means = (membership @ matrix).toarray()
    means /= counts[:, None]

    nt_group = np.asarray([name.startswith("NT-") for name in guide_categories])
    nt_cells = nt_group[guide_codes]
    if not np.any(nt_cells):
        raise ValueError(f"No non-targeting cells in {path.name}")
    baseline = np.asarray(matrix[nt_cells].mean(axis=0)).ravel()
    signatures = means - baseline
    perturbation_groups = np.flatnonzero(~nt_group)

    return {
        "genes": genes,
        "signatures": signatures[perturbation_groups],
        "group_indices": perturbation_groups,
        "guide_categories": guide_categories,
        "group_counts": counts[perturbation_groups],
    }


def align_genes(query, reference):
    query_lookup = {gene: index for index, gene in enumerate(query["genes"])}
    reference_lookup = {gene: index for index, gene in enumerate(reference["genes"])}
    common = sorted(query_lookup.keys() & reference_lookup.keys())
    if not common:
        raise ValueError("Query and reference share no genes")
    query_columns = np.asarray([query_lookup[gene] for gene in common])
    reference_columns = np.asarray([reference_lookup[gene] for gene in common])
    return (
        np.asarray(common),
        query["signatures"][:, query_columns],
        reference["signatures"][:, reference_columns],
    )


def row_normalize(values: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    return values / np.maximum(norms, 1e-12)


def compute_blind_scores(query_signatures, reference_signatures):
    # Variance uses all anonymous guide signatures and no label information.
    pooled = np.vstack((query_signatures, reference_signatures))
    variances = np.var(pooled, axis=0)
    n_features = min(N_FEATURES, pooled.shape[1])
    selected = np.argsort(variances, kind="stable")[-n_features:]
    query_selected = query_signatures[:, selected]
    reference_selected = reference_signatures[:, selected]
    scores = row_normalize(reference_selected) @ row_normalize(query_selected).T
    return scores, selected, query_selected, reference_selected


def read_reference_targets(path: Path, reference_group_indices: np.ndarray):
    """Read labels only after the complete blind score matrix exists."""
    with h5py.File(path, "r") as handle:
        target_categories = decode(handle["obs/target_gene/categories"][:])
        target_codes = handle["obs/target_gene/codes"][:].astype(np.int64)
        guide_codes = handle["obs/guide/codes"][:].astype(np.int64)

    labels = []
    for group_index in reference_group_indices:
        cell_target_codes = np.unique(target_codes[guide_codes == group_index])
        if len(cell_target_codes) != 1 or cell_target_codes[0] < 0:
            raise ValueError(f"Reference guide group {group_index} lacks one target label")
        labels.append(target_categories[cell_target_codes[0]])
    return np.asarray(labels)


def bootstrap_support(
    query_selected: np.ndarray,
    reference_selected: np.ndarray,
    reference_row: int,
    winning_query: int,
) -> float:
    rng = np.random.default_rng(RANDOM_SEED + reference_row)
    wins = 0
    width = query_selected.shape[1]
    for _ in range(N_BOOTSTRAPS):
        columns = rng.integers(0, width, size=width)
        query_boot = row_normalize(query_selected[:, columns])
        reference_boot = row_normalize(reference_selected[reference_row : reference_row + 1, columns])
        if int(np.argmax(reference_boot @ query_boot.T)) == winning_query:
            wins += 1
    return wins / N_BOOTSTRAPS


def main() -> None:
    repo = Path(__file__).resolve().parents[5]
    input_dir = repo / "inputs" / "ls04-perturbseq-reference-map"
    query_path = input_dir / "perturb.seq.align.q1.query.h5ad"
    reference_path = input_dir / "perturb.seq.align.q1.ref.h5ad"

    query = load_blind_signatures(query_path)
    reference = load_blind_signatures(reference_path)
    common_genes, query_signatures, reference_signatures = align_genes(query, reference)
    scores, selected, query_selected, reference_selected = compute_blind_scores(
        query_signatures, reference_signatures
    )

    # Leakage firewall: target metadata is not accessed until all scores and
    # unsupervised selected-feature indices have been computed.
    reference_targets = read_reference_targets(
        reference_path, reference["group_indices"]
    )
    query_ids = query["guide_categories"][query["group_indices"]]

    rows = []
    for target in TARGETS:
        matches = np.flatnonzero(reference_targets == target)
        if len(matches) != 1:
            raise ValueError(f"Expected one reference group for {target}, found {len(matches)}")
        reference_row = int(matches[0])
        order = np.argsort(scores[reference_row], kind="stable")[::-1]
        winner, runner_up = int(order[0]), int(order[1])
        rows.append(
            {
                "target_gene": target,
                "query_guide_id": query_ids[winner],
                "score": f"{scores[reference_row, winner]:.6f}",
                "runner_up_score": f"{scores[reference_row, runner_up]:.6f}",
                "confidence": f"{bootstrap_support(query_selected, reference_selected, reference_row, winner):.3f}",
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
    print(f"query_groups={len(query_ids)} reference_groups={len(reference_targets)}")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
