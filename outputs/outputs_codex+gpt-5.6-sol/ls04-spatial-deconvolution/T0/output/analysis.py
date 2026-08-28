"""scvi-tools-guided count deconvolution for Spot_710-1.

The installed skill routes spatial composition tasks to a raw-count reference
signature model (cell2location/DestVI/Tangram family).  Those packages are not
available here, so this script implements a transparent marker-aware
multinomial abundance model rather than claiming to run them.
"""

from __future__ import annotations

import csv
import gzip
import io
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread


TARGET_SPOT = "Spot_710-1"
PRIMARY_MARKERS_PER_TYPE = 100
MARKER_SENSITIVITY = (50, 75, 100, 150, 200)
REPORT_THRESHOLD = 0.02
N_BOOTSTRAPS = 200
RANDOM_SEED = 121827


def load_input(path: Path):
    with tarfile.open(path, "r:gz") as archive:
        counts = pd.read_csv(archive.extractfile("spatial_q_sc_counts.csv"), index_col=0)
        metadata = pd.read_csv(
            archive.extractfile("spatial_q_sc_metadata.csv"), index_col=0
        )
        with gzip.GzipFile(fileobj=archive.extractfile("barcodes.tsv.gz")) as handle:
            barcodes = [line.decode().strip() for line in handle]
        with gzip.GzipFile(fileobj=archive.extractfile("features.tsv.gz")) as handle:
            features = [line.decode().rstrip("\n").split("\t") for line in handle]
        with gzip.GzipFile(fileobj=archive.extractfile("matrix.mtx.gz")) as handle:
            spatial = mmread(io.BytesIO(handle.read())).tocsc()
    return counts, metadata, barcodes, features, spatial


def reference_signatures(cell_proportions, type_cell_indices):
    signatures = np.column_stack(
        [cell_proportions[indices].mean(axis=0) for indices in type_cell_indices]
    )
    signatures /= signatures.sum(axis=0, keepdims=True)
    return signatures


def select_markers(signatures, markers_per_type):
    selected = set()
    for type_index in range(signatures.shape[1]):
        other_mean = np.delete(signatures, type_index, axis=1).mean(axis=1)
        specificity = np.log2(
            (signatures[:, type_index] + 1e-8) / (other_mean + 1e-8)
        )
        selected.update(np.argsort(specificity, kind="stable")[-markers_per_type:])
    return np.asarray(sorted(selected))


def fit_mixture(observed, signatures):
    weights = np.repeat(1.0 / signatures.shape[1], signatures.shape[1])
    total = observed.sum()
    for _ in range(5_000):
        probability = np.maximum(signatures @ weights, 1e-15)
        updated = weights * (signatures.T @ (observed / probability)) / total
        updated /= updated.sum()
        if np.max(np.abs(updated - weights)) < 1e-10:
            weights = updated
            break
        weights = updated
    return weights


def marker_fit(observed, signatures, markers_per_type):
    selected = select_markers(signatures, markers_per_type)
    marker_signatures = signatures[selected].copy()
    marker_signatures /= marker_signatures.sum(axis=0, keepdims=True)
    weights = fit_mixture(observed[selected], marker_signatures)
    return weights, selected, marker_signatures


def cosine(left, right):
    return float(left @ right / (np.linalg.norm(left) * np.linalg.norm(right)))


def main():
    repo = Path(__file__).resolve().parents[5]
    archive = repo / "inputs" / "ls04-spatial-deconvolution" / "spatial.sim.tar.gz"
    counts, metadata, barcodes, features, spatial = load_input(archive)

    if not counts.index.equals(metadata.index):
        raise ValueError("Reference count and metadata indices differ")
    if [feature[1] for feature in features] != list(counts.columns):
        raise ValueError("Reference/spatial genes are not identically ordered")
    if spatial.shape != (len(features), len(barcodes)):
        raise ValueError("Spatial matrix dimensions disagree with annotations")

    raw_cells = counts.to_numpy()
    if raw_cells.dtype.kind not in "iu" or np.any(raw_cells < 0):
        raise ValueError("Reference is not raw non-negative integer counts")
    cell_library = raw_cells.sum(axis=1, keepdims=True)
    cell_proportions = raw_cells / cell_library
    cell_types = sorted(metadata["cell_type"].unique())
    type_cell_indices = [
        np.flatnonzero(metadata["cell_type"].to_numpy() == cell_type)
        for cell_type in cell_types
    ]
    signatures = reference_signatures(cell_proportions, type_cell_indices)
    observed = np.asarray(spatial[:, barcodes.index(TARGET_SPOT)].todense()).ravel()

    raw_weights, selected, marker_signatures = marker_fit(
        observed, signatures, PRIMARY_MARKERS_PER_TYPE
    )
    sensitivity = np.vstack(
        [marker_fit(observed, signatures, count)[0] for count in MARKER_SENSITIVITY]
    )

    # Joint parametric bootstrap: resample target counts and the 200 reference
    # cells within every type, then reselect markers and refit from scratch.
    fitted_full = signatures @ raw_weights
    rng = np.random.default_rng(RANDOM_SEED)
    bootstrap_weights = []
    for _ in range(N_BOOTSTRAPS):
        target_bootstrap = rng.multinomial(int(observed.sum()), fitted_full)
        sampled_signatures = np.column_stack(
            [
                cell_proportions[
                    rng.choice(indices, size=len(indices), replace=True)
                ].mean(axis=0)
                for indices in type_cell_indices
            ]
        )
        sampled_signatures /= sampled_signatures.sum(axis=0, keepdims=True)
        bootstrap_weights.append(
            marker_fit(
                target_bootstrap, sampled_signatures, PRIMARY_MARKERS_PER_TYPE
            )[0]
        )
    bootstrap_weights = np.vstack(bootstrap_weights)
    intervals = np.quantile(bootstrap_weights, (0.025, 0.975), axis=0)

    retained = raw_weights >= REPORT_THRESHOLD
    reported = raw_weights[retained]
    reported /= reported.sum()

    marker_observed = observed[selected] / observed[selected].sum()
    marker_fit_profile = marker_signatures @ raw_weights
    mixture_cosine = cosine(marker_fit_profile, marker_observed)
    best_pure_cosine = max(
        cosine(marker_signatures[:, index], marker_observed)
        for index in range(len(cell_types))
    )

    rows = []
    for weight, type_index in zip(reported, np.flatnonzero(retained)):
        evidence = (
            f"marker_MLE={raw_weights[type_index]:.3f}; "
            f"joint_bootstrap_95%CI=[{intervals[0, type_index]:.3f},"
            f"{intervals[1, type_index]:.3f}]; "
            f"marker_sensitivity=[{sensitivity[:, type_index].min():.3f},"
            f"{sensitivity[:, type_index].max():.3f}]"
        )
        rows.append(
            {
                "cell_type": cell_types[type_index],
                "weight": f"{weight:.6f}",
                "evidence": evidence,
            }
        )

    output = Path(__file__).with_name("spot_710_composition.csv")
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("cell_type", "weight", "evidence"))
        writer.writeheader()
        writer.writerows(rows)

    print(f"selected_markers={len(selected)} spot_total={int(observed.sum())}")
    print(f"mixture_cosine={mixture_cosine:.6f} best_pure_cosine={best_pure_cosine:.6f}")
    print(f"reported_weight_sum={sum(float(row['weight']) for row in rows):.6f}")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
