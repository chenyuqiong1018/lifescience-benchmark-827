"""Leakage-free composition estimate for Visium spot Spot_710-1."""

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
REPORT_THRESHOLD = 0.02
N_BOOTSTRAPS = 200
RANDOM_SEED = 120827


def read_archive(archive_path: Path):
    with tarfile.open(archive_path, "r:gz") as archive:
        counts = pd.read_csv(archive.extractfile("spatial_q_sc_counts.csv"), index_col=0)
        metadata = pd.read_csv(
            archive.extractfile("spatial_q_sc_metadata.csv"), index_col=0
        )

        barcode_member = archive.extractfile("barcodes.tsv.gz")
        with gzip.GzipFile(fileobj=barcode_member) as compressed:
            barcodes = [line.decode().strip() for line in compressed]

        feature_member = archive.extractfile("features.tsv.gz")
        with gzip.GzipFile(fileobj=feature_member) as compressed:
            features = [line.decode().rstrip("\n").split("\t") for line in compressed]

        matrix_member = archive.extractfile("matrix.mtx.gz")
        with gzip.GzipFile(fileobj=matrix_member) as compressed:
            matrix_bytes = compressed.read()
        spatial = mmread(io.BytesIO(matrix_bytes)).tocsc()

    return counts, metadata, barcodes, features, spatial


def build_reference(counts: pd.DataFrame, metadata: pd.DataFrame):
    if not counts.index.equals(metadata.index):
        raise ValueError("Single-cell count and metadata indices differ")
    cell_types = sorted(metadata["cell_type"].unique())
    library_sizes = counts.sum(axis=1)
    if np.any(library_sizes <= 0):
        raise ValueError("Single-cell reference contains an empty library")
    cell_proportions = counts.div(library_sizes, axis=0)
    signatures = np.column_stack(
        [
            cell_proportions.loc[metadata.index[metadata["cell_type"] == cell_type]]
            .mean(axis=0)
            .to_numpy()
            for cell_type in cell_types
        ]
    )
    signatures /= signatures.sum(axis=0, keepdims=True)
    return cell_types, signatures


def fit_multinomial_mixture(observed: np.ndarray, signatures: np.ndarray):
    """EM maximum likelihood for a convex mixture of type signatures."""
    weights = np.repeat(1.0 / signatures.shape[1], signatures.shape[1])
    total = observed.sum()
    for _ in range(5_000):
        fitted = np.maximum(signatures @ weights, 1e-15)
        updated = weights * (signatures.T @ (observed / fitted)) / total
        updated /= updated.sum()
        if np.max(np.abs(updated - weights)) < 1e-10:
            weights = updated
            break
        weights = updated
    return weights


def cosine(left: np.ndarray, right: np.ndarray):
    return float(left @ right / (np.linalg.norm(left) * np.linalg.norm(right)))


def main():
    repo = Path(__file__).resolve().parents[5]
    archive_path = (
        repo
        / "inputs"
        / "ls04-spatial-deconvolution"
        / "spatial.sim.tar.gz"
    )
    counts, metadata, barcodes, features, spatial = read_archive(archive_path)

    if [row[1] for row in features] != list(counts.columns):
        raise ValueError("Spatial and single-cell gene names/order differ")
    if spatial.shape != (len(features), len(barcodes)):
        raise ValueError("Spatial matrix dimensions do not match annotations")
    if TARGET_SPOT not in barcodes:
        raise ValueError(f"Missing {TARGET_SPOT}")

    cell_types, signatures = build_reference(counts, metadata)
    observed = np.asarray(spatial[:, barcodes.index(TARGET_SPOT)].todense()).ravel()
    if observed.sum() <= 0:
        raise ValueError("Target spot has no counts")
    raw_weights = fit_multinomial_mixture(observed, signatures)
    fitted = signatures @ raw_weights

    rng = np.random.default_rng(RANDOM_SEED)
    bootstrap_weights = np.vstack(
        [
            fit_multinomial_mixture(
                rng.multinomial(int(observed.sum()), fitted), signatures
            )
            for _ in range(N_BOOTSTRAPS)
        ]
    )
    intervals = np.quantile(bootstrap_weights, (0.025, 0.975), axis=0)

    retained = raw_weights >= REPORT_THRESHOLD
    if not np.any(retained):
        retained[np.argmax(raw_weights)] = True
    reported = raw_weights[retained]
    reported /= reported.sum()

    observed_proportion = observed / observed.sum()
    mixture_cosine = cosine(fitted, observed_proportion)
    pure_cosines = np.asarray(
        [cosine(signatures[:, column], observed_proportion) for column in range(len(cell_types))]
    )

    rows = []
    for reported_weight, type_index in zip(reported, np.flatnonzero(retained)):
        evidence = (
            f"unthresholded_MLE={raw_weights[type_index]:.3f}; "
            f"parametric_bootstrap_95%CI=[{intervals[0, type_index]:.3f},"
            f"{intervals[1, type_index]:.3f}]; "
            f"mixture_cosine={mixture_cosine:.3f}; best_pure_cosine={pure_cosines.max():.3f}"
        )
        rows.append(
            {
                "cell_type": cell_types[type_index],
                "weight": f"{reported_weight:.6f}",
                "evidence": evidence,
            }
        )

    output_path = Path(__file__).with_name("spot_710_composition.csv")
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("cell_type", "weight", "evidence"))
        writer.writeheader()
        writer.writerows(rows)

    print(f"spot_total={int(observed.sum())} genes={len(features)}")
    print(f"mixture_cosine={mixture_cosine:.6f} best_pure_cosine={pure_cosines.max():.6f}")
    print(f"reported_weight_sum={sum(float(row['weight']) for row in rows):.6f}")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
