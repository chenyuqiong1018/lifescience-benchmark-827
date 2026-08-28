"""Controlled-skill multi-view deconvolution for Spot_710-1.

The unavailable scVI/cell2location and scGPT runtimes are represented by
explicit count and rank proxies.  The script reports exactly what was run and
creates a publication-style uncertainty diagnostic.
"""

from __future__ import annotations

import csv
import gzip
import io
import tarfile
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.optimize import nnls
from scipy.stats import rankdata


TARGET_SPOT = "Spot_710-1"
MARKERS_PER_TYPE = 100
REPORT_THRESHOLD = 0.02
N_BOOTSTRAPS = 200
RANDOM_SEED = 122827
OKABE_ITO = ("#56B4E9", "#E69F00", "#009E73")


def load_archive(path: Path):
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


def make_signatures(cell_proportions, type_indices):
    signatures = np.column_stack(
        [cell_proportions[indices].mean(axis=0) for indices in type_indices]
    )
    signatures /= signatures.sum(axis=0, keepdims=True)
    return signatures


def multinomial_fit(observed, signatures):
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


def select_markers(signatures):
    selected = set()
    for type_index in range(signatures.shape[1]):
        other = np.delete(signatures, type_index, axis=1).mean(axis=1)
        specificity = np.log2(
            (signatures[:, type_index] + 1e-8) / (other + 1e-8)
        )
        selected.update(
            np.argsort(specificity, kind="stable")[-MARKERS_PER_TYPE:]
        )
    return np.asarray(sorted(selected))


def marker_count_fit(observed, signatures):
    selected = select_markers(signatures)
    subset = signatures[selected].copy()
    subset /= subset.sum(axis=0, keepdims=True)
    return multinomial_fit(observed[selected], subset)


def rank_proxy_fit(observed, signatures):
    """scGPT-guided scale-robust proxy for synthetic, out-of-vocab genes."""
    ranked_signatures = rankdata(signatures, axis=0, method="average")
    ranked_signatures /= ranked_signatures.sum(axis=0, keepdims=True)
    ranked_observed = rankdata(observed, method="average")
    ranked_observed /= ranked_observed.sum()
    weights, _ = nnls(ranked_signatures, ranked_observed)
    if weights.sum() == 0:
        raise ValueError("Rank proxy returned an empty mixture")
    return weights / weights.sum()


def three_view_fit(observed, signatures):
    all_gene = multinomial_fit(observed, signatures)
    marker = marker_count_fit(observed, signatures)
    ranked = rank_proxy_fit(observed, signatures)
    ensemble = np.vstack((all_gene, marker, ranked)).mean(axis=0)
    ensemble /= ensemble.sum()
    return ensemble, np.vstack((all_gene, marker, ranked))


def cosine(left, right):
    return float(left @ right / (np.linalg.norm(left) * np.linalg.norm(right)))


def save_diagnostic(output_dir, names, weights, intervals):
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8,
            "axes.labelsize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "svg.hashsalt": "task12-t1",
        }
    )
    figure, axis = plt.subplots(figsize=(3.5, 2.8), constrained_layout=True)
    x = np.arange(len(names))
    errors = np.vstack((weights - intervals[0], intervals[1] - weights))
    bars = axis.bar(
        x,
        weights,
        color=OKABE_ITO[: len(names)],
        edgecolor="black",
        linewidth=0.7,
        yerr=errors,
        error_kw={"ecolor": "black", "elinewidth": 0.8, "capsize": 3},
    )
    for bar, hatch in zip(bars, ("//", "..", "xx")):
        bar.set_hatch(hatch)
    axis.set_xticks(x, [name.replace("_", "\n") for name in names])
    axis.set_ylabel("Estimated composition (fraction)")
    axis.set_ylim(0, max(0.45, float(intervals[1].max() + 0.06)))
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.text(
        0.01,
        0.97,
        "Error bars: joint bootstrap 95% CI (n=200)",
        transform=axis.transAxes,
        va="top",
        fontsize=6.5,
    )
    figure.savefig(
        output_dir / "composition_diagnostic.png",
        dpi=300,
        metadata={"Software": "matplotlib"},
    )
    figure.savefig(
        output_dir / "composition_diagnostic.svg",
        metadata={"Date": None, "Creator": "matplotlib"},
    )
    plt.close(figure)


def main():
    repo = Path(__file__).resolve().parents[5]
    archive = repo / "inputs" / "ls04-spatial-deconvolution" / "spatial.sim.tar.gz"
    counts, metadata, barcodes, features, spatial = load_archive(archive)

    if not counts.index.equals(metadata.index):
        raise ValueError("Reference count/metadata indices differ")
    if [feature[1] for feature in features] != list(counts.columns):
        raise ValueError("Spatial/reference genes are not identically ordered")
    if spatial.shape != (len(features), len(barcodes)):
        raise ValueError("Spatial matrix dimensions disagree with annotations")

    raw_cells = counts.to_numpy()
    if raw_cells.dtype.kind not in "iu" or np.any(raw_cells < 0):
        raise ValueError("Reference matrix is not raw non-negative integer counts")
    cell_proportions = raw_cells / raw_cells.sum(axis=1, keepdims=True)
    cell_types = sorted(metadata["cell_type"].unique())
    labels = metadata["cell_type"].to_numpy()
    type_indices = [np.flatnonzero(labels == cell_type) for cell_type in cell_types]
    signatures = make_signatures(cell_proportions, type_indices)
    observed = np.asarray(spatial[:, barcodes.index(TARGET_SPOT)].todense()).ravel()

    ensemble, view_weights = three_view_fit(observed, signatures)

    rng = np.random.default_rng(RANDOM_SEED)
    fitted_probability = signatures @ ensemble
    bootstrap = []
    for _ in range(N_BOOTSTRAPS):
        target_sample = rng.multinomial(int(observed.sum()), fitted_probability)
        sampled_signatures = make_signatures(
            cell_proportions,
            [
                rng.choice(indices, size=len(indices), replace=True)
                for indices in type_indices
            ],
        )
        bootstrap.append(three_view_fit(target_sample, sampled_signatures)[0])
    bootstrap = np.vstack(bootstrap)

    retained = ensemble >= REPORT_THRESHOLD
    reported = ensemble[retained]
    reported /= reported.sum()
    bootstrap_reported = bootstrap[:, retained]
    bootstrap_reported /= bootstrap_reported.sum(axis=1, keepdims=True)
    intervals = np.quantile(bootstrap_reported, (0.025, 0.975), axis=0)

    observed_proportion = observed / observed.sum()
    mixture_profile = signatures @ ensemble
    mixture_cosine = cosine(mixture_profile, observed_proportion)
    pure_cosines = [
        cosine(signatures[:, index], observed_proportion)
        for index in range(len(cell_types))
    ]
    mixture_log_likelihood = float(
        np.sum(observed * np.log(np.maximum(mixture_profile, 1e-15)))
    )
    pure_log_likelihood = max(
        float(np.sum(observed * np.log(np.maximum(signatures[:, index], 1e-15))))
        for index in range(len(cell_types))
    )

    rows = []
    for output_index, type_index in enumerate(np.flatnonzero(retained)):
        evidence = (
            f"views_all-marker-rank=[{view_weights[0, type_index]:.3f},"
            f"{view_weights[1, type_index]:.3f},{view_weights[2, type_index]:.3f}]; "
            f"joint_bootstrap_95%CI=[{intervals[0, output_index]:.3f},"
            f"{intervals[1, output_index]:.3f}]"
        )
        rows.append(
            {
                "cell_type": cell_types[type_index],
                "weight": f"{reported[output_index]:.6f}",
                "evidence": evidence,
            }
        )

    output_dir = Path(__file__).parent
    with (output_dir / "spot_710_composition.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=("cell_type", "weight", "evidence"))
        writer.writeheader()
        writer.writerows(rows)

    save_diagnostic(
        output_dir,
        [cell_types[index] for index in np.flatnonzero(retained)],
        reported,
        intervals,
    )

    print(f"spot_total={int(observed.sum())} reported_weight_sum={reported.sum():.6f}")
    print(f"mixture_cosine={mixture_cosine:.6f} best_pure_cosine={max(pure_cosines):.6f}")
    print(f"log_likelihood_gain_vs_best_pure={mixture_log_likelihood-pure_log_likelihood:.3f}")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
