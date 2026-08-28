#!/usr/bin/env python3
"""Reproducible enhancer-promoter distance and lag-association analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


CONTACT_THRESHOLD_NM = 260.0
LAG_VALUES = np.arange(-20, 21, dtype=int)
COLUMNS = [
    "cell_id", "time", "enh_x", "enh_y", "enh_z",
    "prom_x", "prom_y", "prom_z", "transcription",
]


def correlation_from_sums(
    n: float, sx: float, sy: float, sxx: float, syy: float, sxy: float
) -> float:
    numerator = n * sxy - sx * sy
    denominator = np.sqrt((n * sxx - sx * sx) * (n * syy - sy * sy))
    return float(numerator / denominator) if denominator > 0 else float("nan")


def paired_at_lag(x: np.ndarray, y: np.ndarray, lag: int) -> tuple[np.ndarray, np.ndarray]:
    """Pair x(t) with y(t+lag) independently within each matrix row/cell."""
    if lag > 0:
        return x[:, :-lag], y[:, lag:]
    if lag < 0:
        k = -lag
        return x[:, k:], y[:, :-k]
    return x, y


def validate_and_derive(input_csv: Path) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    data = pd.read_csv(input_csv)
    if list(data.columns) != COLUMNS:
        raise ValueError(f"Expected {COLUMNS}, received {list(data.columns)}")
    if data.empty or data.isna().any().any():
        raise ValueError("The input must be nonempty and complete")
    if data.duplicated(["cell_id", "time"]).any():
        raise ValueError("Each (cell_id, time) pair must be unique")
    if not set(data["transcription"].unique()).issubset({0, 1}):
        raise ValueError("transcription must contain only 0 and 1")

    data = data.sort_values(["cell_id", "time"], kind="mergesort").reset_index(drop=True)
    sizes = data.groupby("cell_id", sort=True).size()
    if sizes.nunique() != 1 or int(sizes.iloc[0]) != 250:
        raise ValueError("All cells must have exactly 250 rows")
    for cell_id, frame in data.groupby("cell_id", sort=True):
        t = frame["time"].to_numpy()
        if t.size != 250 or np.unique(t).size != 250 or np.any(np.diff(t) <= 0):
            raise ValueError(f"Times are not unique and increasing in cell {cell_id}")

    enhancer = data[["enh_x", "enh_y", "enh_z"]].to_numpy(dtype=float)
    promoter = data[["prom_x", "prom_y", "prom_z"]].to_numpy(dtype=float)
    data["distance_nm"] = np.linalg.norm(enhancer - promoter, axis=1)
    data["contact"] = (data["distance_nm"] <= CONTACT_THRESHOLD_NM).astype(np.uint8)

    n_cells = sizes.size
    contact = data["contact"].to_numpy().reshape(n_cells, 250)
    transcription = data["transcription"].to_numpy(dtype=np.uint8).reshape(n_cells, 250)
    return data, contact, transcription


def calculate_cell_metrics(data: pd.DataFrame) -> pd.DataFrame:
    result = data.groupby("cell_id", sort=True, as_index=False).agg(
        n_timepoints=("time", "count"),
        mean_distance_nm=("distance_nm", "mean"),
        contact_fraction=("contact", "mean"),
        transcription_fraction=("transcription", "mean"),
    )
    return result[
        ["cell_id", "n_timepoints", "mean_distance_nm", "contact_fraction", "transcription_fraction"]
    ]


def calculate_lag_table(contact: np.ndarray, transcription: np.ndarray) -> pd.DataFrame:
    rows: list[dict[str, float | int]] = []
    for lag in LAG_VALUES:
        x, y = paired_at_lag(contact, transcription, int(lag))
        xf = x.astype(float, copy=False).ravel()
        yf = y.astype(float, copy=False).ravel()
        r = correlation_from_sums(
            float(xf.size), xf.sum(), yf.sum(), np.dot(xf, xf), np.dot(yf, yf), np.dot(xf, yf)
        )
        rows.append({"lag": int(lag), "association": r, "n_observations": int(xf.size)})
    return pd.DataFrame(rows, columns=["lag", "association", "n_observations"])


def cell_cluster_bootstrap(
    contact: np.ndarray, transcription: np.ndarray, lag: int, iterations: int = 2000
) -> tuple[float, float]:
    """Bootstrap whole cells for a conditional interval at one preselected lag."""
    x, y = paired_at_lag(contact, transcription, lag)
    xf = x.astype(float)
    yf = y.astype(float)
    per_cell = np.column_stack(
        [
            np.full(x.shape[0], x.shape[1], dtype=float),
            xf.sum(axis=1),
            yf.sum(axis=1),
            np.square(xf).sum(axis=1),
            np.square(yf).sum(axis=1),
            (xf * yf).sum(axis=1),
        ]
    )
    rng = np.random.default_rng(827)
    estimates = np.empty(iterations, dtype=float)
    for i in range(iterations):
        totals = per_cell[rng.integers(0, x.shape[0], size=x.shape[0])].sum(axis=0)
        estimates[i] = correlation_from_sums(*totals)
    return tuple(float(v) for v in np.nanpercentile(estimates, [2.5, 97.5]))


def make_lag_figure(lag_table: pd.DataFrame, out_dir: Path) -> None:
    """Write a compact, dependency-free, colorblind-safe vector figure."""
    strongest = lag_table.loc[lag_table["association"].abs().idxmax()]
    lag = int(strongest["lag"])
    r = float(strongest["association"])
    width, height = 720, 430
    left, right, top, bottom = 90, 25, 35, 85
    plot_w, plot_h = width - left - right, height - top - bottom
    y_lim = max(0.1, float(lag_table["association"].abs().max()) * 1.18)

    def px(x: float) -> float:
        return left + (x + 20.0) / 40.0 * plot_w

    def py(y: float) -> float:
        return top + (y_lim - y) / (2.0 * y_lim) * plot_h

    points = " ".join(
        f"{px(float(row.lag)):.2f},{py(float(row.association)):.2f}"
        for row in lag_table.itertuples(index=False)
    )
    x_ticks = "".join(
        f'<line x1="{px(x):.2f}" y1="{top + plot_h}" x2="{px(x):.2f}" y2="{top + plot_h + 5}" stroke="#222"/>'
        f'<text x="{px(x):.2f}" y="{top + plot_h + 21}" text-anchor="middle">{x}</text>'
        for x in range(-20, 21, 5)
    )
    y_ticks = "".join(
        f'<line x1="{left - 5}" y1="{py(y):.2f}" x2="{left}" y2="{py(y):.2f}" stroke="#222"/>'
        f'<text x="{left - 9}" y="{py(y) + 4:.2f}" text-anchor="end">{y:.2f}</text>'
        for y in np.linspace(-y_lim, y_lim, 5)
    )
    circles = "".join(
        f'<circle cx="{px(float(row.lag)):.2f}" cy="{py(float(row.association)):.2f}" r="2.7" fill="#0072B2"/>'
        for row in lag_table.itertuples(index=False)
    )
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">
<title id="title">Lag profile of contact-transcription association</title>
<desc id="desc">Pearson association for every lag from minus twenty to plus twenty. The largest absolute association is marked at lag {lag:+d}, r {r:.3f}.</desc>
<rect width="100%" height="100%" fill="white"/>
<g font-family="Arial, Helvetica, sans-serif" font-size="12" fill="#222">
  <line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#222"/>
  <line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#222"/>
  <line x1="{left}" y1="{py(0):.2f}" x2="{left + plot_w}" y2="{py(0):.2f}" stroke="#777" stroke-dasharray="5,4"/>
  {x_ticks}{y_ticks}
  <polyline points="{points}" fill="none" stroke="#0072B2" stroke-width="2"/>
  {circles}
  <rect x="{px(lag) - 5:.2f}" y="{py(r) - 5:.2f}" width="10" height="10" fill="#D55E00" transform="rotate(45 {px(lag):.2f} {py(r):.2f})"/>
  <text x="{left + plot_w / 2:.2f}" y="{height - 28}" text-anchor="middle" font-size="14">Lag (time steps; positive = contact leads transcription)</text>
  <text x="24" y="{top + plot_h / 2:.2f}" text-anchor="middle" font-size="14" transform="rotate(-90 24 {top + plot_h / 2:.2f})">Pooled Pearson association (r)</text>
  <text x="{left + 8}" y="{top + 18}" fill="#D55E00">Largest |r|: lag {lag:+d}, r = {r:.3f}</text>
</g></svg>'''
    (out_dir / "lag_profile.svg").write_text(svg, encoding="utf-8")


def write_report(
    out_path: Path,
    data: pd.DataFrame,
    cell_metrics: pd.DataFrame,
    lag_table: pd.DataFrame,
    bootstrap_ci: tuple[float, float],
) -> None:
    strongest = lag_table.loc[lag_table["association"].abs().idxmax()]
    lag = int(strongest["lag"])
    r = float(strongest["association"])
    n = int(strongest["n_observations"])
    if lag > 0:
        temporal = f"contact precedes transcription by {lag} time step(s)"
    elif lag < 0:
        temporal = f"transcription precedes contact by {-lag} time step(s)"
    else:
        temporal = "contact and transcription are measured at the same time"

    overall_distance = float(data["distance_nm"].mean())
    overall_contact = float(data["contact"].mean())
    overall_tx = float(data["transcription"].mean())
    ci_low, ci_high = bootstrap_ci
    report = f"""# Enhancer-promoter contact and transcription dynamics

## Data and analysis

The dataset has {len(data):,} complete observations: {cell_metrics.shape[0]} cells with 250 ordered time points each. For every observation, enhancer-promoter distance is the Euclidean norm of the x/y/z coordinate difference in nanometers. Contact is 1 when distance is <= {CONTACT_THRESHOLD_NM:.0f} nm, inclusive, and 0 otherwise. Per-cell fractions use all 250 rows.

For each integer lag from -20 through +20, contact(t) is paired only with transcription(t + lag) from the same cell. The requested pooled Pearson correlation and exact number of aligned pairs are reported for every lag. Positive lag means contact leads later transcription. No cross-cell pairs are formed.

The overall mean distance is {overall_distance:.3f} nm, the pooled contact fraction is {overall_contact:.4f}, and the pooled transcription fraction is {overall_tx:.4f}. These pooled summaries do not replace the 600 per-cell values in `cell_metrics.csv`.

## Lag result

The largest absolute association is at lag {lag:+d}: r = {r:.6f}, n = {n:,}, r^2 = {r*r:.6f}. By the supplied sign convention, {temporal}. This is a weak association in magnitude. A 2,000-resample whole-cell bootstrap, conditional on examining this selected lag, gives a descriptive 95% interval of [{ci_low:.6f}, {ci_high:.6f}]. It preserves cell series during resampling, but it does not correct for choosing the maximum absolute value among 41 lags.

The complete numerical profile is in `lag_analysis.csv`. `lag_profile.svg` displays all 41 prespecified lags with a zero reference line, accessible text, and colorblind-safe styling; the plotted point count for each lag is available in the CSV rather than implied to be constant.

## Regulatory annotation scope

The available measurements are microscopy-style Cartesian positions in nanometers. They are not chromosome coordinates and provide no organism, chromosome, assembly, interval, or gene identifier. Consequently, genomic-region overlap, regulatory-element-to-gene, sequence, binding-matrix, or phenotype-region queries are not identifiable from this input and were not fabricated or run. Such annotation would require a real genomic interval plus organism/assembly metadata.

## Association is not causation

The lag pattern is observational. It cannot establish that enhancer-promoter contact causes transcription, that transcription causes subsequent contact, or that either direction is free of shared or time-varying causes. Repeated measurements are serially dependent, and the 149,400 pooled pairs at lag -1 are not 149,400 independent experimental units. The data cannot establish a molecular mechanism, rule out confounding, resolve direction below the sampling interval, or estimate an intervention effect. Causal conclusions require a suitable perturbation or identification design, appropriate controls, and inference that accounts for cells and temporal dependence.
"""
    out_path.write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    data, contact, transcription = validate_and_derive(args.input_csv)
    cell_metrics = calculate_cell_metrics(data)
    lag_table = calculate_lag_table(contact, transcription)
    strongest_lag = int(lag_table.loc[lag_table["association"].abs().idxmax(), "lag"])
    ci = cell_cluster_bootstrap(contact, transcription, strongest_lag)

    cell_metrics.to_csv(args.output_dir / "cell_metrics.csv", index=False, lineterminator="\n")
    lag_table.to_csv(args.output_dir / "lag_analysis.csv", index=False, lineterminator="\n")
    make_lag_figure(lag_table, args.output_dir)
    write_report(args.output_dir / "report.md", data, cell_metrics, lag_table, ci)


if __name__ == "__main__":
    main()
