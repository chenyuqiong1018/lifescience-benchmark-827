#!/usr/bin/env python3
"""Analyze enhancer-promoter distance, contact, and transcription by cell/time."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = [
    "cell_id", "time", "enh_x", "enh_y", "enh_z",
    "prom_x", "prom_y", "prom_z", "transcription",
]
CONTACT_THRESHOLD_NM = 260.0
LAGS = range(-20, 21)


def pearson_binary(x: np.ndarray, y: np.ndarray) -> float:
    """Return Pearson correlation, or NaN when either input has zero variance."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def load_and_validate(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if list(df.columns) != REQUIRED_COLUMNS:
        raise ValueError(f"Unexpected columns: {list(df.columns)}")
    if df.empty or df[REQUIRED_COLUMNS].isna().any().any():
        raise ValueError("Input must be non-empty and contain no missing values")
    if df.duplicated(["cell_id", "time"]).any():
        raise ValueError("Duplicate (cell_id, time) rows are not allowed")
    if not set(pd.unique(df["transcription"])).issubset({0, 1}):
        raise ValueError("transcription must be binary (0/1)")

    df = df.sort_values(["cell_id", "time"], kind="stable").reset_index(drop=True)
    counts = df.groupby("cell_id", sort=True).size()
    if not (counts == 250).all():
        raise ValueError("Every cell must contain exactly 250 time points")
    for cell_id, group in df.groupby("cell_id", sort=True):
        times = group["time"].to_numpy()
        if len(np.unique(times)) != 250 or not np.all(np.diff(times) > 0):
            raise ValueError(f"Invalid time series for cell {cell_id}")

    delta = df[["enh_x", "enh_y", "enh_z"]].to_numpy(float) - df[
        ["prom_x", "prom_y", "prom_z"]
    ].to_numpy(float)
    df["distance_nm"] = np.sqrt(np.sum(delta * delta, axis=1))
    df["contact"] = (df["distance_nm"] <= CONTACT_THRESHOLD_NM).astype(np.int8)
    return df


def make_cell_metrics(df: pd.DataFrame) -> pd.DataFrame:
    metrics = (
        df.groupby("cell_id", sort=True, as_index=False)
        .agg(
            n_timepoints=("time", "size"),
            mean_distance_nm=("distance_nm", "mean"),
            contact_fraction=("contact", "mean"),
            transcription_fraction=("transcription", "mean"),
        )
    )
    return metrics[
        ["cell_id", "n_timepoints", "mean_distance_nm", "contact_fraction", "transcription_fraction"]
    ]


def make_lag_analysis(df: pd.DataFrame) -> pd.DataFrame:
    contacts = []
    transcripts = []
    for _, group in df.groupby("cell_id", sort=True):
        contacts.append(group["contact"].to_numpy(np.int8))
        transcripts.append(group["transcription"].to_numpy(np.int8))
    contact = np.stack(contacts)
    transcription = np.stack(transcripts)

    rows = []
    for lag in LAGS:
        if lag > 0:
            x = contact[:, :-lag]
            y = transcription[:, lag:]
        elif lag < 0:
            k = -lag
            x = contact[:, k:]
            y = transcription[:, :-k]
        else:
            x = contact
            y = transcription
        rows.append(
            {
                "lag": lag,
                "association": pearson_binary(x.ravel(), y.ravel()),
                "n_observations": int(x.size),
            }
        )
    return pd.DataFrame(rows, columns=["lag", "association", "n_observations"])


def cellwise_summary_at_lag(df: pd.DataFrame, lag: int) -> tuple[int, float, float]:
    values = []
    for _, group in df.groupby("cell_id", sort=True):
        contact = group["contact"].to_numpy(np.int8)
        transcription = group["transcription"].to_numpy(np.int8)
        if lag > 0:
            x, y = contact[:-lag], transcription[lag:]
        elif lag < 0:
            k = -lag
            x, y = contact[k:], transcription[:-k]
        else:
            x, y = contact, transcription
        r = pearson_binary(x, y)
        if math.isfinite(r):
            values.append(r)
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return 0, float("nan"), float("nan")
    return int(arr.size), float(np.median(arr)), float(np.mean(arr > 0))


def write_report(
    path: Path, df: pd.DataFrame, metrics: pd.DataFrame, lag_table: pd.DataFrame
) -> None:
    strongest = lag_table.loc[lag_table["association"].abs().idxmax()]
    lag = int(strongest["lag"])
    r = float(strongest["association"])
    n = int(strongest["n_observations"])
    valid_cells, median_cell_r, positive_fraction = cellwise_summary_at_lag(df, lag)
    direction = (
        "contact leads later transcription" if lag > 0 else
        "transcription leads later contact" if lag < 0 else
        "contact and transcription are contemporaneous"
    )
    r2 = r * r

    text = f"""# Enhancer-promoter distance and transcription dynamics

## Methods

The input contains {len(df):,} observations from {metrics.shape[0]:,} cells, with 250 time points per cell and no cross-cell joins. Enhancer-promoter distance is the three-dimensional Euclidean distance in nanometers. A contact is defined by the supplied threshold, distance <= {CONTACT_THRESHOLD_NM:.0f} nm. Cell fractions use all 250 observations in each cell.

For every integer lag from -20 through +20, the analysis pools correctly aligned within-cell pairs and computes Pearson's correlation between contact at time *t* and transcription at time *t + lag*. Positive lag therefore means contact precedes transcription. Correlation itself is the requested effect-size measure; all 41 planned lags are reported without selecting only favorable results.

## Results

The largest absolute pooled temporal association occurs at lag {lag:+d}: r = {r:.6f}, n = {n:,}, r^2 = {r2:.6f}. Under the specified sign convention, this means {direction} by {abs(lag)} time step(s). The magnitude is weak: the variables share only about {100*r2:.3f}% of pooled variance at this selected lag. As a descriptive sensitivity check, {valid_cells} cells had nonconstant paired values at that lag; their median within-cell correlation was {median_cell_r:.6f}, and {100*positive_fraction:.1f}% were positive.

The complete lag profile is in `lag_analysis.csv`; cell-level distance, contact, and transcription summaries are in `cell_metrics.csv`.

## Interpretation and limitations

This is an observational temporal association, not evidence that contact causes transcription (or that transcription causes contact). Selecting the strongest absolute value across 41 lags is exploratory, and adjacent time points are repeated, autocorrelated measurements within cells, so the pooled observation count is not a count of independent experimental units. The data alone cannot establish intervention effects, exclude common or time-varying confounding, determine molecular mechanism, or prove temporal direction at a resolution finer than the sampling interval. A causal claim would require an appropriate perturbation or other identification strategy, controls, and uncertainty analysis that respects cell-level clustering and serial dependence.
"""
    path.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = load_and_validate(args.input_csv)
    metrics = make_cell_metrics(df)
    lag_table = make_lag_analysis(df)
    metrics.to_csv(args.output_dir / "cell_metrics.csv", index=False, lineterminator="\n")
    lag_table.to_csv(args.output_dir / "lag_analysis.csv", index=False, lineterminator="\n")
    write_report(args.output_dir / "report.md", df, metrics, lag_table)


if __name__ == "__main__":
    main()
