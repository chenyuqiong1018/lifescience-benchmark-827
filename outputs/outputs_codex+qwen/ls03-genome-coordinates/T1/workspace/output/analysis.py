#!/usr/bin/env python3
"""Enhancer-promoter 3D distance and transcription dynamics analysis.

Reads single-cell time-series enhancer/promoter coordinates (nm) and binary
transcription states, then computes:

  1. Per-cell metrics -> output/cell_metrics.csv
       cell_id, n_timepoints, mean_distance_nm, contact_fraction,
       transcription_fraction
     - Euclidean distance from the supplied x/y/z coordinates (nm).
     - Contact = distance <= 260 nm (supplied contact threshold).
     - Fractions use all 250 rows per cell.

  2. Lagged association -> output/lag_analysis.csv
       lag, association, n_observations
     - One row per integer lag from -20 through +20.
     - association = pooled Pearson correlation across all cells between
       contact at time t and transcription at time t + lag.
     - Positive lag means contact leads the later transcription value.
     - Observations are never joined across cell boundaries.

Per inputs/ANALYSIS_RULE.md, the strongest absolute association is a temporal
association only and is never interpreted as proof of causation.
"""

import os

import numpy as np
import pandas as pd

CONTACT_THRESHOLD_NM = 260.0
LAG_MIN = -20
LAG_MAX = 20

BASE_DIR = os.path.dirname(os.path.abspath(__file__))          # .../output
WORKSPACE_DIR = os.path.dirname(BASE_DIR)                      # workspace root
INPUT_CSV = os.path.join(WORKSPACE_DIR, "inputs", "single_cell_dynamics_question.csv")
CELL_METRICS_CSV = os.path.join(BASE_DIR, "cell_metrics.csv")
LAG_ANALYSIS_CSV = os.path.join(BASE_DIR, "lag_analysis.csv")


def pooled_pearson(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson correlation of two pooled 1-D arrays; NaN if degenerate."""
    if len(x) < 2:
        return float("nan")
    x = x - x.mean()
    y = y - y.mean()
    denom = float(np.sqrt((x * x).sum()) * np.sqrt((y * y).sum()))
    if denom == 0.0:
        return float("nan")
    return float((x * y).sum() / denom)


def main() -> None:
    df = pd.read_csv(INPUT_CSV)

    required = ["cell_id", "time", "enh_x", "enh_y", "enh_z",
                "prom_x", "prom_y", "prom_z", "transcription"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"missing columns: {missing}")

    # ---- Euclidean enhancer-promoter distance (nm) and contact ----
    delta = (df[["enh_x", "enh_y", "enh_z"]].to_numpy(dtype=np.float64)
             - df[["prom_x", "prom_y", "prom_z"]].to_numpy(dtype=np.float64))
    df["distance_nm"] = np.sqrt((delta ** 2).sum(axis=1))
    df["contact"] = (df["distance_nm"] <= CONTACT_THRESHOLD_NM).astype(np.int64)

    # ---- Per-cell metrics (fractions over all rows of the cell) ----
    g = df.groupby("cell_id", sort=True)
    cell_metrics = pd.DataFrame({
        "cell_id": g["time"].count().index,
        "n_timepoints": g["time"].count().to_numpy(),
        "mean_distance_nm": g["distance_nm"].mean().to_numpy(),
        "contact_fraction": g["contact"].mean().to_numpy(),
        "transcription_fraction": g["transcription"].mean().to_numpy(),
    })
    cell_metrics.to_csv(CELL_METRICS_CSV, index=False, float_format="%.6f")

    # ---- Lagged association (never joining across cell boundaries) ----
    # Build per-cell matrices sorted by (cell_id, time). Each cell has the
    # same number of timepoints with consecutive integer times starting at 0,
    # so pairing time t with t + lag is a positional shift within each cell.
    df_sorted = df.sort_values(["cell_id", "time"], kind="mergesort")
    sizes = df_sorted.groupby("cell_id", sort=True)["time"].count().to_numpy()
    n_cells = len(sizes)
    if sizes.min() != sizes.max():
        raise ValueError("cells have unequal timepoint counts; lag pairing by "
                         "positional shift is invalid")
    T = int(sizes[0])

    # Validate consecutive integer times within each cell.
    first_t = df_sorted.groupby("cell_id", sort=True)["time"].transform("first")
    rank = df_sorted.groupby("cell_id").cumcount()
    if not ((df_sorted["time"] - first_t) == rank).all():
        raise ValueError("timepoints are not consecutive integers within cells")

    contact_mat = df_sorted["contact"].to_numpy(dtype=np.float64).reshape(n_cells, T)
    trans_mat = df_sorted["transcription"].to_numpy(dtype=np.float64).reshape(n_cells, T)

    lag_rows = []
    for lag in range(LAG_MIN, LAG_MAX + 1):
        if lag >= 0:
            x = contact_mat[:, :T - lag].ravel()      # contact at t
            y = trans_mat[:, lag:].ravel()             # transcription at t + lag
        else:
            x = contact_mat[:, -lag:].ravel()          # contact at t
            y = trans_mat[:, :T + lag].ravel()         # transcription at t + lag
        lag_rows.append({
            "lag": lag,
            "association": pooled_pearson(x, y),
            "n_observations": int(len(x)),
        })

    lag_df = pd.DataFrame(lag_rows, columns=["lag", "association", "n_observations"])
    lag_df.to_csv(LAG_ANALYSIS_CSV, index=False, float_format="%.6f")

    # ---- Console summary ----
    print(f"cells={n_cells} timepoints_per_cell={T} total_rows={len(df)}")
    print(f"contact_threshold_nm={CONTACT_THRESHOLD_NM}")
    print(f"overall_mean_distance_nm={df['distance_nm'].mean():.3f}")
    print(f"overall_contact_fraction={df['contact'].mean():.6f}")
    print(f"overall_transcription_fraction={df['transcription'].mean():.6f}")
    print("per-cell contact_fraction min/median/max: "
          f"{cell_metrics['contact_fraction'].min():.4f}/"
          f"{cell_metrics['contact_fraction'].median():.4f}/"
          f"{cell_metrics['contact_fraction'].max():.4f}")
    print("per-cell transcription_fraction min/median/max: "
          f"{cell_metrics['transcription_fraction'].min():.4f}/"
          f"{cell_metrics['transcription_fraction'].median():.4f}/"
          f"{cell_metrics['transcription_fraction'].max():.4f}")
    cell_level_r = pooled_pearson(cell_metrics["contact_fraction"].to_numpy(),
                                  cell_metrics["transcription_fraction"].to_numpy())
    print(f"cell-level corr(contact_fraction, transcription_fraction)={cell_level_r:.6f}")

    a0 = lag_df.loc[lag_df["lag"] == 0, "association"].iloc[0]
    best = lag_df.loc[lag_df["association"].abs().idxmax()]
    print(f"lag0_association={a0:.6f}")
    print(f"strongest |association|: lag={int(best['lag'])} r={best['association']:.6f} "
          f"n={int(best['n_observations'])}")
    for q in (-20, -10, -5, -1, 1, 5, 10, 20):
        val = lag_df.loc[lag_df["lag"] == q, "association"].iloc[0]
        print(f"  lag={q:+d} association={val:.6f}")
    print(f"wrote {CELL_METRICS_CSV} ({len(cell_metrics)} rows)")
    print(f"wrote {LAG_ANALYSIS_CSV} ({len(lag_df)} rows)")


if __name__ == "__main__":
    main()

