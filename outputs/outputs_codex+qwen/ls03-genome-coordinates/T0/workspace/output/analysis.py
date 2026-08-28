"""Enhancer-promoter 3D distance and transcription dynamics analysis.

Reads single-cell time-series data (inputs/single_cell_dynamics_question.csv) and
produces:
  - output/cell_metrics.csv   per-cell distance/contact/transcription summaries
  - output/lag_analysis.csv   pooled contact->transcription lag correlations (-20..+20)

Frozen rules (inputs/ANALYSIS_RULE.md):
  * Euclidean distance uses the supplied x/y/z coordinates in nm.
  * Contact is distance <= 260 nm (supplied contact threshold).
  * Per-cell fractions use all 250 rows per cell.
  * lag_analysis.csv has one row per integer lag -20..+20 with
    lag,association,n_observations where association is the POOLED Pearson
    correlation across all cells between contact at time t and transcription at
    time t+lag; positive lag means contact leads the later transcription value.
  * Observations are never joined across cell boundaries.
  * The strongest absolute association is a TEMPORAL ASSOCIATION only, never
    proof of causation.
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
INPUT_CSV = ROOT / "inputs" / "single_cell_dynamics_question.csv"
OUT_DIR = ROOT / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CONTACT_THRESHOLD_NM = 260.0  # supplied contact threshold
LAG_MIN, LAG_MAX = -20, 20    # inclusive lag range (integer lags)


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson correlation; returns NaN if either variable has zero variance."""
    if x.std() == 0.0 or y.std() == 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def main() -> None:
    df = pd.read_csv(INPUT_CSV)

    required = ["cell_id", "time", "enh_x", "enh_y", "enh_z",
                "prom_x", "prom_y", "prom_z", "transcription"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"input missing columns: {missing}")
    if df[required].isna().any().any():
        raise ValueError("input contains missing values")

    # ---- Euclidean enhancer-promoter distance (nm) and binary contact ----
    d = np.sqrt(
        (df["enh_x"] - df["prom_x"]) ** 2
        + (df["enh_y"] - df["prom_y"]) ** 2
        + (df["enh_z"] - df["prom_z"]) ** 2
    )
    df["distance_nm"] = d
    df["contact"] = (d <= CONTACT_THRESHOLD_NM).astype(int)

    # ---- Sanity checks on the time series structure ----
    counts = df.groupby("cell_id")["time"].size()
    n_cells = counts.size
    if not (counts == 250).all():
        raise ValueError("every cell must contribute exactly 250 rows")
    if not df.groupby("cell_id")["time"].apply(
        lambda s: s.is_monotonic_increasing and s.is_unique
    ).all():
        raise ValueError("time must be strictly increasing and unique within cell")

    # ---- Per-cell metrics (fractions use all 250 rows per cell) ----
    g = df.groupby("cell_id", sort=True)
    metrics = pd.DataFrame(
        {
            "cell_id": g["time"].size().index,
            "n_timepoints": g["time"].size().to_numpy(),
            "mean_distance_nm": g["distance_nm"].mean().to_numpy(),
            "contact_fraction": g["contact"].mean().to_numpy(),
            "transcription_fraction": g["transcription"].mean().to_numpy(),
        }
    )
    metrics.to_csv(OUT_DIR / "cell_metrics.csv", index=False)

    # ---- Lag analysis: pooled Pearson correlation, no cross-cell joins ----
    # Pivot to (cell x time) matrices; alignment is positional within each cell,
    # so pairs (contact at t, transcription at t+lag) never cross cell boundaries.
    contact = df.pivot(index="cell_id", columns="time", values="contact").to_numpy(dtype=float)
    transcr = df.pivot(index="cell_id", columns="time", values="transcription").to_numpy(dtype=float)
    assert contact.shape == transcr.shape
    n_time = contact.shape[1]

    lag_rows = []
    for lag in range(LAG_MIN, LAG_MAX + 1):
        if lag >= 0:
            c = contact[:, : n_time - lag].ravel()
            tr = transcr[:, lag:].ravel()
        else:
            c = contact[:, -lag:].ravel()
            tr = transcr[:, : n_time + lag].ravel()
        valid = np.isfinite(c) & np.isfinite(tr)
        assoc = pearson(c[valid], tr[valid])
        lag_rows.append({"lag": lag, "association": assoc, "n_observations": int(valid.sum())})
    lag_df = pd.DataFrame(lag_rows, columns=["lag", "association", "n_observations"])
    lag_df.to_csv(OUT_DIR / "lag_analysis.csv", index=False)

    # ---- Sensitivity: mean of per-cell (within-cell) Pearson correlations ----
    within_rows = []
    for lag in range(LAG_MIN, LAG_MAX + 1):
        if lag >= 0:
            c_mat, tr_mat = contact[:, : n_time - lag], transcr[:, lag:]
        else:
            c_mat, tr_mat = contact[:, -lag:], transcr[:, : n_time + lag]
        per_cell = []
        for i in range(contact.shape[0]):
            c_i, tr_i = c_mat[i], tr_mat[i]
            if c_i.std() > 0 and tr_i.std() > 0:
                per_cell.append(float(np.corrcoef(c_i, tr_i)[0, 1]))
        within_rows.append(
            {
                "lag": lag,
                "mean_within_cell_association": float(np.mean(per_cell)) if per_cell else float("nan"),
                "n_cells_defined": len(per_cell),
            }
        )
    within_df = pd.DataFrame(within_rows)

    # ---- Console summary (inputs for report.md) ----
    best = lag_df.loc[lag_df["association"].abs().idxmax()]
    print("== dataset ==")
    print(f"cells={n_cells} timepoints_per_cell={n_time} rows={len(df)}")
    print(f"overall mean distance (nm): {df['distance_nm'].mean():.3f}")
    print(f"overall contact fraction (<= {CONTACT_THRESHOLD_NM:g} nm): {df['contact'].mean():.6f}")
    print(f"overall transcription fraction: {df['transcription'].mean():.6f}")
    print(f"transcription values: {sorted(df['transcription'].unique())}")
    print("== cell metrics summary ==")
    print(metrics.describe().to_string())
    print("== lag analysis (pooled) ==")
    print(lag_df.to_string(index=False))
    print(f"strongest |association|: lag={int(best['lag'])} r={best['association']:.6f} "
          f"n={int(best['n_observations'])}")
    print("== sensitivity: within-cell mean correlation ==")
    print(within_df.to_string(index=False))
    wb = within_df.loc[within_df["mean_within_cell_association"].abs().idxmax()]
    print(f"strongest |within-cell mean|: lag={int(wb['lag'])} r={wb['mean_within_cell_association']:.6f}")
    print(f"lag0 (contemporaneous) pooled r: "
          f"{lag_df.loc[lag_df['lag'] == 0, 'association'].iloc[0]:.6f}")


if __name__ == "__main__":
    main()
