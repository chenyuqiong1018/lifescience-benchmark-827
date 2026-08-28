#!/usr/bin/env python3
"""Audit NeuN Cohen's d and solve equal-group t-test power exactly."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from scipy import stats
from statsmodels.stats.power import TTestIndPower


OUT = Path(__file__).resolve().parent
INP = OUT.parents[4] / "inputs" / "ls10-neun-power-analysis"
ALPHA = 0.05
POWER_TARGET = 0.80


def abort(detail: str) -> None:
    raise SystemExit(f"ABORT: {detail}")


def summarize(data: list[float]) -> tuple[float, float]:
    if len(data) < 2:
        abort("a group has fewer than two non-missing observations")
    mean = math.fsum(data) / len(data)
    variance = math.fsum((x - mean) ** 2 for x in data) / (len(data) - 1)
    return mean, math.sqrt(variance)


def independent_t_power(effect: float, n_per_group: int, alpha: float) -> float:
    """Two-sided power from the noncentral-t distribution, equal groups."""
    degrees_freedom = 2 * n_per_group - 2
    critical = stats.t.ppf(1 - alpha / 2, degrees_freedom)
    noncentrality = effect * math.sqrt(n_per_group / 2)
    upper = stats.nct.sf(critical, degrees_freedom, noncentrality)
    lower = stats.nct.cdf(-critical, degrees_freedom, noncentrality)
    return float(upper + lower)


def main() -> None:
    contract = (INP / "README.md").read_text(encoding="utf-8")
    required_contract_text = (
        "Group labels and row order are frozen",
        "sample standard deviations",
        "missing observations remain missing",
    )
    if any(text not in contract for text in required_contract_text):
        abort("README contract is incomplete")

    labels: list[str] = []
    observations: dict[str, list[float]] = {}
    with (INP / "NeuN_quantification.csv").open(
        newline="", encoding="utf-8-sig"
    ) as handle:
        reader = csv.DictReader(handle)
        expected_columns = ["Sample", "Hemispere", "NeuN", "Sex"]
        if reader.fieldnames != expected_columns:
            abort(f"unexpected CSV columns: {reader.fieldnames}")
        for row in reader:
            group = row["Hemispere"].strip()
            if group not in observations:
                labels.append(group)
                observations[group] = []
            raw = row["NeuN"].strip()
            if raw == "":
                continue
            try:
                observations[group].append(float(raw))
            except ValueError:
                abort(f"invalid NeuN measurement {raw!r}")

    if len(labels) != 2 or any(label == "" for label in labels):
        abort(f"expected two non-empty groups, obtained {labels}")
    g1, g2 = labels
    n_each = {g1: len(observations[g1]), g2: len(observations[g2])}
    mean_1, sd_1 = summarize(observations[g1])
    mean_2, sd_2 = summarize(observations[g2])
    means = {g1: mean_1, g2: mean_2}
    sds = {g1: sd_1, g2: sd_2}
    df = n_each[g1] + n_each[g2] - 2
    if df <= 0:
        abort("non-positive pooled degrees of freedom")
    pooled_sd = math.sqrt(
        ((n_each[g1] - 1) * sd_1**2 + (n_each[g2] - 1) * sd_2**2) / df
    )
    if pooled_sd == 0:
        abort("pooled standard deviation is zero")
    cohens_d = (mean_1 - mean_2) / pooled_sd
    effect = abs(cohens_d)
    if effect == 0:
        abort("zero effect cannot reach target power at a finite sample size")

    required_n = None
    for candidate in range(2, 1_000_001):
        if independent_t_power(effect, candidate, ALPHA) >= POWER_TARGET:
            required_n = candidate
            break
    if required_n is None:
        abort("sample-size search exceeded its explicit bound")
    power_at_n = independent_t_power(effect, required_n, ALPHA)
    power_below = independent_t_power(effect, required_n - 1, ALPHA)
    continuous_n = float(
        TTestIndPower().solve_power(
            effect_size=effect,
            nobs1=None,
            alpha=ALPHA,
            power=POWER_TARGET,
            ratio=1.0,
            alternative="two-sided",
        )
    )
    if required_n != math.ceil(continuous_n):
        abort("noncentral-t integer search disagrees with statsmodels continuous solution")

    # Diagnostics describe the observed pilot data but do not replace the
    # independent equal-variance design explicitly requested in the prompt.
    shapiro_1 = stats.shapiro(observations[g1])
    shapiro_2 = stats.shapiro(observations[g2])
    levene = stats.levene(observations[g1], observations[g2], center="median")

    result = {
        "group_labels": labels,
        "n_each": n_each,
        "means": means,
        "sds": sds,
        "pooled_sd": pooled_sd,
        "cohens_d": cohens_d,
        "alpha": ALPHA,
        "power": POWER_TARGET,
        "alternative": "two-sided",
        "required_n_per_group": required_n,
    }
    (OUT / "power_result.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n"
    )

    report = f"""# NeuN power-analysis audit

## Frozen analysis choice

The prompt specifies a two-sided independent t-test, equal group sizes, alpha = {ALPHA:.2f}, and power = {POWER_TARGET:.2f}. The repeated Sample identifiers were therefore not used to silently substitute a paired design. All NeuN cells are observed; no missing value was imputed.

## Descriptive statistics and effect size

- {g1}: n = {n_each[g1]}, mean = {mean_1:.3f}, sample SD = {sd_1:.6f}
- {g2}: n = {n_each[g2]}, mean = {mean_2:.3f}, sample SD = {sd_2:.6f}

The equal-variance pooled SD is **{pooled_sd:.6f}**. In frozen label order ({g1} minus {g2}), Cohen's d is **{cohens_d:.6f}**. The sign records direction; the requested two-sided power calculation uses |d|.

For transparency, pilot-data diagnostics were Shapiro-Wilk p = {shapiro_1.pvalue:.6f} ({g1}) and {shapiro_2.pvalue:.6f} ({g2}), with median-centered Levene p = {levene.pvalue:.6f}. These small-sample diagnostics do not redefine the test model fixed by the prompt.

## Power result

Direct noncentral-t enumeration finds **{required_n} observations per group**. At {required_n - 1}, power is {power_below:.6f}; at {required_n}, it is {power_at_n:.6f}. The statsmodels continuous solution is {continuous_n:.6f}, whose upward rounding independently confirms {required_n}.

## Skill use

`statistical-analysis` supplied the independent-group Cohen's d convention, sample-SD requirement, assumption transparency, and noncentral-t planning model. `code_execution_analysis` was asked for a fixed-data arithmetic cross-check, but its endpoint only returned the code without executing it; the rerunnable local script is therefore the executed source of record. `biomarker_discovery` was opened as required but its external biomarker databases were not queried because they are irrelevant to the supplied measurements and external data are forbidden.
"""
    (OUT / "report.md").write_text(report, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
