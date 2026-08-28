#!/usr/bin/env python3
"""Reproducible independent-groups Cohen-d and power calculation."""

from __future__ import annotations

import csv
import json
from math import ceil, sqrt
from pathlib import Path

from statsmodels.stats.power import tt_ind_solve_power


OUTPUT_DIR = Path(__file__).resolve().parent
INPUT_DIR = OUTPUT_DIR.parents[4] / "inputs" / "ls10-neun-power-analysis"


def stop(reason: str) -> None:
    raise SystemExit(f"ABORT: {reason}")


def sample_summary(numbers: list[float]) -> tuple[int, float, float]:
    n = len(numbers)
    if n < 2:
        stop("each group must contain at least two observed values")
    mean = sum(numbers) / n
    sd = sqrt(sum((value - mean) ** 2 for value in numbers) / (n - 1))
    return n, mean, sd


def main() -> None:
    contract = (INPUT_DIR / "README.md").read_text(encoding="utf-8")
    if "missing observations remain missing" not in contract:
        stop("missing-data rule is absent from the input contract")

    groups: dict[str, list[float]] = {}
    labels: list[str] = []
    with (INPUT_DIR / "NeuN_quantification.csv").open(
        newline="", encoding="utf-8-sig"
    ) as handle:
        for row in csv.DictReader(handle):
            label = row["Hemispere"].strip()
            if not label:
                stop("encountered an empty group label")
            if label not in groups:
                labels.append(label)
                groups[label] = []
            measurement = row["NeuN"].strip()
            if measurement == "":
                continue
            try:
                groups[label].append(float(measurement))
            except ValueError:
                stop(f"non-numeric NeuN measurement: {measurement!r}")

    if len(labels) != 2:
        stop(f"independent two-group calculation requires 2 labels, found {len(labels)}")
    label_1, label_2 = labels
    n_1, mean_1, sd_1 = sample_summary(groups[label_1])
    n_2, mean_2, sd_2 = sample_summary(groups[label_2])
    pooled_sd = sqrt(
        ((n_1 - 1) * sd_1**2 + (n_2 - 1) * sd_2**2) / (n_1 + n_2 - 2)
    )
    if pooled_sd <= 0:
        stop("pooled standard deviation must be positive")
    d = (mean_1 - mean_2) / pooled_sd
    if d == 0:
        stop("zero observed effect has no finite sample-size solution")

    alpha = 0.05
    target_power = 0.80
    alternative = "two-sided"
    continuous_n = float(
        tt_ind_solve_power(
            effect_size=abs(d),
            nobs1=None,
            alpha=alpha,
            power=target_power,
            ratio=1.0,
            alternative=alternative,
        )
    )
    required_n = ceil(continuous_n)
    power_at_required = float(
        tt_ind_solve_power(
            effect_size=abs(d),
            nobs1=required_n,
            alpha=alpha,
            power=None,
            ratio=1.0,
            alternative=alternative,
        )
    )
    power_one_less = float(
        tt_ind_solve_power(
            effect_size=abs(d),
            nobs1=required_n - 1,
            alpha=alpha,
            power=None,
            ratio=1.0,
            alternative=alternative,
        )
    )
    if not (power_one_less < target_power <= power_at_required):
        stop("integer sample-size boundary did not bracket target power")

    payload = {
        "group_labels": labels,
        "n_each": {label_1: n_1, label_2: n_2},
        "means": {label_1: mean_1, label_2: mean_2},
        "sds": {label_1: sd_1, label_2: sd_2},
        "pooled_sd": pooled_sd,
        "cohens_d": d,
        "alpha": alpha,
        "power": target_power,
        "alternative": alternative,
        "required_n_per_group": required_n,
    }
    (OUTPUT_DIR / "power_result.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n"
    )

    report = f"""# NeuN standardized difference and power

## Data handling

The frozen file contains two labels in first-occurrence order: **{label_1}** and **{label_2}**. Missing NeuN cells would be omitted rather than converted to zero; none are missing, leaving n = {n_1} and n = {n_2}. Although sample identifiers repeat across labels, the requested design is explicitly a two-sided **independent** t-test, so no paired analysis was substituted.

## Effect size

{label_1}: M = {mean_1:.3f}, sample SD = {sd_1:.6f}.  
{label_2}: M = {mean_2:.3f}, sample SD = {sd_2:.6f}.

The equal-variance pooled SD is **{pooled_sd:.6f}**. Cohen's d, defined in the frozen label order as `({label_1} mean - {label_2} mean) / pooled SD`, is **{d:.6f}**. Its sign therefore records direction; the two-sided power calculation uses its absolute magnitude.

## Required sample size

Using the noncentral-t independent-samples model with equal group sizes, alpha = {alpha:.2f}, target power = {target_power:.2f}, and a two-sided alternative, the continuous solution is {continuous_n:.6f} per group. Rounding upward yields **{required_n} per group**. Power is {power_at_required:.6f} at {required_n}, while {required_n - 1} gives {power_one_less:.6f}; therefore {required_n} is the minimum integer meeting 0.80 under the specified model.

The `statistical-analysis` skill guided the independent-group Cohen's d definition, use of sample SDs, statsmodels power solver, transparent direction convention, and upward/minimality check. No inferential group-comparison p-value was requested or added.
"""
    (OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
