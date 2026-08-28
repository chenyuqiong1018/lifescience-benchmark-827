#!/usr/bin/env python3
"""Compute Cohen's d and equal-group independent-t-test sample size."""

from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path

from statsmodels.stats.power import TTestIndPower


HERE = Path(__file__).resolve().parent
TASK_INPUT = HERE.parents[4] / "inputs" / "ls10-neun-power-analysis"
ALPHA = 0.05
TARGET_POWER = 0.80
ALTERNATIVE = "two-sided"


def abort(message: str) -> None:
    raise SystemExit(f"ABORT: {message}")


def main() -> None:
    readme = (TASK_INPUT / "README.md").read_text(encoding="utf-8")
    if "sample standard deviations" not in readme:
        abort("input contract does not specify sample standard deviations")

    with (TASK_INPUT / "NeuN_quantification.csv").open(
        newline="", encoding="utf-8-sig"
    ) as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        abort("no observations")

    group_labels: list[str] = []
    values: dict[str, list[float]] = {}
    for row in rows:
        label = row["Hemispere"].strip()
        raw = row["NeuN"].strip()
        if not label:
            abort("missing group label")
        if label not in values:
            group_labels.append(label)
            values[label] = []
        if raw:
            try:
                values[label].append(float(raw))
            except ValueError:
                abort(f"invalid NeuN value: {raw!r}")
    if len(group_labels) != 2:
        abort(f"expected two groups, found {len(group_labels)}")
    if any(len(values[label]) < 2 for label in group_labels):
        abort("each group needs at least two non-missing observations")

    first, second = group_labels
    means = {label: statistics.fmean(values[label]) for label in group_labels}
    sds = {label: statistics.stdev(values[label]) for label in group_labels}
    ns = {label: len(values[label]) for label in group_labels}
    pooled_variance = (
        (ns[first] - 1) * sds[first] ** 2
        + (ns[second] - 1) * sds[second] ** 2
    ) / (ns[first] + ns[second] - 2)
    pooled_sd = math.sqrt(pooled_variance)
    if pooled_sd == 0:
        abort("pooled standard deviation is zero")
    cohens_d = (means[first] - means[second]) / pooled_sd
    if cohens_d == 0:
        abort("observed effect size is zero; finite target sample size is undefined")

    power_model = TTestIndPower()
    solved_n = power_model.solve_power(
        effect_size=abs(cohens_d),
        nobs1=None,
        alpha=ALPHA,
        power=TARGET_POWER,
        ratio=1.0,
        alternative=ALTERNATIVE,
    )
    required_n = math.ceil(solved_n)
    achieved = power_model.power(
        effect_size=abs(cohens_d),
        nobs1=required_n,
        alpha=ALPHA,
        ratio=1.0,
        alternative=ALTERNATIVE,
    )
    previous = power_model.power(
        effect_size=abs(cohens_d),
        nobs1=required_n - 1,
        alpha=ALPHA,
        ratio=1.0,
        alternative=ALTERNATIVE,
    )
    if achieved < TARGET_POWER or previous >= TARGET_POWER:
        abort("rounded sample size failed the minimal-power check")

    result = {
        "group_labels": group_labels,
        "n_each": ns,
        "means": means,
        "sds": sds,
        "pooled_sd": pooled_sd,
        "cohens_d": cohens_d,
        "alpha": ALPHA,
        "power": TARGET_POWER,
        "alternative": ALTERNATIVE,
        "required_n_per_group": required_n,
    }
    (HERE / "power_result.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n"
    )

    report = f"""# NeuN effect-size and power analysis

The two frozen groups are **{first}** and **{second}**, with {ns[first]} and {ns[second]} non-missing observations. Their means are {means[first]:.3f} and {means[second]:.3f}; their sample standard deviations (`ddof=1`) are {sds[first]:.6f} and {sds[second]:.6f}.

Using the usual equal-variance pooled standard deviation,

`s_p = sqrt(((n1-1)s1^2 + (n2-1)s2^2) / (n1+n2-2)) = {pooled_sd:.6f}`.

With the frozen group order {first} minus {second}, Cohen's d is **{cohens_d:.6f}**. The positive sign means the observed {first} mean is higher; power uses `abs(d)` because the requested test is two-sided.

For an equal-size, two-sided independent t-test with alpha {ALPHA:.2f} and target power {TARGET_POWER:.2f}, the continuous solution is {solved_n:.6f} observations per group. Rounding upward gives **{required_n} observations per group**. The resulting power is {achieved:.6f}; {required_n - 1} per group would give {previous:.6f}, so {required_n} is the smallest integer meeting the target under this model.
"""
    (HERE / "report.md").write_text(report, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
