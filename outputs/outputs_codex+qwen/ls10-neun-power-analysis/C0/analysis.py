from pathlib import Path
import json
import math
import shutil
import sys

import pandas as pd
from statsmodels.stats.power import TTestIndPower


def main(input_dir: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(input_dir / "NeuN_quantification.csv")
    label_col = "Hemispere"
    labels = list(dict.fromkeys(df[label_col].astype(str)))
    if len(labels) != 2:
        raise ValueError("Expected exactly two groups")
    groups = [df.loc[df[label_col].astype(str).eq(label), "NeuN"].astype(float) for label in labels]
    ns = [len(x) for x in groups]
    means = [float(x.mean()) for x in groups]
    sds = [float(x.std(ddof=1)) for x in groups]
    pooled = math.sqrt(((ns[0] - 1) * sds[0] ** 2 + (ns[1] - 1) * sds[1] ** 2) / (ns[0] + ns[1] - 2))
    d = (means[0] - means[1]) / pooled
    alpha, power = 0.05, 0.80
    required = math.ceil(TTestIndPower().solve_power(effect_size=abs(d), alpha=alpha, power=power, ratio=1, alternative="two-sided"))
    result = {
        "group_labels": labels,
        "n_each": dict(zip(labels, ns)),
        "means": dict(zip(labels, means)),
        "sds": dict(zip(labels, sds)),
        "pooled_sd": pooled,
        "cohens_d": d,
        "alpha": alpha,
        "power": power,
        "alternative": "two-sided",
        "required_n_per_group": required,
    }
    (output_dir / "power_result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    shutil.copy2(Path(__file__), output_dir / "analysis.py")
    report = f"""# NeuN independent-groups power analysis

The two supplied groups were analyzed as independent samples, as required by the task contract. The standardized mean difference is defined as `{labels[0]} - {labels[1]}` divided by the unbiased pooled within-group standard deviation. A two-sided equal-size independent t-test at alpha {alpha} and power {power} requires {required} observations per group after rounding upward.

- {labels[0]}: n={ns[0]}, mean={means[0]:.6f}, SD={sds[0]:.6f}
- {labels[1]}: n={ns[1]}, mean={means[1]:.6f}, SD={sds[1]:.6f}
- pooled SD: {pooled:.6f}
- Cohen's d ({labels[0]} - {labels[1]}): {d:.6f}
"""
    (output_dir / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
