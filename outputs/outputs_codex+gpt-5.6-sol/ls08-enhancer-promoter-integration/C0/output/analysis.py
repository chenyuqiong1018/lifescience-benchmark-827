"""Integrate frozen Hi-C and CRISPR-expression evidence."""

from __future__ import annotations

import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path


OUTPUT_DIR = Path(__file__).resolve().parent
REPO_DIR = Path(__file__).resolve().parents[5]
INPUT_DIR = REPO_DIR / "inputs" / "ls08-enhancer-promoter-integration"
HIC_PATH = INPUT_DIR / "ep.interactions.q1.hic.csv"
EXPR_PATH = INPUT_DIR / "ep.interactions.q1.expr.csv"
ELIGIBILITY_THRESHOLD = 1.645
EXPECTED_PAIRS = [f"EP{number}" for number in range(1, 8)]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def fit_ols(x: list[float], y: list[float]) -> tuple[float, float]:
    x_mean = mean(x)
    y_mean = mean(y)
    denominator = sum((value - x_mean) ** 2 for value in x)
    if denominator == 0:
        raise ValueError("Background distances do not vary")
    slope = sum((xi - x_mean) * (yi - y_mean) for xi, yi in zip(x, y)) / denominator
    intercept = y_mean - slope * x_mean
    return intercept, slope


def main() -> None:
    hic_rows = read_rows(HIC_PATH)
    expr_rows = read_rows(EXPR_PATH)
    background_rows = [row for row in hic_rows if row["set"] == "background"]
    candidate_rows = [row for row in hic_rows if row["set"] == "candidate"]
    if len(background_rows) != 200:
        raise ValueError(f"Expected 200 background rows, observed {len(background_rows)}")
    if [row["pair_id"] for row in candidate_rows] != EXPECTED_PAIRS:
        raise ValueError("Candidate rows must be exactly EP1 through EP7")

    mean_counts: dict[str, float] = {}
    distances: dict[str, float] = {}
    for row in hic_rows:
        values = [float(row[f"count_rep{rep}"]) for rep in (1, 2, 3)]
        mean_counts[row["pair_id"]] = mean(values)
        distances[row["pair_id"]] = float(row["distance_bp"])

    background_x = [math.log10(distances[row["pair_id"]]) for row in background_rows]
    background_y = [math.log1p(mean_counts[row["pair_id"]]) for row in background_rows]
    intercept, slope = fit_ols(background_x, background_y)

    residuals: dict[str, float] = {}
    for row in hic_rows:
        pair_id = row["pair_id"]
        fitted = intercept + slope * math.log10(distances[pair_id])
        residuals[pair_id] = math.log1p(mean_counts[pair_id]) - fitted
    background_residuals = [residuals[row["pair_id"]] for row in background_rows]
    residual_median = statistics.median(background_residuals)
    mad = statistics.median(abs(value - residual_median) for value in background_residuals)
    robust_scale = 1.4826 * mad
    if robust_scale <= 0:
        raise ValueError("Background residual MAD must be positive")
    contact_evidence = {
        pair_id: (residual - residual_median) / robust_scale
        for pair_id, residual in residuals.items()
    }

    expression: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    expr_pairs: set[str] = set()
    for row in expr_rows:
        pair_id = row["pair_id"]
        expr_pairs.add(pair_id)
        expression[(pair_id, row["guide_id"], row["condition"])].append(
            float(row["rna_count"])
        )
    if sorted(expr_pairs) != EXPECTED_PAIRS:
        raise ValueError("Expression pairs must be exactly EP1 through EP7")

    guide_effects: dict[str, list[float]] = defaultdict(list)
    guide_ids = sorted({row["guide_id"] for row in expr_rows})
    for pair_id in EXPECTED_PAIRS:
        for guide_id in guide_ids:
            control = expression[(pair_id, guide_id, "control")]
            perturbed = expression[(pair_id, guide_id, "perturbed")]
            if len(control) != 3 or len(perturbed) != 3:
                raise ValueError(f"Expected three replicates per condition for {pair_id}/{guide_id}")
            effect = math.log2((mean(perturbed) + 0.5) / (mean(control) + 0.5))
            guide_effects[pair_id].append(effect)
    perturbation_effect = {
        pair_id: statistics.median(effects) for pair_id, effects in guide_effects.items()
    }

    records: list[dict[str, float | int | str | bool]] = []
    for pair_id in EXPECTED_PAIRS:
        eligible = contact_evidence[pair_id] >= ELIGIBILITY_THRESHOLD
        support = (
            contact_evidence[pair_id] * abs(perturbation_effect[pair_id]) if eligible else 0.0
        )
        records.append(
            {
                "pair_id": pair_id,
                "contact_evidence": contact_evidence[pair_id],
                "perturbation_effect": perturbation_effect[pair_id],
                "combined_support": support,
                "rank": 0,
                "eligible": eligible,
            }
        )

    eligible_records = sorted(
        (record for record in records if record["eligible"]),
        key=lambda record: (-float(record["combined_support"]), str(record["pair_id"])),
    )
    for rank, record in enumerate(eligible_records, start=1):
        record["rank"] = rank
    if not eligible_records:
        raise ValueError("No physically eligible causal candidates")
    least_supported = max(
        eligible_records, key=lambda record: (int(record["rank"]), str(record["pair_id"]))
    )

    records.sort(key=lambda record: str(record["pair_id"]))
    with (OUTPUT_DIR / "pair_evidence.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pair_id",
                "contact_evidence",
                "perturbation_effect",
                "combined_support",
                "rank",
            ],
        )
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "pair_id": record["pair_id"],
                    "contact_evidence": format(float(record["contact_evidence"]), ".12g"),
                    "perturbation_effect": format(
                        float(record["perturbation_effect"]), ".12g"
                    ),
                    "combined_support": format(float(record["combined_support"]), ".12g"),
                    "rank": record["rank"],
                }
            )

    audit = {
        "least_supported_causal_pair": least_supported["pair_id"],
        "contact_evidence": least_supported["contact_evidence"],
        "perturbation_effect": least_supported["perturbation_effect"],
        "combined_support": least_supported["combined_support"],
        "rank": least_supported["rank"],
        "physical_eligibility_threshold": ELIGIBILITY_THRESHOLD,
        "eligible_pair_count": len(eligible_records),
        "candidate_pair_count": len(records),
        "candidate_pairs": EXPECTED_PAIRS,
        "ep8_invented": False,
        "contact_model": {
            "background_row_count": len(background_rows),
            "response": "log1p(mean_count)",
            "predictor": "log10(distance_bp)",
            "intercept": intercept,
            "slope": slope,
            "background_residual_median": residual_median,
            "background_residual_mad": mad,
            "robust_scale": robust_scale,
        },
        "perturbation_model": {
            "guide_effect": "log2((mean_perturbed + 0.5)/(mean_control + 0.5))",
            "pair_aggregation": "median across guides",
            "guides_per_pair": len(guide_ids),
            "replicates_per_guide_condition": 3,
        },
        "interpretation_limit": (
            "Eligibility and rank follow the frozen local integration policy; they do not "
            "by themselves prove a biological enhancer-promoter causal relationship."
        ),
    }
    with (OUTPUT_DIR / "least_supported.json").open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")

    report = f"""# Enhancer-promoter evidence integration

All supplied rows were integrated by `pair_id`, preserving physical contact and CRISPR-expression evidence as distinct measurements. The Hi-C model used {len(background_rows)} background rows: ordinary least squares of `log1p(mean_count)` on `log10(distance_bp)`. Residuals were standardized with the background residual median and 1.4826 times MAD to produce `contact_evidence`.

For each guide, the perturbation effect was `log2((mean_perturbed + 0.5)/(mean_control + 0.5))`; each pair's `perturbation_effect` is the median across {len(guide_ids)} guides. Pairs with contact evidence below {ELIGIBILITY_THRESHOLD} were ineligible and assigned support and rank zero. Eligible pairs were ranked from strongest to weakest by `contact_evidence × abs(perturbation_effect)`, with `pair_id` as the deterministic tie-breaker.

The least-supported eligible causal candidate under this frozen rule is **{least_supported['pair_id']}** (contact evidence {float(least_supported['contact_evidence']):.4f}, perturbation effect {float(least_supported['perturbation_effect']):.4f}, combined support {float(least_supported['combined_support']):.4f}, rank {least_supported['rank']}). The files contain exactly EP1-EP7; no EP8 row was invented.

“Causal candidate” is a label defined by the supplied integration rule. Physical proximity and expression response are complementary evidence, but this calculation alone does not establish a biological causal enhancer-promoter link or exclude indirect perturbation effects.
"""
    (OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
