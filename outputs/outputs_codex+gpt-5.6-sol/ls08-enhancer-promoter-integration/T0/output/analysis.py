"""Integrate physical and perturbational regulatory-pair evidence."""

from __future__ import annotations

import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path


SKILL = "regulatory_region_analysis"
THRESHOLD = 1.645
PAIRS = [f"EP{i}" for i in range(1, 8)]
OUTPUT_DIR = Path(__file__).resolve().parent
REPO_DIR = Path(__file__).resolve().parents[5]
INPUT_DIR = REPO_DIR / "inputs" / "ls08-enhancer-promoter-integration"


def load_csv(filename: str) -> list[dict[str, str]]:
    with (INPUT_DIR / filename).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def average(values: list[float]) -> float:
    return sum(values) / len(values)


def ols(x: list[float], y: list[float]) -> tuple[float, float]:
    xbar, ybar = average(x), average(y)
    denominator = sum((value - xbar) ** 2 for value in x)
    if denominator == 0:
        raise ValueError("OLS predictor has no variance")
    slope = sum((xi - xbar) * (yi - ybar) for xi, yi in zip(x, y)) / denominator
    return ybar - slope * xbar, slope


def main() -> None:
    hic = load_csv("ep.interactions.q1.hic.csv")
    expr = load_csv("ep.interactions.q1.expr.csv")
    background = [row for row in hic if row["set"] == "background"]
    candidates = [row for row in hic if row["set"] == "candidate"]
    if len(background) != 200 or [row["pair_id"] for row in candidates] != PAIRS:
        raise ValueError("Hi-C input must contain 200 backgrounds and EP1-EP7 candidates")

    counts: dict[str, float] = {}
    distance: dict[str, float] = {}
    for row in hic:
        pair_id = row["pair_id"]
        counts[pair_id] = average([float(row[f"count_rep{i}"]) for i in (1, 2, 3)])
        distance[pair_id] = float(row["distance_bp"])
    x = [math.log10(distance[row["pair_id"]]) for row in background]
    y = [math.log1p(counts[row["pair_id"]]) for row in background]
    intercept, slope = ols(x, y)
    residual = {
        row["pair_id"]: math.log1p(counts[row["pair_id"]])
        - (intercept + slope * math.log10(distance[row["pair_id"]]))
        for row in hic
    }
    background_residual = [residual[row["pair_id"]] for row in background]
    residual_median = statistics.median(background_residual)
    mad = statistics.median(abs(value - residual_median) for value in background_residual)
    scale = 1.4826 * mad
    if scale <= 0:
        raise ValueError("Robust residual scale must be positive")
    contact = {
        pair_id: (value - residual_median) / scale for pair_id, value in residual.items()
    }

    grouped: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for row in expr:
        grouped[(row["pair_id"], row["guide_id"], row["condition"])].append(
            float(row["rna_count"])
        )
    guides = sorted({row["guide_id"] for row in expr})
    if sorted({row["pair_id"] for row in expr}) != PAIRS:
        raise ValueError("Expression input must contain exactly EP1-EP7")
    pair_effect: dict[str, float] = {}
    for pair_id in PAIRS:
        effects: list[float] = []
        for guide in guides:
            control = grouped[(pair_id, guide, "control")]
            perturbed = grouped[(pair_id, guide, "perturbed")]
            if len(control) != 3 or len(perturbed) != 3:
                raise ValueError(f"Replicate mismatch for {pair_id}/{guide}")
            effects.append(
                math.log2((average(perturbed) + 0.5) / (average(control) + 0.5))
            )
        pair_effect[pair_id] = statistics.median(effects)

    records: list[dict[str, object]] = []
    for pair_id in PAIRS:
        eligible = contact[pair_id] >= THRESHOLD
        records.append(
            {
                "pair_id": pair_id,
                "contact_evidence": contact[pair_id],
                "perturbation_effect": pair_effect[pair_id],
                "combined_support": contact[pair_id] * abs(pair_effect[pair_id])
                if eligible
                else 0.0,
                "rank": 0,
                "eligible": eligible,
            }
        )
    eligible_records = sorted(
        [record for record in records if record["eligible"]],
        key=lambda record: (-float(record["combined_support"]), str(record["pair_id"])),
    )
    for rank, record in enumerate(eligible_records, 1):
        record["rank"] = rank
    if not eligible_records:
        raise ValueError("No physically eligible pair")
    least = max(eligible_records, key=lambda row: (int(row["rank"]), str(row["pair_id"])))

    records.sort(key=lambda row: str(row["pair_id"]))
    with (OUTPUT_DIR / "pair_evidence.csv").open("w", encoding="utf-8", newline="") as handle:
        columns = [
            "pair_id",
            "contact_evidence",
            "perturbation_effect",
            "combined_support",
            "rank",
        ]
        writer = csv.DictWriter(handle, fieldnames=columns)
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

    result = {
        "least_supported_causal_pair": least["pair_id"],
        "contact_evidence": least["contact_evidence"],
        "perturbation_effect": least["perturbation_effect"],
        "combined_support": least["combined_support"],
        "rank": least["rank"],
        "physical_eligibility_threshold": THRESHOLD,
        "eligible_pair_count": len(eligible_records),
        "candidate_pair_count": len(records),
        "candidate_pairs": PAIRS,
        "ep8_invented": False,
        "contact_model": {
            "background_row_count": len(background),
            "response": "log1p(mean_count)",
            "predictor": "log10(distance_bp)",
            "intercept": intercept,
            "slope": slope,
            "background_residual_median": residual_median,
            "background_residual_mad": mad,
            "robust_scale": scale,
        },
        "perturbation_model": {
            "guide_effect": "log2((mean_perturbed + 0.5)/(mean_control + 0.5))",
            "pair_aggregation": "median across guides",
            "guides_per_pair": len(guides),
            "replicates_per_guide_condition": 3,
        },
        "skill_usage": {
            "selected_skill": SKILL,
            "applied_as": "separate physical-region and functional-perturbation evidence",
            "external_ensembl_ucsc_calls_used": False,
            "reason": "The frozen inputs provide pair keys but no genomic coordinates for external lookup.",
        },
        "interpretation_limit": (
            "The frozen-rule causal-candidate label integrates contact and perturbation evidence "
            "but does not establish a biological causal enhancer-promoter link."
        ),
    }
    with (OUTPUT_DIR / "least_supported.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")

    report = f"""# Enhancer-promoter evidence integration

Following the `{SKILL}` evidence-separation perspective, all supplied rows were joined only by `pair_id`. Physical Hi-C contact and CRISPR-expression response remained distinct. No Ensembl/UCSC overlap, sequence, binding-matrix, or phenotype endpoint was called: the frozen pair-level files provide no genomic coordinates and the supplied integration rule does not permit external annotation.

For Hi-C, the mean of three counts was transformed with `log1p`; ordinary least squares against `log10(distance_bp)` used {len(background)} background rows. Residuals were standardized by the background median and 1.4826×MAD. For each guide, expression effect was `log2((mean_perturbed + 0.5)/(mean_control + 0.5))`; pair effects are medians across {len(guides)} guides.

Pairs required contact evidence ≥ {THRESHOLD}. Eligible pairs were ranked from strongest to weakest by contact evidence times absolute perturbation effect; ineligible pairs received support and rank zero. The least-supported eligible pair is **{least['pair_id']}** (contact {float(least['contact_evidence']):.4f}, perturbation {float(least['perturbation_effect']):.4f}, support {float(least['combined_support']):.4f}, rank {least['rank']}). Inputs contain EP1-EP7 only; EP8 was not invented.

The term “causal candidate” follows the frozen local rule. Contact plus expression response is stronger than either modality alone, but it does not itself prove a direct causal enhancer-promoter relationship or rule out indirect perturbation effects.
"""
    (OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
