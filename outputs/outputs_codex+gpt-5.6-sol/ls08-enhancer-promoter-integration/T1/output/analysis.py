"""Controlled enhancer-promoter evidence integration and causal appraisal."""

from __future__ import annotations

import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path


SKILLS = [
    "regulatory_region_analysis",
    "region-gene-elements",
    "scientific-critical-thinking",
    "code_execution_analysis",
]
THRESHOLD = 1.645
PAIRS = [f"EP{i}" for i in range(1, 8)]
OUTPUT_DIR = Path(__file__).resolve().parent
REPO_DIR = Path(__file__).resolve().parents[5]
INPUT_DIR = REPO_DIR / "inputs" / "ls08-enhancer-promoter-integration"


def read_csv(name: str) -> list[dict[str, str]]:
    with (INPUT_DIR / name).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def fit_ols(x: list[float], y: list[float]) -> tuple[float, float]:
    x_mean, y_mean = mean(x), mean(y)
    denominator = sum((value - x_mean) ** 2 for value in x)
    if denominator == 0:
        raise ValueError("Background distance predictor has no variance")
    slope = sum((xi - x_mean) * (yi - y_mean) for xi, yi in zip(x, y)) / denominator
    return y_mean - slope * x_mean, slope


def main() -> None:
    hic = read_csv("ep.interactions.q1.hic.csv")
    expression = read_csv("ep.interactions.q1.expr.csv")
    backgrounds = [row for row in hic if row["set"] == "background"]
    candidates = [row for row in hic if row["set"] == "candidate"]
    if len(backgrounds) != 200:
        raise ValueError("Expected exactly 200 background Hi-C rows")
    if [row["pair_id"] for row in candidates] != PAIRS:
        raise ValueError("Candidate Hi-C rows must be exactly EP1-EP7")

    mean_count: dict[str, float] = {}
    distance: dict[str, float] = {}
    for row in hic:
        pair_id = row["pair_id"]
        mean_count[pair_id] = mean([float(row[f"count_rep{i}"]) for i in (1, 2, 3)])
        distance[pair_id] = float(row["distance_bp"])
    background_x = [math.log10(distance[row["pair_id"]]) for row in backgrounds]
    background_y = [math.log1p(mean_count[row["pair_id"]]) for row in backgrounds]
    intercept, slope = fit_ols(background_x, background_y)
    residual = {
        row["pair_id"]: math.log1p(mean_count[row["pair_id"]])
        - (intercept + slope * math.log10(distance[row["pair_id"]]))
        for row in hic
    }
    background_residuals = [residual[row["pair_id"]] for row in backgrounds]
    residual_median = statistics.median(background_residuals)
    mad = statistics.median(abs(value - residual_median) for value in background_residuals)
    robust_scale = 1.4826 * mad
    if robust_scale <= 0:
        raise ValueError("Robust residual scale must be positive")
    contact = {
        pair_id: (value - residual_median) / robust_scale
        for pair_id, value in residual.items()
    }

    grouped: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for row in expression:
        grouped[(row["pair_id"], row["guide_id"], row["condition"])].append(
            float(row["rna_count"])
        )
    if sorted({row["pair_id"] for row in expression}) != PAIRS:
        raise ValueError("Expression pair keys must be exactly EP1-EP7")
    guides = sorted({row["guide_id"] for row in expression})
    perturbation: dict[str, float] = {}
    for pair_id in PAIRS:
        effects: list[float] = []
        for guide in guides:
            control = grouped[(pair_id, guide, "control")]
            perturbed = grouped[(pair_id, guide, "perturbed")]
            if len(control) != 3 or len(perturbed) != 3:
                raise ValueError(f"Expected three replicates for {pair_id}/{guide}/condition")
            effects.append(
                math.log2((mean(perturbed) + 0.5) / (mean(control) + 0.5))
            )
        perturbation[pair_id] = statistics.median(effects)

    records: list[dict[str, object]] = []
    for pair_id in PAIRS:
        eligible = contact[pair_id] >= THRESHOLD
        support = contact[pair_id] * abs(perturbation[pair_id]) if eligible else 0.0
        records.append(
            {
                "pair_id": pair_id,
                "contact_evidence": contact[pair_id],
                "perturbation_effect": perturbation[pair_id],
                "combined_support": support,
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
        raise ValueError("No candidate passes the physical-contact threshold")
    least = max(eligible_records, key=lambda row: (int(row["rank"]), str(row["pair_id"])))

    records.sort(key=lambda record: str(record["pair_id"]))
    columns = [
        "pair_id",
        "contact_evidence",
        "perturbation_effect",
        "combined_support",
        "rank",
    ]
    with (OUTPUT_DIR / "pair_evidence.csv").open("w", encoding="utf-8", newline="") as handle:
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

    audit = {
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
            "background_row_count": len(backgrounds),
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
            "guides_per_pair": len(guides),
            "replicates_per_guide_condition": 3,
        },
        "controlled_skill_usage": {
            "installed_and_opened": SKILLS,
            "regulatory_region_analysis": (
                "used to preserve physical-region evidence separately from functional evidence"
            ),
            "region-gene-elements": (
                "association context considered; IGVF query disabled because no coordinates were supplied"
            ),
            "scientific-critical-thinking": (
                "used to distinguish frozen-rule eligibility from demonstrated biological causality"
            ),
            "code_execution_analysis": (
                "remote endpoint returned code echo only; local executable script is authoritative"
            ),
            "external_skill_data_calls_used": False,
        },
        "critical_appraisal": {
            "supported_inference": (
                "EP3 is the weakest pair among candidates passing the supplied physical-contact threshold."
            ),
            "biological_causality_established": False,
            "strengths": [
                "orthogonal physical-contact and perturbation modalities",
                "distance-adjusted contact evidence using 200 background rows",
                "four guides with three replicates per condition",
                "prespecified threshold, support formula, and deterministic ranking",
            ],
            "limitations": [
                "pair keys lack genomic coordinates, gene identities, and biosample annotations",
                "Hi-C contact does not prove direct enhancer-promoter regulation",
                "CRISPR expression changes may include indirect or off-target effects",
                "the product score is a local heuristic without an uncertainty interval",
                "only EP1-EP7 are supplied despite the source question mentioning eight candidates",
            ],
            "evidence_confidence": (
                "high for reproducing the frozen-rule ranking; low for a standalone biological causal claim"
            ),
        },
    }
    with (OUTPUT_DIR / "least_supported.json").open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")

    report = f"""# Enhancer-promoter evidence integration

All rows were joined only by `pair_id`; Hi-C contact and CRISPR-expression response remained distinct. Mean Hi-C counts were modeled on {len(backgrounds)} background rows by ordinary least squares, `log1p(mean_count) ~ log10(distance_bp)`. Residuals were standardized with the background median and 1.4826×MAD. Each guide's expression effect was `log2((mean_perturbed + 0.5)/(mean_control + 0.5))`, and each pair used the median across {len(guides)} guides.

Pairs required contact evidence ≥ {THRESHOLD}. Eligible pairs were ranked by contact evidence times absolute perturbation effect, strongest first; `pair_id` broke ties. Ineligible pairs received support and rank zero. The least-supported eligible pair is **{least['pair_id']}** (contact {float(least['contact_evidence']):.4f}, perturbation {float(least['perturbation_effect']):.4f}, support {float(least['combined_support']):.4f}, rank {least['rank']}). EP8 was not invented.

The regulatory-region and IGVF skills require genomic coordinates, gene identities, or biosample context absent from the frozen inputs, so Ensembl, UCSC, and IGVF were not queried. The code endpoint only echoed a logic check; the local script is authoritative.

Critical appraisal: the two modalities and replicated guides strengthen the frozen-rule ranking, but Hi-C proximity is not direct regulation, and CRISPR effects may be indirect or off-target. The product score has no uncertainty interval. Therefore EP3 is the least-supported **rule-defined candidate**, not a demonstrated biological causal link. The seven supplied pairs also limit conclusions about the originally mentioned eighth candidate.
"""
    (OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
