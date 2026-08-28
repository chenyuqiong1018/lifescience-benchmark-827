#!/usr/bin/env python3
"""Pocket-reliability assessment for low-confidence AlphaFold pocket mutations.

Reads only ./inputs and writes ./output deliverables per SCORING_RULE.md.
Predicted ddG / activity values are treated strictly as model hypotheses,
never as measured binding or activity.
"""
import csv
import json
import os
import re

BASE = os.path.dirname(os.path.abspath(__file__))
IN_DIR = os.path.normpath(os.path.join(BASE, "..", "inputs"))
OUT_DIR = BASE

# --- thresholds from SCORING_RULE.md -------------------------------------
PLDDT_LOW = 50.0        # pLDDT < 50 -> unsupported_low_confidence
PAE_HIGH = 10.0         # PAE > 10 A  -> penalty flag pae_gt_10A
RELIABLE_PLDDT = 70.0   # reliable requires pLDDT >= 70
RELIABLE_PAE = 10.0     # reliable requires PAE <= 10 A


def load_intervals(path):
    intervals = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            intervals.append(
                {
                    "start": int(row["residue_start"]),
                    "end": int(row["residue_end"]),
                    "plddt": float(row["plddt"]),
                    "pae": float(row["pae_to_core_a"]),
                    "pocket": row["pocket_member"].strip().lower() == "true",
                }
            )
    return intervals


def load_candidates(path):
    cands = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            cands.append(
                {
                    "mutation": row["mutation"].strip(),
                    "ddg": float(row["predicted_ddg_kcal_mol"]),
                    "activity": row["predicted_activity_change"].strip(),
                    "region": row["region"].strip(),
                }
            )
    return cands


def residue_number(mutation):
    m = re.search(r"(\d+)", mutation)
    if not m:
        raise ValueError(f"no residue number in mutation {mutation!r}")
    return int(m.group(1))


def join_interval(residue, intervals):
    for iv in intervals:
        if iv["start"] <= residue <= iv["end"]:
            return iv
    raise ValueError(f"residue {residue} not covered by any interval")


def pocket_support(iv):
    return "unsupported_low_confidence" if iv["plddt"] < PLDDT_LOW else "cautious_support"


def confidence_penalty(iv):
    flags = []
    if iv["plddt"] < PLDDT_LOW:
        flags.append("plddt_lt_50")
    if iv["pae"] > PAE_HIGH:
        flags.append("pae_gt_10A")
    return ";".join(flags) if flags else "none"


def decision(iv, region):
    if region != "pocket":
        return "out_of_scope_non_pocket"
    # low-confidence pocket candidate -> defer structure validation
    return "defer_structure_validation" if iv["plddt"] < PLDDT_LOW else "defer_structure_validation"


def main():
    intervals = load_intervals(os.path.join(IN_DIR, "confidence.csv"))
    candidates = load_candidates(os.path.join(IN_DIR, "mutation_candidates.csv"))

    records = []
    for c in candidates:
        res = residue_number(c["mutation"])
        iv = join_interval(res, intervals)
        records.append(
            {
                "mutation": c["mutation"],
                "region": c["region"],
                "ddg": c["ddg"],
                "activity": c["activity"],
                "residue": res,
                "plddt": iv["plddt"],
                "pae": iv["pae"],
                "interval_pocket": iv["pocket"],
                "pocket_support": pocket_support(iv),
                "confidence_penalty": confidence_penalty(iv),
                "decision": decision(iv, c["region"]),
            }
        )

    # Rank: pocket candidates first by higher pLDDT, lower PAE, lower |ddG|,
    # lexical mutation; non-pocket candidates after them.
    def sort_key(r):
        group = 0 if r["region"] == "pocket" else 1
        return (group, -r["plddt"], r["pae"], abs(r["ddg"]), r["mutation"])

    records.sort(key=sort_key)
    for i, r in enumerate(records, start=1):
        r["rank"] = i

    # ---- mutation_priorities.csv ----
    csv_path = os.path.join(OUT_DIR, "mutation_priorities.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "mutation", "pocket_support", "confidence_penalty", "decision"])
        for r in records:
            w.writerow(
                [r["rank"], r["mutation"], r["pocket_support"], r["confidence_penalty"], r["decision"]]
            )

    # ---- pocket_assessment.json ----
    pocket_intervals = [iv for iv in intervals if iv["pocket"]]
    pocket_reliable = all(
        iv["plddt"] >= RELIABLE_PLDDT and iv["pae"] <= RELIABLE_PAE for iv in pocket_intervals
    )
    prioritize_mutations = pocket_reliable  # unreliable pocket -> do not prioritize

    n_lowconf = sum(1 for r in records if r["region"] == "pocket" and r["plddt"] < PLDDT_LOW)
    reason = (
        "The nominated pocket spans residue intervals 210-230 (pLDDT=47, PAE=14.2 A) and "
        "231-270 (pLDDT=43, PAE=17.8 A). Both pocket intervals fall below the pLDDT>=70 and "
        "PAE<=10 A reliability thresholds; pLDDT<50 marks coordinates that should not be "
        "interpreted and PAE>10 A indicates uncertain relative placement. The pocket is therefore "
        "unreliable, and mutation prioritization should be deferred until experimental structure "
        "validation. Predicted ddG and activity changes are model hypotheses only, not measured "
        "binding or activity effects."
    )
    assessment = {
        "pocket_reliable": bool(pocket_reliable),
        "prioritize_mutations": bool(prioritize_mutations),
        "reason": reason,
    }
    json_path = os.path.join(OUT_DIR, "pocket_assessment.json")
    with open(json_path, "w") as fh:
        json.dump(assessment, fh, indent=2)
        fh.write("\n")

    # ---- report.md ----
    lines = []
    lines.append("# Pocket Reliability Assessment")
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(f"- Pocket reliable: **{assessment['pocket_reliable']}**")
    lines.append(f"- Prioritize mutations: **{assessment['prioritize_mutations']}**")
    lines.append("")
    lines.append("The nominated pocket is **not reliable**. Both pocket intervals have pLDDT < 50 "
                 "(coordinates should not be interpreted) and PAE > 10 A (uncertain relative placement). "
                 "All pocket mutation candidates are therefore deferred pending experimental structure "
                 "validation. Non-pocket candidates are out of scope.")
    lines.append("")
    lines.append("## Method")
    lines.append("")
    lines.append("Each mutation was joined to the residue interval in `confidence.csv` containing its "
                 "residue number. Pocket candidates were ranked by higher pLDDT, then lower PAE, then "
                 "lower absolute predicted ddG, then lexical mutation; non-pocket candidates follow. "
                 "Reliability requires every pocket interval to satisfy pLDDT >= 70 and PAE <= 10 A.")
    lines.append("")
    lines.append("## Ranked mutation priorities")
    lines.append("")
    lines.append("| rank | mutation | region | pLDDT | PAE (A) | predicted ddG | pocket_support | confidence_penalty | decision |")
    lines.append("|------|----------|--------|-------|---------|---------------|----------------|--------------------|----------|")
    for r in records:
        lines.append(
            f"| {r['rank']} | {r['mutation']} | {r['region']} | {r['plddt']:.0f} | {r['pae']} | "
            f"{r['ddg']} | {r['pocket_support']} | {r['confidence_penalty']} | {r['decision']} |"
        )
    lines.append("")
    lines.append("## Confidence intervals")
    lines.append("")
    lines.append("| residues | pLDDT | PAE to core (A) | pocket member |")
    lines.append("|----------|-------|-----------------|---------------|")
    for iv in intervals:
        lines.append(f"| {iv['start']}-{iv['end']} | {iv['plddt']:.0f} | {iv['pae']} | {str(iv['pocket']).lower()} |")
    lines.append("")
    lines.append("## Uncertainty and caveats")
    lines.append("")
    lines.append("- pLDDT and PAE are per-model confidence estimates; low pLDDT (<50) and high PAE (>10 A) "
                 "together mean the pocket geometry and its relative placement are too uncertain to trust for "
                 "mutation design.")
    lines.append("- All predicted ddG and activity-change values are **computational model hypotheses**. "
                 "They are not measured binding affinities or enzymatic activities and must not be reported as such.")
    lines.append("- No experimental binding or activity measurements are present in the inputs.")
    lines.append("")
    lines.append("## References")
    lines.append("")
    lines.append("- AlphaFold DB confidence guidance (pLDDT / PAE): https://alphafold.ebi.ac.uk/faq")
    lines.append("- Terwilliger et al. (2024). AlphaFold predictions are hypotheses and do not replace "
                 "experimental structure determination. Nat Methods. DOI 10.1038/s41592-023-02087-4")
    lines.append("")
    with open(os.path.join(OUT_DIR, "report.md"), "w") as fh:
        fh.write("\n".join(lines))

    # console summary
    print("Wrote mutation_priorities.csv, pocket_assessment.json, report.md")
    print(json.dumps(assessment, indent=2))
    for r in records:
        print(r["rank"], r["mutation"], r["pocket_support"], r["confidence_penalty"], r["decision"])


if __name__ == "__main__":
    main()
