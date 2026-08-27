from __future__ import annotations

import csv
import json
from pathlib import Path

ACCEPTED = True


def _failure_mode_is_tip_pickup_before_aspirate(value):
    text = str(value or "").strip().lower().replace("-", " ").replace("_", " ")
    tokens = set(text.split())
    tip_pickup = "tip" in tokens and ("pickup" in tokens or {"pick", "up"}.issubset(tokens))
    failed = bool(tokens & {"failed", "failure", "error"})
    before_aspirate = "before" in tokens and bool(tokens & {"aspirate", "aspiration"})
    return tip_pickup and failed and before_aspirate


def check(workspace: Path):
    output = workspace / "output"
    failures = []
    criteria = {}
    try:
        root = json.loads((output / "root_cause.json").read_text(encoding="utf-8"))
    except Exception:
        root = {}
    root_ok = (
        root.get("failed_well") == "B2"
        and _failure_mode_is_tip_pickup_before_aspirate(root.get("failure_mode"))
        and root.get("liquid_moved") is False
        and set(root.get("completed_wells", [])) == {"A1", "A2", "A3", "B1"}
        and set(root.get("recovery_wells", [])) == {"B2", "B3"}
    )

    try:
        with (output / "recovery_plan.csv").open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except Exception:
        rows = []
    required = {"step", "source", "destination", "transfer_uL", "transfer_pipette", "diluent_source", "diluent_uL", "diluent_pipette", "final_concentration", "final_volume_uL"}
    schema_ok = len(rows) == 2 and required.issubset(rows[0] if rows else {})
    plan_ok = schema_ok
    concentration_ok = schema_ok
    pipette_ok = schema_ok
    no_overdraw = schema_ok
    seen = set()
    used = 0.0
    if schema_ok:
        for row in rows:
            try:
                transfer = float(row["transfer_uL"])
                diluent = float(row["diluent_uL"])
                final_c = float(row["final_concentration"])
                final_v = float(row["final_volume_uL"])
            except (TypeError, ValueError):
                plan_ok = concentration_ok = pipette_ok = no_overdraw = False
                continue
            well = row["destination"]
            seen.add(well)
            used += transfer
            if row["source"] != "source:A2" or row["diluent_source"] != "diluent:R1" or well not in {"B2", "B3"} or abs(transfer - 2) > 1e-6 or abs(diluent - 98) > 1e-6:
                plan_ok = False
            if abs(final_v - 100) > 1e-6 or abs(final_c - 0.5) > 1e-6 or abs(25 * transfer / final_v - final_c) > 1e-6:
                concentration_ok = False
            if row["transfer_pipette"] != "P20" or row["diluent_pipette"] != "P300" or not (2 <= transfer <= 20) or not (20 <= diluent <= 300):
                pipette_ok = False
        plan_ok = plan_ok and seen == {"B2", "B3"}
        no_overdraw = no_overdraw and used <= 1998

    report_text = ""
    try:
        report_text = (output / "report.md").read_text(encoding="utf-8").lower()
    except Exception:
        pass
    summary_ok = all(term in report_text for term in ["b2", "b3", "tip pickup", "0.5", "100"])

    criteria.update({"root_cause_from_log": root_ok, "two_well_recovery_contract": plan_ok,
                     "dilution_mass_balance": concentration_ok, "pipette_feasibility": pipette_ok,
                     "source_not_overdrawn": no_overdraw, "report_consistency": summary_ok})
    for name, ok in criteria.items():
        if not ok:
            failures.append("LS09_RECOVERY_" + name.upper())
    core = (14 if root_ok else 0) + (10 if plan_ok else 0) + (10 if concentration_ok else 0) + (6 if pipette_ok and no_overdraw else 0)
    direction = 15 if root_ok and plan_ok and concentration_ok and pipette_ok and no_overdraw else 0
    summary = 5 if summary_ok else 0
    return {"core_science": core, "direction": direction, "summary": summary,
            "hardgate_pass": root_ok and plan_ok and concentration_ok and pipette_ok and no_overdraw,
            "failure_codes": failures, "criteria": criteria}
