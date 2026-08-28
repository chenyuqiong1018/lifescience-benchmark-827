#!/usr/bin/env python3
"""Reconstruct and validate only the physically necessary dilution recovery."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from decimal import Decimal
from pathlib import Path


HERE = Path(__file__).resolve().parent
TASK_INPUT = HERE.parents[4] / "inputs" / "ls09-plate-dilution-recovery"


def fail(message: str) -> None:
    raise SystemExit(f"ABORT: {message}")


def read_csv(name: str) -> list[dict[str, str]]:
    path = TASK_INPUT / name
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def as_decimal(value: str, label: str) -> Decimal:
    try:
        return Decimal(value)
    except Exception:
        fail(f"invalid decimal for {label}: {value!r}")


def volume_is_usable(volume: Decimal, pipette: dict[str, str]) -> bool:
    minimum = as_decimal(pipette["min_ul"], "pipette minimum")
    maximum = as_decimal(pipette["max_ul"], "pipette maximum")
    increment = as_decimal(pipette["increment_ul"], "pipette increment")
    if not minimum <= volume <= maximum or increment <= 0:
        return False
    return volume % increment == 0


def unique_pipette(volume: Decimal, pipettes: list[dict[str, str]]) -> str:
    choices = [p["pipette"] for p in pipettes if volume_is_usable(volume, p)]
    if len(choices) != 1:
        fail(f"{volume} uL has {len(choices)} frozen-range pipette choices: {choices}")
    return choices[0]


def main() -> None:
    # README is part of the input contract and is intentionally read as well.
    readme = (TASK_INPUT / "README.md").read_text(encoding="utf-8")
    if "numeric `event_id` order" not in readme:
        fail("input contract does not define event ordering")

    requests = read_csv("dilution_request.csv")
    plate = read_csv("plate_map.csv")
    events = read_csv("run_log.csv")
    inventory = read_csv("source_inventory.csv")
    pipettes = read_csv("pipettes.csv")

    request_by_compound = {row["compound_id"]: row for row in requests}
    if len(request_by_compound) != len(requests):
        fail("duplicate compound request")
    plate_counts = Counter(row["compound_id"] for row in plate)
    for compound, request in request_by_compound.items():
        if plate_counts[compound] != int(request["replicates"]):
            fail(f"plate map replicate count differs for {compound}")

    try:
        events = sorted(events, key=lambda row: int(row["event_id"]))
    except (KeyError, ValueError):
        fail("event_id must be numeric")
    if len({row["event_id"] for row in events}) != len(events):
        fail("duplicate event_id")

    dispense = [row for row in events if row["operation"] == "dispense_final"]
    completed_wells = [row["well"] for row in dispense if row["status"] == "completed"]
    if len(set(completed_wells)) != len(completed_wells):
        fail("a well has multiple completed dispense events")

    failures = [row for row in dispense if row["status"] != "completed"]
    if len(failures) != 1:
        fail(f"expected one stopped dispense event, found {len(failures)}")
    failure = failures[0]
    before_aspiration = failure["status"] == "failed_before_aspirate"
    no_liquid_moved = "no liquid moved" in failure["detail"].lower()
    if not (before_aspiration and no_liquid_moved):
        fail("failed dispense does not prove that no liquid moved before aspiration")

    requested_wells = [row["well"] for row in plate]
    unknown_completed = set(completed_wells) - set(requested_wells)
    if unknown_completed:
        fail(f"completed wells absent from plate map: {sorted(unknown_completed)}")
    recovery_wells = [well for well in requested_wells if well not in set(completed_wells)]
    if failure["well"] not in recovery_wells:
        fail("failed well is not recoverable without redoing a completed well")

    plate_by_well = {row["well"]: row for row in plate}
    compound_sources: dict[str, list[dict[str, str]]] = {}
    diluents: list[dict[str, str]] = []
    for row in inventory:
        if row["compound_id"]:
            compound_sources.setdefault(row["compound_id"], []).append(row)
        elif as_decimal(row["concentration_uM"], "diluent concentration") == 0:
            diluents.append(row)
    if len(diluents) != 1:
        fail(f"expected one explicit zero-concentration diluent, found {len(diluents)}")
    diluent = diluents[0]
    if diluent["solvent"] != "media":
        fail("recorded diluent solvent identity is not media")

    plan: list[dict[str, str]] = []
    demand: Counter[str] = Counter()
    for step, well in enumerate(recovery_wells, start=1):
        target = plate_by_well[well]
        compound = target["compound_id"]
        if compound not in request_by_compound:
            fail(f"no dilution request for {compound}")
        request = request_by_compound[compound]
        target_c = as_decimal(target["requested_final_uM"], f"{well} final concentration")
        final_v = as_decimal(target["requested_final_volume_uL"], f"{well} final volume")
        if target_c != as_decimal(request["final_um"], f"{compound} request concentration"):
            fail(f"plate/request concentration mismatch for {well}")
        if final_v != as_decimal(request["final_volume_ul"], f"{compound} request volume"):
            fail(f"plate/request volume mismatch for {well}")

        sources = compound_sources.get(compound, [])
        if len(sources) != 1:
            fail(f"expected one recorded intermediate source for {compound}, found {len(sources)}")
        source = sources[0]
        if source["solvent"] != "DMSO_0.5pct_in_media":
            fail(f"unexpected source solvent identity for {source['source']}")
        source_c = as_decimal(source["concentration_uM"], f"{source['source']} concentration")
        if source_c <= target_c:
            fail(f"source concentration cannot produce requested dilution for {well}")
        transfer_v = target_c * final_v / source_c
        diluent_v = final_v - transfer_v
        if source_c * transfer_v != target_c * final_v:
            fail(f"mass-balance equation failed for {well}")
        if transfer_v <= 0 or diluent_v <= 0:
            fail(f"missing or non-positive transfer for {well}")

        transfer_pipette = unique_pipette(transfer_v, pipettes)
        diluent_pipette = unique_pipette(diluent_v, pipettes)
        demand[source["source"]] += transfer_v
        demand[diluent["source"]] += diluent_v
        plan.append(
            {
                "step": str(step),
                "source": source["source"],
                "destination": well,
                "transfer_uL": format(transfer_v, "f"),
                "transfer_pipette": transfer_pipette,
                "diluent_source": diluent["source"],
                "diluent_uL": format(diluent_v, "f"),
                "diluent_pipette": diluent_pipette,
                "final_concentration": f"{format(target_c, 'f')} uM",
                "final_volume_uL": format(final_v, "f"),
            }
        )

    inventory_by_source = {row["source"]: row for row in inventory}
    for source_name, required in demand.items():
        available = as_decimal(inventory_by_source[source_name]["available_volume_uL"], "available volume")
        if required > available:
            fail(f"insufficient inventory at {source_name}: need {required} uL, have {available} uL")

    root_cause = {
        "failed_well": failure["well"],
        "failure_mode": "tip pickup failed before aspiration",
        "liquid_moved": False,
        "completed_wells": completed_wells,
        "recovery_wells": recovery_wells,
    }
    (HERE / "root_cause.json").write_text(
        json.dumps(root_cause, indent=2) + "\n", encoding="utf-8", newline="\n"
    )

    fieldnames = [
        "step",
        "source",
        "destination",
        "transfer_uL",
        "transfer_pipette",
        "diluent_source",
        "diluent_uL",
        "diluent_pipette",
        "final_concentration",
        "final_volume_uL",
    ]
    with (HERE / "recovery_plan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(plan)

    report = f"""# Dilution-run recovery report

## Root cause

The run stopped at well **{failure['well']}** because tip pickup failed before aspiration. The event log explicitly states that no liquid moved, so the failed well remains physically unfilled. This is a categorical mechanical failure, not evidence of a transferred-volume measurement error.

Completed final-dispense wells, in numeric event order: **{', '.join(completed_wells)}**. They are excluded from recovery. The only wells still requiring work are **{', '.join(recovery_wells)}**.

## Recovery

For each of {', '.join(recovery_wells)}, transfer **2 uL** from **source:A2** with **P20**, then add **98 uL** from **diluent:R1** with **P300**. Each well finishes at **100 uL of 0.5 uM CMPD_B**. No completed well is repeated.

The concentration check is `25 uM * 2 uL = 0.5 uM * 100 uL = 50 uM*uL`. Frozen ranges select exactly one physical instrument for each movement: P20 for 2 uL and P300 for 98 uL. Total remaining demand is **4 uL** from source:A2 and **196 uL** from diluent:R1, below the recorded current inventories of 1998 uL and 5000 uL. The recorded source solvent (`DMSO_0.5pct_in_media`) and diluent solvent (`media`) are preserved without substitution.

## Skill-assisted audit

The controlled protocol-generation result agreed with the two fixed liquid movements and introduced no extra operation. The code-execution capability returned the requested exact-arithmetic snippet but did not execute it, so this rerunnable script independently performs the Decimal calculation and all validations. Measurement-error guidance was applied by treating the pre-aspiration tip-pickup failure as a discrete run-state event rather than inventing an uncertainty correction. Bioassay database lookups and unrelated unit-conversion endpoints were not used because the task forbids external data and already fixes concentration and volume units.
"""
    (HERE / "report.md").write_text(report, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
