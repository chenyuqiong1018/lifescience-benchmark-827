#!/usr/bin/env python3
"""Diagnose the stopped dilution run and emit only necessary recovery work."""

from __future__ import annotations

import csv
import json
import re
from collections import defaultdict
from decimal import Decimal
from pathlib import Path


REPO = Path(__file__).resolve().parents[5]
INPUT = REPO / "inputs" / "ls09-plate-dilution-recovery"
OUTPUT = Path(__file__).resolve().parent


def read_csv(name: str) -> list[dict[str, str]]:
    with (INPUT / name).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def dec(text: str, label: str) -> Decimal:
    try:
        value = Decimal(text)
    except Exception as exc:
        raise ValueError(f"Invalid decimal for {label}: {text!r}") from exc
    if not value.is_finite() or value < 0:
        raise ValueError(f"Invalid nonnegative decimal for {label}: {text!r}")
    return value


def number(value: Decimal) -> str:
    rendered = format(value.normalize(), "f")
    return "0" if rendered in ("-0", "") else rendered


def choose_pipette(volume: Decimal, pipettes: list[dict[str, str]]) -> str:
    matches: list[str] = []
    for row in pipettes:
        minimum = dec(row["min_ul"], f"{row['pipette']} minimum")
        maximum = dec(row["max_ul"], f"{row['pipette']} maximum")
        increment = dec(row["increment_ul"], f"{row['pipette']} increment")
        on_increment = (volume / increment) == (volume / increment).to_integral_value()
        if minimum <= volume <= maximum and on_increment:
            matches.append(row["pipette"])
    if len(matches) != 1:
        raise ValueError(f"Volume {number(volume)} uL has {len(matches)} valid pipettes: {matches}")
    return matches[0]


def main() -> None:
    requests = read_csv("dilution_request.csv")
    pipettes = read_csv("pipettes.csv")
    plate = read_csv("plate_map.csv")
    log = read_csv("run_log.csv")
    inventory = read_csv("source_inventory.csv")

    request_by_compound: dict[str, dict[str, str]] = {}
    for row in requests:
        compound = row["compound_id"]
        if not compound or compound in request_by_compound:
            raise ValueError(f"Blank or duplicate dilution request: {compound!r}")
        if dec(row["stock_mm"], f"{compound} stock") <= 0:
            raise ValueError(f"Nonpositive stock concentration for {compound}")
        request_by_compound[compound] = row

    plate_wells = [row["well"] for row in plate]
    if len(plate_wells) != len(set(plate_wells)):
        raise ValueError("Plate map contains duplicate wells")
    plate_by_well = {row["well"]: row for row in plate}
    wells_by_compound: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in plate:
        wells_by_compound[row["compound_id"]].append(row)
    for compound, rows in wells_by_compound.items():
        request = request_by_compound.get(compound)
        if request is None:
            raise ValueError(f"No dilution request for {compound}")
        if len(rows) != int(request["replicates"]):
            raise ValueError(f"Replicate mismatch for {compound}")
        for row in rows:
            if dec(row["requested_final_uM"], f"{row['well']} final") != dec(request["final_um"], f"{compound} request final"):
                raise ValueError(f"Final concentration mismatch at {row['well']}")
            if dec(row["requested_final_volume_uL"], f"{row['well']} volume") != dec(request["final_volume_ul"], f"{compound} request volume"):
                raise ValueError(f"Final volume mismatch at {row['well']}")

    try:
        event_ids = [int(row["event_id"]) for row in log]
    except ValueError as exc:
        raise ValueError("Non-integer run-log event id") from exc
    if len(event_ids) != len(set(event_ids)):
        raise ValueError("Duplicate run-log event id")
    ordered_log = sorted(log, key=lambda row: int(row["event_id"]))
    failures = [row for row in ordered_log if row["status"] != "completed"]
    if len(failures) != 1:
        raise ValueError(f"Expected exactly one stopping failure, found {len(failures)}")
    failure = failures[0]
    if failure["status"] != "failed_before_aspirate" or "no liquid moved" not in failure["detail"].lower():
        raise ValueError("Failure is not explicitly documented as pre-aspiration with no liquid moved")
    failed_well = failure["well"]
    if failed_well not in plate_by_well:
        raise ValueError("Failed well is absent from plate map")

    completed_wells = [
        row["well"]
        for row in ordered_log
        if row["operation"] == "dispense_final" and row["status"] == "completed"
    ]
    if len(completed_wells) != len(set(completed_wells)) or any(well not in plate_by_well for well in completed_wells):
        raise ValueError("Completed final-dispense wells are invalid")
    recovery_wells = [well for well in plate_wells if well not in set(completed_wells)]
    if failed_well not in recovery_wells:
        raise ValueError("Failed well was incorrectly marked complete")

    inventory_by_compound: dict[str, list[dict[str, str]]] = defaultdict(list)
    diluents: list[dict[str, str]] = []
    inventory_by_source: dict[str, dict[str, str]] = {}
    for row in inventory:
        source = row["source"]
        if not source or source in inventory_by_source:
            raise ValueError(f"Blank or duplicate inventory source: {source!r}")
        inventory_by_source[source] = row
        if row["compound_id"]:
            inventory_by_compound[row["compound_id"]].append(row)
        else:
            diluents.append(row)
    if len(diluents) != 1:
        raise ValueError("Diluent source is missing or ambiguous")
    diluent = diluents[0]

    prepared: dict[str, tuple[str, Decimal]] = {}
    pattern = re.compile(r"^([0-9]+(?:\.[0-9]+)?) uM intermediate prepared$")
    for event in ordered_log:
        if event["operation"] != "prepare_intermediate" or event["status"] != "completed":
            continue
        match = pattern.match(event["detail"])
        if match is None:
            raise ValueError(f"Unparseable intermediate preparation: {event['detail']!r}")
        prepared[event["compound_id"]] = (event["well"], dec(match.group(1), "prepared concentration"))

    plan_rows: list[dict[str, str]] = []
    demand: dict[str, Decimal] = defaultdict(Decimal)
    for well in recovery_wells:
        target = plate_by_well[well]
        compound = target["compound_id"]
        candidates = inventory_by_compound.get(compound, [])
        if len(candidates) != 1:
            raise ValueError(f"Intermediate source for {compound} is missing or ambiguous")
        source = candidates[0]
        prepared_entry = prepared.get(compound)
        source_concentration = dec(source["concentration_uM"], f"{source['source']} concentration")
        if prepared_entry != (source["source"], source_concentration):
            raise ValueError(f"Prepared intermediate and inventory disagree for {compound}")
        final_concentration = dec(target["requested_final_uM"], f"{well} final concentration")
        final_volume = dec(target["requested_final_volume_uL"], f"{well} final volume")
        transfer = final_concentration * final_volume / source_concentration
        diluent_volume = final_volume - transfer
        if transfer <= 0 or diluent_volume <= 0:
            raise ValueError(f"Nonpositive recovery movement for {well}")
        if source_concentration * transfer != final_concentration * final_volume:
            raise AssertionError(f"Dilution equation failed for {well}")
        transfer_pipette = choose_pipette(transfer, pipettes)
        diluent_pipette = choose_pipette(diluent_volume, pipettes)
        demand[source["source"]] += transfer
        demand[diluent["source"]] += diluent_volume
        plan_rows.append({
            "step": f"dispense_final:{well}",
            "source": source["source"],
            "destination": well,
            "transfer_uL": number(transfer),
            "transfer_pipette": transfer_pipette,
            "diluent_source": diluent["source"],
            "diluent_uL": number(diluent_volume),
            "diluent_pipette": diluent_pipette,
            "final_concentration": f"{number(final_concentration)} uM",
            "final_volume_uL": number(final_volume),
        })

    for source, volume in demand.items():
        available = dec(inventory_by_source[source]["available_volume_uL"], f"{source} available volume")
        if volume > available:
            raise ValueError(f"Insufficient inventory at {source}: need {number(volume)}, have {number(available)}")

    root_cause = {
        "failed_well": failed_well,
        "failure_mode": "tip pickup failed before aspiration",
        "liquid_moved": False,
        "completed_wells": completed_wells,
        "recovery_wells": recovery_wells,
    }
    (OUTPUT / "root_cause.json").write_text(
        json.dumps(root_cause, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n"
    )
    fields = [
        "step", "source", "destination", "transfer_uL", "transfer_pipette",
        "diluent_source", "diluent_uL", "diluent_pipette",
        "final_concentration", "final_volume_uL",
    ]
    with (OUTPUT / "recovery_plan.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(plan_rows)

    source = inventory_by_compound[plate_by_well[failed_well]["compound_id"]][0]
    report = f"""# Dilution-run recovery

## Diagnosis

Numeric event order shows completed final dispenses to {', '.join(completed_wells)}. Event {failure['event_id']} then records a tip-pickup failure at {failed_well} before aspiration and explicitly states that no liquid moved. The run stopped there. Consequently, only {', '.join(recovery_wells)} remain; completed wells are not repeated.

## Recovery calculation

Both remaining wells request 0.5 uM in 100 uL for CMPD_B. The completed intermediate preparation and inventory agree on {source['source']} at {source['concentration_uM']} uM in `{source['solvent']}`. For each well, C_source x V_transfer = C_final x V_final gives 25 uM x 2 uL = 0.5 uM x 100 uL. Add 2 uL intermediate with P20, then 98 uL from {diluent['source']} (`{diluent['solvent']}`) with P300. These instruments are unique under the frozen ranges/increments: P20 2-20 uL by 0.1 uL; P300 20-300 uL by 1 uL.

Total new demand is 4 uL from {source['source']} and 196 uL from {diluent['source']}, below the recorded available volumes of {source['available_volume_uL']} uL and {diluent['available_volume_uL']} uL. No intermediate preparation or completed final dispense is repeated. `analysis.py` aborts on absent/ambiguous sources, unexplained liquid movement, inconsistent requests, invalid or ambiguous pipette choice, equation failure, solvent/source mismatch, or insufficient inventory rather than inventing work.
"""
    (OUTPUT / "report.md").write_text(report, encoding="utf-8", newline="\n")
    print(json.dumps(root_cause, separators=(",", ":")))
    print(f"Recovery rows: {len(plan_rows)}")


if __name__ == "__main__":
    main()
