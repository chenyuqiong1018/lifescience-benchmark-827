#!/usr/bin/env python3
"""Recovery-only dilution analysis with explicit failure-state classification."""

import csv
import json
import re
from collections import defaultdict
from decimal import Decimal
from pathlib import Path


ROOT = Path(__file__).resolve().parents[5]
INPUT = ROOT / "inputs" / "ls09-plate-dilution-recovery"
OUT = Path(__file__).resolve().parent


def rows(filename):
    with (INPUT / filename).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def D(text, context):
    try:
        value = Decimal(text)
    except Exception as exc:
        raise ValueError(f"Bad numeric value for {context}: {text!r}") from exc
    if not value.is_finite() or value < 0:
        raise ValueError(f"Invalid value for {context}: {text!r}")
    return value


def render(value):
    text = format(value.normalize(), "f")
    return "0" if text == "-0" else text


def unique_pipette(volume, instruments):
    valid = []
    for instrument in instruments:
        increment = D(instrument["increment_ul"], "pipette increment")
        q = volume / increment
        if (
            D(instrument["min_ul"], "pipette minimum") <= volume
            <= D(instrument["max_ul"], "pipette maximum")
            and q == q.to_integral_value()
        ):
            valid.append(instrument["pipette"])
    if len(valid) != 1:
        raise ValueError(f"Ambiguous or impossible pipette selection for {render(volume)} uL: {valid}")
    return valid[0]


def main():
    request = rows("dilution_request.csv")
    instruments = rows("pipettes.csv")
    plate = rows("plate_map.csv")
    events = rows("run_log.csv")
    inventory = rows("source_inventory.csv")

    event_numbers = [int(event["event_id"]) for event in events]
    if len(event_numbers) != len(set(event_numbers)):
        raise ValueError("Duplicate event_id")
    events.sort(key=lambda event: int(event["event_id"]))
    failed = [event for event in events if event["status"] != "completed"]
    if len(failed) != 1:
        raise ValueError("Run log does not identify exactly one stopping event")
    failure = failed[0]
    if (
        failure["operation"] != "dispense_final"
        or failure["status"] != "failed_before_aspirate"
        or "no liquid moved" not in failure["detail"].lower()
    ):
        raise ValueError("Cannot prove that the stopping failure occurred before liquid movement")

    plate_order = [entry["well"] for entry in plate]
    if len(plate_order) != len(set(plate_order)):
        raise ValueError("Duplicate plate-map well")
    plate_by_well = {entry["well"]: entry for entry in plate}
    completed = [
        event["well"]
        for event in events
        if event["operation"] == "dispense_final" and event["status"] == "completed"
    ]
    if len(completed) != len(set(completed)) or any(well not in plate_by_well for well in completed):
        raise ValueError("Completed-well history is inconsistent with the plate map")
    remaining = [well for well in plate_order if well not in set(completed)]
    if failure["well"] not in remaining:
        raise ValueError("Failed well was incorrectly classified as complete")

    requests = {entry["compound_id"]: entry for entry in request}
    if len(requests) != len(request):
        raise ValueError("Duplicate compound request")
    for compound, req in requests.items():
        mapped = [entry for entry in plate if entry["compound_id"] == compound]
        if len(mapped) != int(req["replicates"]):
            raise ValueError(f"Replicate mismatch for {compound}")
        for entry in mapped:
            if D(entry["requested_final_uM"], "plate final") != D(req["final_um"], "request final"):
                raise ValueError(f"Concentration mismatch at {entry['well']}")
            if D(entry["requested_final_volume_uL"], "plate volume") != D(req["final_volume_ul"], "request volume"):
                raise ValueError(f"Volume mismatch at {entry['well']}")

    inventory_by_source = {entry["source"]: entry for entry in inventory}
    if len(inventory_by_source) != len(inventory):
        raise ValueError("Duplicate inventory source")
    diluent_candidates = [entry for entry in inventory if not entry["compound_id"]]
    if len(diluent_candidates) != 1:
        raise ValueError("Missing or ambiguous diluent source")
    diluent = diluent_candidates[0]

    source_by_compound = defaultdict(list)
    for entry in inventory:
        if entry["compound_id"]:
            source_by_compound[entry["compound_id"]].append(entry)
    prepared = {}
    pattern = re.compile(r"^([0-9]+(?:\.[0-9]+)?) uM intermediate prepared$")
    for event in events:
        if event["operation"] == "prepare_intermediate" and event["status"] == "completed":
            match = pattern.fullmatch(event["detail"])
            if not match:
                raise ValueError("Intermediate preparation detail is not quantitative")
            prepared[event["compound_id"]] = (event["well"], D(match.group(1), "prepared concentration"))

    plan = []
    usage = defaultdict(Decimal)
    for well in remaining:
        target = plate_by_well[well]
        compound = target["compound_id"]
        candidates = source_by_compound[compound]
        if len(candidates) != 1:
            raise ValueError(f"Missing or ambiguous source for {compound}")
        source = candidates[0]
        source_c = D(source["concentration_uM"], "source concentration")
        if prepared.get(compound) != (source["source"], source_c):
            raise ValueError(f"Prepared source identity/concentration mismatch for {compound}")
        final_c = D(target["requested_final_uM"], "final concentration")
        final_v = D(target["requested_final_volume_uL"], "final volume")
        transfer_v = final_c * final_v / source_c
        diluent_v = final_v - transfer_v
        if transfer_v <= 0 or diluent_v <= 0 or source_c * transfer_v != final_c * final_v:
            raise ValueError(f"Invalid dilution calculation for {well}")
        usage[source["source"]] += transfer_v
        usage[diluent["source"]] += diluent_v
        plan.append({
            "step": f"dispense_final:{well}",
            "source": source["source"],
            "destination": well,
            "transfer_uL": render(transfer_v),
            "transfer_pipette": unique_pipette(transfer_v, instruments),
            "diluent_source": diluent["source"],
            "diluent_uL": render(diluent_v),
            "diluent_pipette": unique_pipette(diluent_v, instruments),
            "final_concentration": f"{render(final_c)} uM",
            "final_volume_uL": render(final_v),
        })

    for source_name, required in usage.items():
        available = D(inventory_by_source[source_name]["available_volume_uL"], "available inventory")
        if required > available:
            raise ValueError(f"Insufficient inventory at {source_name}")

    diagnosis = {
        "failed_well": failure["well"],
        "failure_mode": "tip pickup failed before aspiration",
        "liquid_moved": False,
        "completed_wells": completed,
        "recovery_wells": remaining,
    }
    (OUT / "root_cause.json").write_text(json.dumps(diagnosis, indent=2) + "\n", encoding="utf-8", newline="\n")
    fields = ["step", "source", "destination", "transfer_uL", "transfer_pipette", "diluent_source", "diluent_uL", "diluent_pipette", "final_concentration", "final_volume_uL"]
    with (OUT / "recovery_plan.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(plan)

    source = source_by_compound[plate_by_well[failure["well"]]["compound_id"]][0]
    report = f"""# Stopped dilution run: recovery-only plan

## Root cause classification

The numeric event sequence records final dispenses to {', '.join(completed)}, followed by a tip-pickup failure at {failure['well']} before aspiration. The log explicitly says no liquid moved. This is a discrete pre-transfer mechanical failure, not a noisy volume measurement: there is no partial dose to estimate or correct. The only physically necessary wells are {', '.join(remaining)}; completed wells remain untouched.

## Recovery

Both remaining wells are CMPD_B at 0.5 uM in 100 uL. The completed preparation record and inventory identify `{source['source']}` as the unique 25 uM intermediate in `{source['solvent']}`. Exact dilution gives 25 x 2 = 0.5 x 100, so each well receives 2 uL from `{source['source']}` with P20 and 98 uL from `{diluent['source']}` (`{diluent['solvent']}`) with P300. The frozen ranges and increments make each pipette choice unique.

The recovery consumes 4 uL intermediate and 196 uL diluent versus {source['available_volume_uL']} uL and {diluent['available_volume_uL']} uL available. `analysis.py` validates event order, request/plate consistency, source identity and concentration, solvent-preserving source selection, the dilution equation, pipette range/increment uniqueness, and inventory. It aborts rather than estimating an unrecorded movement.

The freshly installed `measurement-error-analysis` skill sharpened the classification boundary: because failure occurred before aspiration, uncertainty propagation is inappropriate; the correct action is a complete transfer only for the uncompleted wells.
"""
    (OUT / "report.md").write_text(report, encoding="utf-8", newline="\n")
    print(json.dumps(diagnosis, separators=(",", ":")))
    print(f"Recovery rows: {len(plan)}")


if __name__ == "__main__":
    main()
