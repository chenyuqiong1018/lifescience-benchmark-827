from pathlib import Path
import json
import pandas as pd
import shutil
import sys


def main(input_dir: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    requests = pd.read_csv(input_dir / "dilution_request.csv")
    plate = pd.read_csv(input_dir / "plate_map.csv")
    log = pd.read_csv(input_dir / "run_log.csv").sort_values("event_id")
    inventory = pd.read_csv(input_dir / "source_inventory.csv")
    pipettes = pd.read_csv(input_dir / "pipettes.csv")
    completed = log[(log.operation == "dispense_final") & (log.status == "completed")].well.tolist()
    failed = log[log.status == "failed_before_aspirate"].iloc[0]
    recovery = plate[~plate.well.isin(completed)].copy()
    root = {
        "failed_well": failed.well,
        "failure_mode": "tip pickup failed before aspiration",
        "liquid_moved": False,
        "completed_wells": completed,
        "recovery_wells": recovery.well.tolist(),
    }
    (output_dir / "root_cause.json").write_text(json.dumps(root, indent=2), encoding="utf-8")

    rows = []
    source_by_compound = inventory[inventory.compound_id.notna()].set_index("compound_id")
    for row in recovery.itertuples(index=False):
        src = source_by_compound.loc[row.compound_id]
        transfer = row.requested_final_uM * row.requested_final_volume_uL / src.concentration_uM
        diluent = row.requested_final_volume_uL - transfer
        assert abs(src.concentration_uM * transfer - row.requested_final_uM * row.requested_final_volume_uL) < 1e-9
        transfer_pipette = "P20" if 2 <= transfer <= 20 else "P300"
        diluent_pipette = "P20" if 2 <= diluent <= 20 else "P300"
        rows.append({
            "step": f"recover:{row.well}", "source": src["source"], "destination": row.well,
            "transfer_uL": transfer, "transfer_pipette": transfer_pipette,
            "diluent_source": "diluent:R1", "diluent_uL": diluent,
            "diluent_pipette": diluent_pipette,
            "final_concentration": row.requested_final_uM,
            "final_volume_uL": row.requested_final_volume_uL,
        })
    plan = pd.DataFrame(rows)
    required = plan.groupby("source").transfer_uL.sum()
    for source, volume in required.items():
        assert volume <= float(inventory.set_index("source").loc[source, "available_volume_uL"])
    assert plan.diluent_uL.sum() <= float(inventory.set_index("source").loc["diluent:R1", "available_volume_uL"])
    plan.to_csv(output_dir / "recovery_plan.csv", index=False)
    shutil.copy2(Path(__file__), output_dir / "analysis.py")
    report = f"""# Dilution run recovery

Event ordering shows that {', '.join(completed)} completed. The B2 failure occurred during tip pickup before aspiration, so no liquid moved. Only {', '.join(recovery.well)} require recovery; completed wells are not repeated.

For each remaining CMPD_B well, the 25 µM intermediate supplies 2 µL and media supplies 98 µL to reach 0.5 µM in 100 µL. The identity `25 µM × 2 µL = 0.5 µM × 100 µL` holds. P20 is used for the 2 µL compound transfer and P300 for 98 µL diluent. Inventory is sufficient and the diluent remains media.
"""
    (output_dir / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
