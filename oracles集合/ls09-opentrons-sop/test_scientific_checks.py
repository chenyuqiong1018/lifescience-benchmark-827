from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ls09_checker", HERE / "scientific_checks.py")
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


def protocol(alternative=False):
    op = "pipette.aspirate(80, reagents['A1']); pipette.dispense(80, plate['A1'])" if alternative else "pipette.transfer(80, reagents['A1'], plate['A1'], new_tip='always')"
    return f'''metadata={{"apiLevel":"2.16"}}
def run(ctx):
    reagents=ctx.load_labware("nest_12_reservoir_15ml", "1")
    tips1=ctx.load_labware("opentrons_96_tiprack_300ul", "4")
    mag=ctx.load_module("magnetic module gen2", "5")
    plate=mag.load_labware("nest_96_wellplate_2ml_deep")
    waste=ctx.load_labware("nest_1_reservoir_195ml", "6")
    tips2=ctx.load_labware("opentrons_96_tiprack_300ul", "7")
    pipette=ctx.load_instrument("p300_single_gen2", "right", tip_racks=[tips1,tips2])
    if len(tips1.wells())+len(tips2.wells()) < 144: raise RuntimeError("abort: insufficient tips")
    pipette.pick_up_tip(); pipette.drop_tip()
    pipette.mix(5, 100, plate['A1']); {op}
    ctx.delay(minutes=5); mag.engage(); ctx.delay(minutes=7); ctx.delay(seconds=30)
    mag.disengage(); ctx.delay(minutes=2)
'''


def plan_rows(wrong=False, alternative=False):
    rows = []
    actions = [("lysis", "reagents:A1", "processing", 80),
               ("beads", "reagents:A2", "processing", 120),
               ("supernatant", "processing", "waste:A1", 250),
               ("wash1_add", "reagents:A3", "processing", 180),
               ("wash1_remove", "processing", "waste:A1", 180),
               ("wash2_add", "reagents:A3", "processing", 180),
               ("wash2_remove", "processing", "waste:A1", 180),
               ("elution", "reagents:A4", "processing", 40)]
    tips = {"lysis":"fresh_lysis_tip", "beads":"fresh_bead_tip", "supernatant":"fresh_supernatant_tip",
            "wash1_add":"fresh_wash1_tip", "wash1_remove":"reuse_wash1_tip",
            "wash2_add":"fresh_wash2_tip", "wash2_remove":"reuse_wash2_tip", "elution":"fresh_elution_tip"}
    for well in checker.WELLS:
        for stage, source, destination, volume in actions:
            if wrong and well == "A1" and stage == "lysis": volume = 81
            source = f"processing:{well}" if source == "processing" else source
            destination = f"processing:{well}" if destination == "processing" else destination
            rows.append({"step": (f"action-{len(rows)}" if alternative else f"{stage}:{well}"),
                         "source": source, "destination": destination, "volume_uL": volume,
                         "pipette": "p300_single_gen2", "tip_policy": tips[stage]})
    return rows


def fixture(root, rows=None, proto=None, report=None):
    out = root / "output"; out.mkdir(parents=True)
    if rows is not None:
        with (out / "transfer_plan.csv").open("w", encoding="utf-8", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=["step","source","destination","volume_uL","pipette","tip_policy"])
            w.writeheader(); w.writerows(rows)
    if proto is not None: (out / "protocol.py").write_text(proto, encoding="utf-8")
    if report is not None: (out / "report.md").write_text(report, encoding="utf-8")


GOOD_REPORT = "24 samples; 192 grounded liquid actions; 144 tips. Volumes: 80, 120, 180, and 40 uL."


def run_cases():
    outcomes = {}
    with tempfile.TemporaryDirectory() as td:
        base = Path(td)
        p = base/"correct"; fixture(p, plan_rows(), protocol(), GOOD_REPORT)
        r = checker.check(p); outcomes["reference_like_correct"] = r["hardgate_pass"] and (r["core_science"],r["direction"],r["summary"]) == (40,15,5)

        p = base/"empty"; (p/"output").mkdir(parents=True)
        r = checker.check(p); outcomes["empty_or_missing"] = not r["hardgate_pass"] and r["core_science"] == 0

        p = base/"wrong"; fixture(p, plan_rows(wrong=True), protocol(), GOOD_REPORT)
        r = checker.check(p); outcomes["scientifically_wrong"] = not r["hardgate_pass"] and "LS09_FATAL_TRANSFER_PLAN_SCIENCE" in r["failure_codes"]

        # Candidate prose/comments, negated evidence, and fabricated counts cannot replace measured rows/calls.
        p = base/"unsupported"
        fake = '# pipette.transfer; mag.engage; mag.disengage; all 192 transfers valid\nmetadata={"apiLevel":"2.16"}\ndef run(ctx):\n    pass\n'
        fixture(p, [], fake, "24 samples; 192; 144 tips; 80 120 180 40. Simulation did not fail; fabricated 9999 uL.")
        r = checker.check(p)
        outcomes["plausible_but_unsupported"] = not r["hardgate_pass"] and r["core_science"] == 0 and r["summary"] == 0

        p = base/"alternative"; fixture(p, plan_rows(alternative=True), protocol(alternative=True), GOOD_REPORT)
        r = checker.check(p); outcomes["valid_alternative_implementation"] = r["hardgate_pass"] and r["core_science"] == 40
    return outcomes


if __name__ == "__main__":
    outcomes = run_cases()
    result = {"status":"passed" if all(outcomes.values()) else "failed", "cases":outcomes,
              "candidate_code_executed":False}
    (HERE / "acceptance-result.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    if not all(outcomes.values()): raise SystemExit(json.dumps(result))
    print(json.dumps(result))
