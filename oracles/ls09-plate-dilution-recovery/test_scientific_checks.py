import csv
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("checker", HERE / "scientific_checks.py")
checker = importlib.util.module_from_spec(spec); spec.loader.exec_module(checker)


def make_reference(root: Path, wrong=False, semantic_phrase=False, ambiguous_combined=False):
    out = root / "output"; out.mkdir(parents=True, exist_ok=True)
    root_cause = {"failed_well":"B2", "failure_mode":"Tip pickup failed before aspiration" if semantic_phrase else "tip_pickup_failed_before_aspirate", "liquid_moved":False,
                  "completed_wells":["A1","A2","A3","B1"], "recovery_wells":["B2","B3"]}
    (out / "root_cause.json").write_text(json.dumps(root_cause), encoding="utf-8")
    fields = (["step","source","destination","transfer_uL","diluent_uL","final_concentration","final_volume_uL","pipette"]
              if ambiguous_combined else ["step","source","destination","transfer_uL","transfer_pipette","diluent_source","diluent_uL","diluent_pipette","final_concentration","final_volume_uL"])
    with (out / "recovery_plan.csv").open("w", encoding="utf-8", newline="") as h:
        w=csv.DictWriter(h, fieldnames=fields); w.writeheader()
        for i, well in enumerate(("B2","B3"), 1):
            row={"step":i,"source":"source:A2","destination":well,"transfer_uL":3 if wrong and well=="B2" else 2,
                 "diluent_uL":98,"final_concentration":0.5,"final_volume_uL":100}
            if ambiguous_combined:
                row["pipette"]="P20+P300"
            else:
                row.update({"transfer_pipette":"P20","diluent_source":"diluent:R1","diluent_pipette":"P300"})
            w.writerow(row)
    (out / "report.md").write_text("Tip pickup failed before aspirate at B2. Recover B2 and B3 at 0.5 uM in 100 uL.\n", encoding="utf-8")
    (out / "analysis.py").write_text("# Frozen reference test artifact; submission code is never imported by the oracle.\n", encoding="utf-8")


def test_acceptance_three_runs_and_negative_controls():
    reference = HERE / "testdata" / "reference"
    wrong = HERE / "testdata" / "wrong"
    empty = HERE / "testdata" / "empty"
    semantic = HERE / "testdata" / "semantic_regression"
    ambiguous = HERE / "testdata" / "ambiguous_combined_pipette"
    for _ in range(3):
        make_reference(reference); result=checker.check(reference)
        assert result["hardgate_pass"] and result["core_science"] == 40 and result["direction"] == 15 and result["summary"] == 5
    make_reference(wrong, wrong=True)
    assert not checker.check(wrong)["hardgate_pass"]
    assert not checker.check(empty)["hardgate_pass"]
    make_reference(semantic, semantic_phrase=True)
    semantic_result = checker.check(semantic)
    assert semantic_result["hardgate_pass"] and semantic_result["core_science"] == 40
    make_reference(ambiguous, ambiguous_combined=True)
    ambiguous_result = checker.check(ambiguous)
    assert not ambiguous_result["hardgate_pass"]
    assert "LS09_RECOVERY_TWO_WELL_RECOVERY_CONTRACT" in ambiguous_result["failure_codes"]


if __name__ == "__main__":
    test_acceptance_three_runs_and_negative_controls()
    print("PASS: 3/3 reference, semantic regression, empty/wrong, and ambiguous-pipette control")
