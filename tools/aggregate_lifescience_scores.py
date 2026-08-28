import csv
import json
import re
from pathlib import Path

DET = Path(r"C:\Users\chenyuqiong\Documents\Codex\2026-08-27\lifescience-benchmark-827-determin-scorer\outputs\outputs_codex+gpt-5.6-sol")
RUN1 = Path(r"C:\Users\chenyuqiong\Documents\Codex\2026-08-27\lifescience-benchmark-827-judger\outputs")
RUN2 = Path(r"C:\Users\chenyuqiong\Documents\Codex\2026-08-28\lifescience-benchmark-827-judger-run2\outputs\judging_progress.md")
RUN3 = Path(r"C:\Users\chenyuqiong\Documents\Codex\2026-08-28\lifescience-benchmark-827-judger-run3\outputs\blind_judge_run3_scores.csv")
OUT = Path(__file__).resolve().parents[1] / "evaluation-results" / "outputs_codex+gpt-5.6-sol"
ARMS = ("C0", "T0", "T1")
DIMS = ("Evidence", "Method", "Restraint", "Readability")


def load_deterministic():
    records = {}
    bundle = json.loads((DET / "deterministic_scores_tasks_2_4.json").read_text(encoding="utf-8"))
    docs = []
    for path in DET.glob("*/deterministic_scores.json"):
        docs.append(json.loads(path.read_text(encoding="utf-8")))
    for task_no, task in bundle["tasks"].items():
        docs.append({"task_number": int(task_no), "task_id": task["task_id"], "arms": task["arms"]})
    for doc in docs:
        n, slug = int(doc["task_number"]), doc["task_id"]
        for arm, data in doc["arms"].items():
            score = data.get("deterministic_score_original", data.get("deterministic_score"))
            records[(slug, arm)] = {
                "TaskNumber": n,
                "Task": slug,
                "Arm": arm,
                "DeterministicScore": int(score),
                "CoreArtifactMissing": bool(data.get("core_artifact_missing", False)),
                "HardGatePass": bool(data["hardgate_pass"]),
                "FailureCodes": ";".join(data.get("failure_codes", [])),
            }
    return records


def load_run1():
    out = {}
    paths = list(RUN1.glob("task*/run1_results.json")) + list(RUN1.glob("task*/results.json"))
    for path in paths:
        doc = json.loads(path.read_text(encoding="utf-8"))
        for arm, a in doc["arms"].items():
            scores = a.get("scores", a)
            lower = {str(k).lower(): v for k, v in scores.items()}
            vals = {d: float(lower[d.lower()]) for d in DIMS}
            vals["JudgeScore"] = float(a.get("JudgeScore", a.get("judge_score", sum(vals.values()))))
            out[(doc["task_id"], arm)] = vals
    return out


def load_run2():
    text = RUN2.read_text(encoding="utf-8")
    out = {}
    task_matches = list(re.finditer(r"(?m)^##\s+\d+\.\s+([^\s]+)\s*$", text))
    for i, tm in enumerate(task_matches):
        slug = tm.group(1)
        block = text[tm.end(): task_matches[i + 1].start() if i + 1 < len(task_matches) else len(text)]
        arm_matches = list(re.finditer(r"(?m)^###\s+(C0|T0|T1)\s*$", block))
        for j, am in enumerate(arm_matches):
            arm = am.group(1)
            ab = block[am.end(): arm_matches[j + 1].start() if j + 1 < len(arm_matches) else len(block)]
            vals = {}
            for d in DIMS:
                m = re.search(rf"(?m)^\|\s*{d}\s*\|\s*([035])\s*\|", ab)
                if not m:
                    raise ValueError(f"run2 missing {slug} {arm} {d}")
                vals[d] = float(m.group(1))
            m = re.search(r"\*\*JudgeScore:\s*([0-9.]+)/20\.\*\*", ab)
            vals["JudgeScore"] = float(m.group(1))
            out[(slug, arm)] = vals
    return out


def load_run3():
    out = {}
    with RUN3.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            out[(row["task_id"], row["arm"])] = {
                "Evidence": float(row["evidence"]), "Method": float(row["method"]),
                "Restraint": float(row["restraint"]), "Readability": float(row["readability"]),
                "JudgeScore": float(row["judge_score"]),
            }
    return out


def main():
    det, r1, r2, r3 = load_deterministic(), load_run1(), load_run2(), load_run3()
    expected = set(det)
    assert len(expected) == 75, len(expected)
    assert set(r1) == expected, (len(r1), len(expected - set(r1)), len(set(r1) - expected))
    assert set(r2) == expected, (len(r2), len(expected - set(r2)), len(set(r2) - expected))
    assert set(r3) == expected, (len(r3), len(expected - set(r3)), len(set(r3) - expected))
    rows = []
    for key, base in sorted(det.items(), key=lambda kv: (kv[1]["TaskNumber"], ARMS.index(kv[1]["Arm"]))):
        runs = [r1[key], r2[key], r3[key]]
        avg = {d: round(sum(x[d] for x in runs) / 3, 2) for d in DIMS}
        judge = round(sum(x["JudgeScore"] for x in runs) / 3, 2)
        det_score = 0 if base["CoreArtifactMissing"] else base["DeterministicScore"]
        raw = round(det_score + judge, 2)
        display = f"{raw:.2f}" if base["HardGatePass"] else ("FAIL（49）" if raw > 49 else f"FAIL（{raw:.2f}）")
        row = dict(base)
        row.update({
            "Run1_E": runs[0]["Evidence"], "Run1_M": runs[0]["Method"], "Run1_R": runs[0]["Restraint"], "Run1_Read": runs[0]["Readability"], "Run1_Total": runs[0]["JudgeScore"],
            "Run2_E": runs[1]["Evidence"], "Run2_M": runs[1]["Method"], "Run2_R": runs[1]["Restraint"], "Run2_Read": runs[1]["Readability"], "Run2_Total": runs[1]["JudgeScore"],
            "Run3_E": runs[2]["Evidence"], "Run3_M": runs[2]["Method"], "Run3_R": runs[2]["Restraint"], "Run3_Read": runs[2]["Readability"], "Run3_Total": runs[2]["JudgeScore"],
            "AvgEvidence": avg["Evidence"], "AvgMethod": avg["Method"], "AvgRestraint": avg["Restraint"], "AvgReadability": avg["Readability"],
            "JudgeScore": judge, "RawTotal": raw, "DisplayResult": display,
        })
        rows.append(row)
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "scores.csv").open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    (OUT / "scores.json").write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    summary = ["# Codex + gpt-5.6-sol-high 生命科学评测汇总", "", f"- 题目：25", f"- 实验臂：{len(rows)}", f"- Hard gate PASS：{sum(r['HardGatePass'] for r in rows)}", f"- Hard gate FAIL：{sum(not r['HardGatePass'] for r in rows)}", ""]
    (OUT / "README.md").write_text("\n".join(summary), encoding="utf-8")
    print(json.dumps({"rows": len(rows), "pass": sum(r["HardGatePass"] for r in rows), "fail": sum(not r["HardGatePass"] for r in rows), "out": str(OUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
