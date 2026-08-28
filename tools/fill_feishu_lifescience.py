import csv
import json
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

DOC = "https://aicarrier.feishu.cn/wiki/OMGgwu5UyinGImkh0dEcnVEinre"
TABLE_ID = "doxcn6h5ns2o8Bm29C2d6LriSsW"
LARK = [r"C:\Users\chenyuqiong\AppData\Roaming\npm\lark-cli.cmd"]
ROOT = Path(__file__).resolve().parents[1]
SCORES = ROOT / "evaluation-results" / "outputs_codex+gpt-5.6-sol" / "scores.csv"


def run(args, stdin=None):
    p = subprocess.run(args, input=stdin, text=True, encoding="utf-8", capture_output=True)
    if p.returncode:
        raise RuntimeError(p.stderr or p.stdout)
    return p.stdout


def fmt(v):
    return f"{float(v):.2f}"


def main():
    by_task = {}
    with SCORES.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            by_task.setdefault(row["Task"], {})[row["Arm"]] = row
    assert len(by_task) == 25 and all(set(v) == {"C0", "T0", "T1"} for v in by_task.values())

    fetched = json.loads(run(LARK + [
        "docs", "+fetch", "--doc", DOC, "--scope", "range",
        "--start-block-id", TABLE_ID, "--end-block-id", TABLE_ID, "--detail", "with-ids"
    ]))
    fragment = ET.fromstring(fetched["data"]["document"]["content"])
    table = fragment.find("table")
    assert table is not None and table.attrib.get("id") == TABLE_ID
    updated = 0
    for tr in table.findall("./tbody/tr"):
        cells = tr.findall("td")
        task_idx = None
        task = None
        for i, td in enumerate(cells):
            text = "".join(td.itertext()).strip()
            if text in by_task:
                task_idx, task = i, text
                break
        if task is None:
            continue
        targets = cells[task_idx + 1:]
        assert len(targets) >= 14, (task, len(targets))
        r = by_task[task]
        c0, t0, t1 = r["C0"], r["T0"], r["T1"]
        values = [
            fmt(c0["DeterministicScore"]), fmt(c0["JudgeScore"]), fmt(c0["RawTotal"]),
            fmt(t0["DeterministicScore"]), fmt(t0["JudgeScore"]), fmt(t0["RawTotal"]),
            fmt(t1["DeterministicScore"]), fmt(t1["JudgeScore"]), fmt(t1["RawTotal"]),
            "", fmt(float(t0["RawTotal"]) - float(c0["RawTotal"])),
            fmt(float(t1["RawTotal"]) - float(c0["RawTotal"])),
            "三次独立 blind judge 已按维度取平均；此表展示原始分，不应用 hard-gate 49 分封顶。",
            f"三臂原始总分：C0 {fmt(c0['RawTotal'])}，T0 {fmt(t0['RawTotal'])}，T1 {fmt(t1['RawTotal'])}。",
        ]
        for td, value in zip(targets[:14], values):
            ps = td.findall("p")
            assert ps, task
            p = ps[0]
            for child in list(p):
                p.remove(child)
            p.text = value
        updated += 1
    assert updated == 25, updated
    table_xml = ET.tostring(table, encoding="unicode", short_empty_elements=True)
    result = json.loads(run(LARK + [
        "docs", "+update", "--doc", DOC, "--command", "block_replace",
        "--block-id", TABLE_ID, "--content", "-"
    ], stdin=table_xml))
    assert result.get("ok") is True and result["data"].get("result") == "success", result
    print(json.dumps({"updated_rows": updated, "revision_id": result["data"]["document"]["revision_id"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
