#!/usr/bin/env python3
"""Construct audit analysis.

Audits each construct in ``inputs/constructs.csv`` for start/stop codons,
reading frame, and tag compatibility, following the frozen rules in
``inputs/AUDIT_RULE.md``. Only fields represented in the input fixture are
checked; no vector sequence or other features are inferred.

Frozen rules (inputs/AUDIT_RULE.md):
  * start_ok : the insert begins with ``ATG``.
  * stop_ok  : the insert ends in an in-frame stop codon (TAA/TAG/TGA).
  * frame_ok : insert length is divisible by three AND ``claimed_frame`` is
               ``in_frame``.
  * tag_ok   : for ``C_terminal_*`` fusions the insert must not contain a
               terminal stop codon before the downstream tag; other tag
               strings are unsupported and fail closed.
  * overall_status = ``pass`` only when all four checks are true; otherwise
    ``fail`` and every failed check is listed in ``issues`` using the labels
    START / STOP / FRAME / TAG.
"""

import csv
from pathlib import Path

STOP_CODONS = {"TAA", "TAG", "TGA"}

# Ordered (check_key, issue_label) pairs; issues are reported in this order.
CHECK_LABELS = [
    ("start_ok", "START"),
    ("stop_ok", "STOP"),
    ("frame_ok", "FRAME"),
    ("tag_ok", "TAG"),
]

CSV_FIELDS = [
    "construct_id",
    "frame_ok",
    "start_ok",
    "stop_ok",
    "tag_ok",
    "overall_status",
    "issues",
]


def audit_construct(row):
    """Apply the frozen audit rules to a single construct row."""
    construct_id = (row.get("construct_id") or "").strip()
    insert = (row.get("insert_sequence") or "").strip().upper()
    tag = (row.get("tag") or "").strip()
    claimed_frame = (row.get("claimed_frame") or "").strip()
    promoter = (row.get("promoter") or "").strip()

    n = len(insert)

    # start_ok: the insert begins with ATG.
    start_ok = insert.startswith("ATG")

    # stop_ok: the insert ends in an in-frame stop codon. The reading frame
    # starts at the first nucleotide of the supplied insert (start codon), so
    # an in-frame terminal codon requires the length to be divisible by 3.
    stop_ok = (n % 3 == 0) and (n >= 3) and (insert[-3:] in STOP_CODONS)

    # frame_ok: insert length divisible by three AND claimed_frame is in_frame.
    frame_ok = (n % 3 == 0) and (claimed_frame == "in_frame")

    # tag_ok: for C_terminal_* fusions the insert must not contain a terminal
    # stop codon before the downstream tag; other tag strings fail closed.
    if tag.startswith("C_terminal_"):
        terminal_stop = (n >= 3) and (insert[-3:] in STOP_CODONS)
        tag_ok = not terminal_stop
    else:
        tag_ok = False

    checks = {
        "start_ok": start_ok,
        "stop_ok": stop_ok,
        "frame_ok": frame_ok,
        "tag_ok": tag_ok,
    }
    issues = [label for key, label in CHECK_LABELS if not checks[key]]
    overall_status = "pass" if not issues else "fail"

    return {
        "construct_id": construct_id,
        "frame_ok": frame_ok,
        "start_ok": start_ok,
        "stop_ok": stop_ok,
        "tag_ok": tag_ok,
        "overall_status": overall_status,
        "issues": ";".join(issues),
        # Extra metadata used only for the report (not written to the CSV).
        "_insert": insert,
        "_length": n,
        "_tag": tag,
        "_claimed_frame": claimed_frame,
        "_promoter": promoter,
    }


def explain(r):
    """Human-readable explanation for each failed check of one construct."""
    out = []
    insert = r["_insert"]
    n = r["_length"]
    if not r["start_ok"]:
        out.append("START: insert does not begin with `ATG`.")
    if not r["stop_ok"]:
        if n % 3 != 0:
            out.append(
                "STOP: insert length %d is not divisible by 3, so the insert "
                "does not end in an in-frame stop codon." % n
            )
        else:
            out.append(
                "STOP: terminal codon `%s` is not a stop codon (TAA/TAG/TGA)."
                % insert[-3:]
            )
    if not r["frame_ok"]:
        reasons = []
        if n % 3 != 0:
            reasons.append("insert length %d is not divisible by 3" % n)
        if r["_claimed_frame"] != "in_frame":
            reasons.append("claimed_frame is `%s`, not `in_frame`" % r["_claimed_frame"])
        out.append("FRAME: " + " and ".join(reasons) + ".")
    if not r["tag_ok"]:
        if r["_tag"].startswith("C_terminal_"):
            out.append(
                "TAG: `%s` fusion but the insert ends in stop codon `%s`, "
                "which would terminate translation before the downstream tag."
                % (r["_tag"], insert[-3:] if n >= 3 else "n/a")
            )
        else:
            out.append("TAG: tag `%s` is unsupported and fails closed." % r["_tag"])
    if not out:
        out.append("All four checks passed.")
    return out


def write_csv(results, output_csv):
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for r in results:
            writer.writerow({k: r[k] for k in CSV_FIELDS})


def write_report(results, report_md):
    total = len(results)
    n_pass = sum(1 for r in results if r["overall_status"] == "pass")
    n_fail = total - n_pass

    L = []
    L.append("# Construct Audit Report")
    L.append("")
    L.append(
        "Audit of the constructs in `inputs/constructs.csv` against the frozen "
        "rules in `inputs/AUDIT_RULE.md`. Only fields represented in the input "
        "fixture are checked; no vector sequence or other sequence features are "
        "inferred."
    )
    L.append("")
    L.append("## Rules applied")
    L.append("")
    L.append("| Check | Rule |")
    L.append("| --- | --- |")
    L.append("| `start_ok` | Insert begins with `ATG`. |")
    L.append("| `stop_ok` | Insert ends in an in-frame stop codon (`TAA`/`TAG`/`TGA`). |")
    L.append("| `frame_ok` | Insert length divisible by 3 **and** `claimed_frame` is `in_frame`. |")
    L.append(
        "| `tag_ok` | For `C_terminal_*` fusions the insert must not contain a "
        "terminal stop codon before the downstream tag; other tag strings are "
        "unsupported and fail closed. |"
    )
    L.append(
        "| `overall_status` | `pass` only when all four checks are true; otherwise "
        "`fail`, listing every failed check in `issues` (labels START/STOP/FRAME/TAG). |"
    )
    L.append("")
    L.append("## Summary")
    L.append("")
    L.append("- Constructs audited: %d" % total)
    L.append("- Pass: %d" % n_pass)
    L.append("- Fail: %d" % n_fail)
    L.append("")
    L.append("## Results")
    L.append("")
    L.append("| construct_id | frame_ok | start_ok | stop_ok | tag_ok | overall_status | issues |")
    L.append("| --- | --- | --- | --- | --- | --- | --- |")
    for r in results:
        L.append(
            "| %s | %s | %s | %s | %s | %s | %s |"
            % (
                r["construct_id"],
                r["frame_ok"],
                r["start_ok"],
                r["stop_ok"],
                r["tag_ok"],
                r["overall_status"],
                r["issues"] if r["issues"] else "-",
            )
        )
    L.append("")
    L.append("## Per-construct detail")
    for r in results:
        L.append("")
        L.append("### %s" % r["construct_id"])
        L.append("")
        L.append("- promoter: `%s`" % r["_promoter"])
        L.append("- tag: `%s`" % r["_tag"])
        L.append("- claimed_frame: `%s`" % r["_claimed_frame"])
        L.append("- insert_sequence: `%s`" % r["_insert"])
        L.append("- insert length: %d nt (length mod 3 = %d)" % (r["_length"], r["_length"] % 3))
        L.append("- start_ok: %s" % r["start_ok"])
        L.append("- stop_ok: %s" % r["stop_ok"])
        L.append("- frame_ok: %s" % r["frame_ok"])
        L.append("- tag_ok: %s" % r["tag_ok"])
        L.append("- overall_status: **%s**" % r["overall_status"])
        L.append("- issues: %s" % (r["issues"] if r["issues"] else "(none)"))
        L.append("")
        L.append("Explanation:")
        L.append("")
        for e in explain(r):
            L.append("- %s" % e)
    L.append("")
    L.append("## Scope notes")
    L.append("")
    L.append(
        "- The fixture supplies only `construct_id`, `promoter`, `insert_sequence`, "
        "`tag`, and `claimed_frame`. The frozen rules define checks over the insert "
        "sequence, tag, and claimed frame; `promoter` is recorded but has no audit "
        "rule attached."
    )
    L.append(
        "- No linker sequence and no cloning flags (e.g., restriction sites) are "
        "represented in the input, so none were evaluated or inferred, per the rule "
        "that only fields represented by the fixture are checked."
    )
    L.append(
        "- For `C_terminal_*` fusions the rules make `stop_ok` (insert ends in an "
        "in-frame stop) and `tag_ok` (no terminal stop before the downstream tag) "
        "mutually exclusive, so a C-terminal fusion construct cannot pass both; each "
        "check is reported independently as specified."
    )
    L.append("")

    with open(report_md, "w", encoding="utf-8") as f:
        f.write("\n".join(L))


def main():
    script_dir = Path(__file__).resolve().parent
    workspace_root = script_dir.parent
    input_csv = workspace_root / "inputs" / "constructs.csv"
    output_csv = script_dir / "construct_audit.csv"
    report_md = script_dir / "report.md"

    with open(input_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        results = [audit_construct(row) for row in reader]

    write_csv(results, output_csv)
    write_report(results, report_md)

    for r in results:
        print("%s: %s issues=[%s]" % (r["construct_id"], r["overall_status"], r["issues"]))
    print("Wrote: %s" % output_csv)
    print("Wrote: %s" % report_md)


if __name__ == "__main__":
    main()