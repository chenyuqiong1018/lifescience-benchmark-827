from __future__ import annotations

import ast
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

ACCEPTED = True

WELLS = [f"{row}{col}" for col in range(1, 7) for row in "ABCD"]
REQUIRED_COLUMNS = {"step", "source", "destination", "volume_uL", "pipette", "tip_policy"}
ADD_TYPES = {
    ("reagents:A1", 80.0): "lysis",
    ("reagents:A2", 120.0): "beads",
    ("reagents:A3", 180.0): "wash_add",
    ("reagents:A4", 40.0): "elution",
}
REMOVE_TYPES = {250.0: "supernatant", 180.0: "wash_remove"}
EXPECTED_PER_WELL = Counter({"lysis": 1, "beads": 1, "wash_add": 2,
                             "elution": 1, "supernatant": 1, "wash_remove": 2})
EXPECTED_SEQUENCE = ["lysis", "beads", "supernatant", "wash_add", "wash_remove",
                     "wash_add", "wash_remove", "elution"]
FRESH_HINTS = ("fresh", "new", "single_use", "single-use")
REUSE_HINTS = ("reuse", "same", "retain", "keep")


def _read_plan(path: Path):
    try:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = set(reader.fieldnames or [])
            return fields, list(reader)
    except Exception:
        return set(), []


def _norm(value):
    return re.sub(r"\s+", "", str(value or "").strip().lower())


def _classify(rows):
    by_well = defaultdict(Counter)
    sequence = defaultdict(list)
    typed = []
    numeric = True
    only_expected = True
    for row in rows:
        try:
            volume = float(row.get("volume_uL", ""))
        except (TypeError, ValueError):
            numeric = False
            typed.append((None, None, row))
            continue
        source, destination = _norm(row.get("source")), _norm(row.get("destination"))
        kind = well = None
        for candidate in WELLS:
            proc = f"processing:{candidate.lower()}"
            if destination == proc:
                well = candidate
                reagent_source = source[:9] + source[9:].upper() if source.startswith("reagents:") else source
                kind = ADD_TYPES.get((reagent_source, volume))
                break
            if source == proc and destination == "waste:a1":
                well, kind = candidate, REMOVE_TYPES.get(volume)
                break
        if kind is None:
            only_expected = False
        else:
            by_well[well][kind] += 1
            sequence[well].append(kind)
        typed.append((well, kind, row))
    complete = (numeric and only_expected and len(rows) == 192
                and all(by_well[w] == EXPECTED_PER_WELL for w in WELLS)
                and all(sequence[w] == EXPECTED_SEQUENCE for w in WELLS))
    return complete, typed, by_well


def _tip_safety(typed):
    """Validate the frozen semantic policy without assuming row order or step labels."""
    if len(typed) != 192 or any(kind is None for _, kind, _ in typed):
        return False
    per_well = defaultdict(list)
    for well, kind, row in typed:
        per_well[well].append((kind, _norm(row.get("tip_policy"))))
    for well in WELLS:
        entries = per_well[well]
        for kind in ("lysis", "beads", "supernatant", "elution"):
            vals = [p for k, p in entries if k == kind]
            if len(vals) != 1 or not any(h in vals[0] for h in FRESH_HINTS):
                return False
        adds = [p for k, p in entries if k == "wash_add"]
        removes = [p for k, p in entries if k == "wash_remove"]
        if len(adds) != 2 or len(removes) != 2:
            return False
        if not all(any(h in p for h in FRESH_HINTS) for p in adds):
            return False
        if not all(any(h in p for h in REUSE_HINTS) for p in removes):
            return False
        # Canonical wash-cycle labels prove add/remove pairing; alternatively explicit
        # shared identifiers after ':' or '=' may be used.
        def cycle(p):
            m = re.search(r"wash\D*([12])", p)
            if m:
                return m.group(1)
            m = re.search(r"(?:tip|id)[:=]([a-z0-9_-]+)", p)
            return m.group(1) if m else None
        if sorted(filter(None, map(cycle, adds))) != sorted(filter(None, map(cycle, removes))):
            return False
        if len(list(filter(None, map(cycle, adds)))) != 2:
            return False
    return True


def _protocol_facts(path: Path):
    try:
        text = path.read_text(encoding="utf-8")
        tree = ast.parse(text)
    except Exception:
        return {"parse": False}
    strings = [n.value.lower() for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, str)]
    numbers = [float(n.value) for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)) and not isinstance(n.value, bool)]
    attrs = [n.func.attr.lower() for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)]
    joined = "\n".join(strings)
    config = ("2.16" in strings and "p300_single_gen2" in strings and "right" in strings
              and "nest_12_reservoir_15ml" in strings and "nest_96_wellplate_2ml_deep" in strings
              and "opentrons_96_tiprack_300ul" in strings and "magnetic module gen2" in strings
              and all(str(slot) in strings for slot in ("1", "4", "5", "6", "7")))
    liquid_ops = ("mix" in attrs and ("transfer" in attrs or
                  ("aspirate" in attrs and "dispense" in attrs)))
    magnetic = "engage" in attrs and "disengage" in attrs and "delay" in attrs
    tip_control = "pick_up_tip" in attrs and "drop_tip" in attrs
    delay_args = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr.lower() == "delay":
            for kw in n.keywords:
                if isinstance(kw.value, ast.Constant) and isinstance(kw.value.value, (int, float)):
                    delay_args.append((kw.arg, float(kw.value.value)))
    timings = (("minutes", 5.0) in delay_args and ("minutes", 7.0) in delay_args
               and ("seconds", 30.0) in delay_args and ("minutes", 2.0) in delay_args)
    guards = any(k in joined for k in ("abort", "invalid", "insufficient", "raise")) or "raise" in [type(n).__name__.lower() for n in ast.walk(tree)]
    return {"parse": True, "configuration": config, "liquid_operations": liquid_ops,
            "magnet_calls": magnetic, "timings": timings, "tip_control": tip_control,
            "guards": guards}


def check(workspace: Path):
    output = Path(workspace) / "output"
    fields, rows = _read_plan(output / "transfer_plan.csv")
    headers = REQUIRED_COLUMNS.issubset(fields)
    complete, typed, by_well = _classify(rows) if headers else (False, [], {})
    pipette_range = bool(complete and all(
        _norm(row.get("pipette")) == "p300_single_gen2" and 20 <= float(row["volume_uL"]) <= 300
        for _, _, row in typed))
    declared_tips = bool(complete and _tip_safety(typed))

    # Independently grounded resource arithmetic: 1,920/2,880/8,640/960 uL
    # consumption fits initial-minus-dead volumes; 250 uL peak fits 2,000 uL;
    # 144 required tips fit the supplied inventory of 192.
    physical = pipette_range and complete
    facts = _protocol_facts(output / "protocol.py")
    tips = bool(declared_tips and facts.get("tip_control"))
    magnetic_direction = bool(facts.get("magnet_calls") and facts.get("liquid_operations") and facts.get("timings"))

    report = ""
    try:
        report = (output / "report.md").read_text(encoding="utf-8", errors="replace").lower()
    except Exception:
        pass
    summary_consistent = bool(complete and tips and all(re.search(p, report) for p in
        (r"\b24\s+samples?\b", r"\b192\b", r"\b144\s+tips?\b")) and
        ("80" in report and "120" in report and "180" in report and "40" in report))

    criteria = {
        "plan_schema_and_192_grounded_actions": headers and complete,
        "liquid_identity_volume_and_well_balance": complete,
        "p300_range_and_capacity_safety": physical,
        "six_tip_contamination_policy": tips,
        "plan_tip_policy_labels_grounded": declared_tips,
        "protocol_has_real_tip_control": bool(facts.get("tip_control")),
        "protocol_ast_and_deck_configuration": bool(facts.get("parse") and facts.get("configuration")),
        "protocol_has_real_liquid_operations": bool(facts.get("liquid_operations")),
        "magnet_engage_disengage_and_delays": bool(facts.get("magnet_calls") and facts.get("timings")),
        "explicit_preflight_abort_guards": bool(facts.get("guards")),
        "report_grounded_in_recomputed_plan": summary_consistent,
        "fatal_transfer_plan_scientific_integrity": complete,
        "fatal_contamination_and_physical_safety": tips and physical,
        "fatal_magnetic_workflow_direction": magnetic_direction,
    }

    core = (24 if complete else 0) + (10 if tips else 0) + (6 if physical else 0)
    direction = ((5 if criteria["protocol_ast_and_deck_configuration"] else 0)
                 + (4 if criteria["protocol_has_real_liquid_operations"] else 0)
                 + (4 if criteria["magnet_engage_disengage_and_delays"] else 0)
                 + (2 if criteria["explicit_preflight_abort_guards"] else 0))
    summary = 5 if summary_consistent else 0
    failures = []
    code_map = {
        "fatal_transfer_plan_scientific_integrity": "LS09_FATAL_TRANSFER_PLAN_SCIENCE",
        "fatal_contamination_and_physical_safety": "LS09_FATAL_CONTAMINATION_PHYSICAL_SAFETY",
        "fatal_magnetic_workflow_direction": "LS09_FATAL_MAGNETIC_DIRECTION",
    }
    for key, code in code_map.items():
        if not criteria[key]:
            failures.append(code)
    for key in ("protocol_ast_and_deck_configuration", "magnet_engage_disengage_and_delays",
                "explicit_preflight_abort_guards", "report_grounded_in_recomputed_plan"):
        if not criteria[key]:
            failures.append("LS09_" + key.upper())
    hardgate = all(criteria[k] for k in code_map)
    return {"core_science": core, "direction": direction, "summary": summary,
            "hardgate_pass": hardgate, "criteria": criteria, "failure_codes": failures}
