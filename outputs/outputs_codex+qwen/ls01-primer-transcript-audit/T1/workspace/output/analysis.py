#!/usr/bin/env python3
"""Audit primer pairs against the supplied transcript isoforms.

Inputs (the ONLY data used by this audit):
    inputs/primer_candidates.csv
        columns: pair_id, forward, reverse, expected_transcript, expected_product_bp
    inputs/transcripts.fa
        FASTA; headers carry metadata such as CDS=<start>-<end>.

Outputs:
    output/primer_audit.csv
        columns: pair_id, transcripts_matched, amplicon_length, cds_compatible,
                 status, reason
    output/report.md
        human-readable audit report.

Rules implemented here:
    * Only the supplied sequences are used; nothing is fetched from, or inferred
      against, any external database.
    * Primer matching is EXACT (no mismatches, no gaps). The reverse primer is
      matched as its reverse complement against the sense strand of each
      transcript.
    * A primer pair "matches" a transcript only when both primers bind in the
      correct orientation and order (forward upstream of the reverse-complement
      site), producing an amplicon.
    * Malformed or internally inconsistent sequence metadata (for example a CDS
      range that extends past the end of the sequence) is REPORTED, never
      silently repaired, clamped, or ignored.
"""

from __future__ import annotations

import csv
import re
from datetime import datetime
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parent.parent
INPUT_DIR = WORKSPACE / "inputs"
OUTPUT_DIR = WORKSPACE / "output"

PRIMER_CSV = INPUT_DIR / "primer_candidates.csv"
FASTA_PATH = INPUT_DIR / "transcripts.fa"
AUDIT_CSV = OUTPUT_DIR / "primer_audit.csv"
REPORT_MD = OUTPUT_DIR / "report.md"

COMPLEMENT = str.maketrans("ACGT", "TGCA")
CDS_RE = re.compile(r"CDS=(\d+)-(\d+)")
CSV_COLUMNS = ["pair_id", "transcripts_matched", "amplicon_length",
               "cds_compatible", "status", "reason"]


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


# ---------------------------------------------------------------------
# FASTA parsing and metadata validation
# ---------------------------------------------------------------------

def parse_fasta(path: Path):
    records = []
    current = None
    chunks = []

    def close():
        nonlocal current, chunks
        if current is not None:
            current["seq"] = "".join(chunks).upper()
            records.append(current)
        current, chunks = None, []

    with path.open(encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                close()
                header = line[1:].strip()
                ident = header.split()[0] if header else "NOID_LINE%d" % lineno
                current = {
                    "id": ident,
                    "header": header,
                    "seq": "",
                    "cds": None,
                    "seq_issues": [],
                    "cds_issues": [],
                }
                m = CDS_RE.search(header)
                if m:
                    current["cds"] = (int(m.group(1)), int(m.group(2)))
            else:
                if current is None:
                    raise ValueError(
                        "%s line %d: sequence data before any header" % (path.name, lineno))
                seq = line.upper()
                bad = sorted(set(seq) - set("ACGT"))
                if bad:
                    current["seq_issues"].append(
                        "invalid sequence characters %s at line %d" % (bad, lineno))
                chunks.append(seq)
    close()
    return records


def validate_record(rec):
    """Flag malformed / internally inconsistent metadata. Never repair it."""
    n = len(rec["seq"])
    if n == 0:
        rec["seq_issues"].append("sequence is empty")
    cds = rec["cds"]
    if cds is None:
        rec["cds_issues"].append("no CDS=<start>-<end> annotation present in header")
    else:
        s, e = cds
        if s < 1:
            rec["cds_issues"].append("CDS start %d is < 1" % s)
        if s > e:
            rec["cds_issues"].append("CDS start %d > CDS end %d" % (s, e))
        elif e > n:
            rec["cds_issues"].append(
                "annotated CDS=%d-%d extends past the end of the %d bp sequence "
                "(annotated CDS span of %d bp cannot fit inside a %d bp transcript)"
                % (s, e, n, e - s + 1, n))
        elif (e - s + 1) % 3 != 0:
            rec["cds_issues"].append(
                "CDS span %d bp is not a multiple of 3" % (e - s + 1))
    rec["consistent"] = not rec["seq_issues"] and not rec["cds_issues"]
    return rec


# ---------------------------------------------------------------------
# Primer matching
# ---------------------------------------------------------------------

def find_all(haystack: str, needle: str):
    """All 1-based start positions of exact occurrences of needle in haystack."""
    hits = []
    i = haystack.find(needle)
    while i != -1:
        hits.append(i + 1)
        i = haystack.find(needle, i + 1)
    return hits


def audit_pair(pair, records):
    fwd = pair["forward"]
    rev_rc = revcomp(pair["reverse"])
    audit = {"fwd_hits": {}, "rev_hits": {}, "amplicons": {}}
    for rec in records:
        f_hits = find_all(rec["seq"], fwd)
        r_hits = find_all(rec["seq"], rev_rc)
        audit["fwd_hits"][rec["id"]] = f_hits
        audit["rev_hits"][rec["id"]] = r_hits
        amps = []
        for f in f_hits:
            for r in r_hits:
                r_end = r + len(pair["reverse"]) - 1
                if r > f:  # reverse site downstream of forward primer
                    amps.append({"start": f, "end": r_end,
                                 "length": r_end - f + 1})
        if amps:
            audit["amplicons"][rec["id"]] = amps
    return audit


def cds_compatible_verdict(rec, amps):
    """Return (value, explanation) for the cds_compatible column."""
    if not amps:
        return "na", "no amplicon, nothing to evaluate against the CDS"
    if rec["cds"] is None:
        return "na", "transcript header carries no CDS annotation"
    if rec["cds_issues"]:
        s, e = rec["cds"]
        a = amps[0]
        if a["start"] >= s and a["end"] <= e:
            face = ("even taking the annotation at face value the amplicon would sit "
                    "inside the CDS, but the annotation itself contradicts the sequence length")
        else:
            face = ("even taking the annotation at face value, amplicon %d..%d is not "
                    "contained within CDS %d..%d" % (a["start"], a["end"], s, e))
        expl = ("CDS annotation CDS=%d-%d is internally inconsistent with the %d bp "
                "sequence (%s); compatibility cannot be assessed; metadata reported, "
                "not repaired" % (s, e, len(rec["seq"]), "; ".join(rec["cds_issues"])))
        return "inconsistent_metadata", expl + "; " + face
    s, e = rec["cds"]
    if all(a["start"] >= s and a["end"] <= e for a in amps):
        return "yes", "amplicon lies within annotated CDS %d..%d" % (s, e)
    return "no", "amplicon extends outside annotated CDS %d..%d" % (s, e)


def evaluate_pair(pair, audit, records):
    rec_by_id = {r["id"]: r for r in records}
    expected = pair["expected_transcript"]
    expected_bp = pair["expected_product_bp"]
    matched = list(audit["amplicons"].keys())
    exp_rec = rec_by_id.get(expected)
    exp_amps = audit["amplicons"].get(expected, []) if exp_rec else []

    sentences = []
    amp_len = ""

    # --- primer binding summary on the expected transcript -------------
    if exp_rec is not None:
        lf, lr = len(pair["forward"]), len(pair["reverse"])
        fhits = audit["fwd_hits"].get(expected, [])
        rhits = audit["rev_hits"].get(expected, [])
        if fhits:
            coord = ", ".join("%d..%d" % (p, p + lf - 1) for p in fhits)
            sentences.append("forward primer binds %s %sat %s"
                             % (expected, "uniquely " if len(fhits) == 1 else "", coord))
        else:
            sentences.append("forward primer has no exact match in %s" % expected)
        if rhits:
            coord = ", ".join("%d..%d" % (p, p + lr - 1) for p in rhits)
            sentences.append("reverse primer (reverse-complemented) binds %s %sat %s"
                             % (expected, "uniquely " if len(rhits) == 1 else "", coord))
        else:
            sentences.append("reverse primer (reverse-complemented) has no exact match in %s"
                             % expected)
    else:
        sentences.append("expected transcript %s is not present among the supplied sequences"
                         % expected)

    # --- status ----------------------------------------------------------
    if pair.get("input_issue"):
        status = "malformed_input"
        sentences.append(pair["input_issue"])
    elif exp_rec is None:
        status = "expected_transcript_missing"
    elif not matched:
        status = "no_binding"
        if not any(audit["fwd_hits"].values()):
            sentences.append("forward primer matches none of the %d supplied transcripts"
                             % len(records))
        if not any(audit["rev_hits"].values()):
            sentences.append("reverse-complemented reverse primer matches none of the %d "
                             "supplied transcripts" % len(records))
        sentences.append("no amplicon can be formed on any supplied isoform")
    elif not exp_amps:
        status = "off_target"
        amp_len = audit["amplicons"][matched[0]][0]["length"]
        sentences.append("amplicon formed only on non-expected transcript(s): %s"
                         % "; ".join(matched))
    elif (len(exp_amps) > 1 or len(audit["fwd_hits"][expected]) > 1
          or len(audit["rev_hits"][expected]) > 1):
        status = "ambiguous"
        amp_len = exp_amps[0]["length"]
        sentences.append("%d candidate amplicons on %s; primer binding is not unique"
                         % (len(exp_amps), expected))
    else:
        amp = exp_amps[0]
        amp_len = amp["length"]
        sentences.append("predicted amplicon %s:%d..%d, length %d bp"
                         % (expected, amp["start"], amp["end"], amp["length"]))
        if expected_bp is not None and amp["length"] != expected_bp:
            status = "length_mismatch"
            sentences.append("amplicon length %d bp differs from expected_product_bp %d "
                             "(delta %+d bp)"
                             % (amp["length"], expected_bp, amp["length"] - expected_bp))
        else:
            status = None  # resolved below via CDS verdict

    # --- CDS compatibility ------------------------------------------------
    if exp_rec is not None and exp_amps:
        cds_value, cds_expl = cds_compatible_verdict(exp_rec, exp_amps)
        sentences.append("cds_compatible=%s: %s" % (cds_value, cds_expl))
    else:
        cds_value, cds_expl = "na", "no amplicon on the expected transcript"
        if status not in (None,):
            pass

    if status is None:
        # length matched expectation; decide via CDS verdict
        if cds_value == "yes":
            status = "pass"
        elif cds_value == "inconsistent_metadata":
            status = "metadata_conflict"
        else:
            status = "cds_incompatible"

    # --- specificity against the other supplied transcripts ----------------
    others = [r["id"] for r in records if r["id"] != expected]
    if exp_rec is not None and others:
        off = [t for t in others if t in audit["amplicons"]]
        if off:
            sentences.append("WARNING: also amplifies %s (not isoform-specific)"
                             % "; ".join(off))
        else:
            sentences.append("no amplicon on other supplied transcript(s): %s"
                             % "; ".join(others))

    # --- metadata issues of the expected transcript -------------------------
    if exp_rec is not None and (exp_rec["seq_issues"] or exp_rec["cds_issues"]):
        for issue in exp_rec["seq_issues"] + exp_rec["cds_issues"]:
            sentences.append("%s header metadata issue (reported, not repaired): %s"
                             % (expected, issue))

    reason = ". ".join(s.strip().rstrip(".") for s in sentences) + "."
    return {
        "pair": pair,
        "audit": audit,
        "csv_row": {
            "pair_id": pair["pair_id"],
            "transcripts_matched": ";".join(matched) if matched else "none",
            "amplicon_length": amp_len if amp_len != "" else "na",
            "cds_compatible": cds_value,
            "status": status,
            "reason": reason,
        },
        "expected_record": exp_rec,
        "exp_amps": exp_amps,
        "matched": matched,
    }


# ---------------------------------------------------------------------
# Input loading
# ---------------------------------------------------------------------

def load_primers(path: Path):
    pairs = []
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        required = ["pair_id", "forward", "reverse",
                    "expected_transcript", "expected_product_bp"]
        missing = [c for c in required if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError("primer CSV missing columns: %s" % missing)
        for row in reader:
            pair = {(k or "").strip(): (v or "").strip() for k, v in row.items()}
            issues = []
            pair["forward"] = pair["forward"].upper()
            pair["reverse"] = pair["reverse"].upper()
            if not re.fullmatch(r"[ACGT]+", pair["forward"]):
                issues.append("forward primer '%s' is not valid DNA (ACGT only)"
                              % pair["forward"])
            if not re.fullmatch(r"[ACGT]+", pair["reverse"]):
                issues.append("reverse primer '%s' is not valid DNA (ACGT only)"
                              % pair["reverse"])
            try:
                pair["expected_product_bp"] = int(pair["expected_product_bp"])
            except ValueError:
                issues.append("expected_product_bp '%s' is not an integer"
                              % pair["expected_product_bp"])
                pair["expected_product_bp"] = None
            pair["input_issue"] = "; ".join(issues) if issues else None
            pairs.append(pair)
    return pairs


# ---------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------

def build_report(records, pairs, evaluations):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M")
    L = []
    L.append("# Primer audit report")
    L.append("")
    L.append("Generated: %s by `output/analysis.py`" % ts)
    L.append("")
    L.append("Scope: %d primer pair(s) audited against %d supplied transcript isoform(s). "
             "Only the supplied files `inputs/primer_candidates.csv` and "
             "`inputs/transcripts.fa` were used; no external sequences or databases "
             "were consulted. Matching is exact (no mismatches/gaps); the reverse primer "
             "is matched as its reverse complement against the sense strand."
             % (len(pairs), len(records)))
    L.append("")

    # ---------------- metadata integrity ----------------
    L.append("## 1. Transcript metadata integrity check")
    L.append("")
    bad = [r for r in records if not r["consistent"]]
    if bad:
        L.append("**Internally inconsistent metadata was detected and is reported below. "
                 "Per the audit policy it was NOT repaired, clamped, or ignored.**")
    else:
        L.append("All supplied transcript metadata is internally consistent.")
    L.append("")
    L.append("| transcript | header | seq length (bp) | CDS annotation | metadata verdict |")
    L.append("|---|---|---|---|---|")
    for r in records:
        cds = ("%d..%d (span %d bp)" % (r["cds"][0], r["cds"][1],
                                        r["cds"][1] - r["cds"][0] + 1)) if r["cds"] else "none"
        verdict = "consistent" if r["consistent"] else "INCONSISTENT"
        L.append("| %s | `%s` | %d | %s | %s |"
                 % (r["id"], r["header"], len(r["seq"]), cds, verdict))
    L.append("")
    for r in records:
        if r["consistent"]:
            continue
        L.append("- **%s**: %s. The sequence was used exactly as supplied."
                 % (r["id"], "; ".join(r["seq_issues"] + r["cds_issues"])))
    if bad:
        L.append("")
        L.append("Consequence: for any amplicon mapping to a transcript with an inconsistent "
                 "CDS annotation, `cds_compatible` is reported as `inconsistent_metadata` "
                 "instead of yes/no, because the annotation cannot support a reliable "
                 "compatibility call.")
    L.append("")

    # ---------------- methods ----------------
    L.append("## 2. Method")
    L.append("")
    L.append("- Forward primer: exact substring search on the sense strand of each transcript.")
    L.append("- Reverse primer: reverse-complemented, then exact substring search on the sense strand.")
    L.append("- A pair matches a transcript only if a forward site lies upstream of a "
             "reverse-complement site; amplicon = forward_start .. rev_site_end, "
             "length = end - start + 1.")
    L.append("- `cds_compatible` = `yes` only if the amplicon lies fully inside a "
             "self-consistent annotated CDS; `inconsistent_metadata` if the CDS annotation "
             "contradicts the sequence; `na` if there is no amplicon / no CDS annotation.")
    L.append("- Status vocabulary: `pass`, `length_mismatch`, `no_binding`, `off_target`, "
             "`ambiguous`, `cds_incompatible`, `metadata_conflict`, "
             "`expected_transcript_missing`, `malformed_input`.")
    L.append("")

    # ---------------- per-pair results ----------------
    L.append("## 3. Per-pair audit results")
    for ev in evaluations:
        pair = ev["pair"]
        audit = ev["audit"]
        row = ev["csv_row"]
        L.append("")
        L.append("### %s (expected transcript: %s, expected product: %s bp)"
                 % (pair["pair_id"], pair["expected_transcript"],
                    pair["expected_product_bp"] if pair["expected_product_bp"] is not None else "?"))
        L.append("")
        L.append("- forward: `%s` (%d nt)" % (pair["forward"], len(pair["forward"])))
        L.append("- reverse: `%s` (%d nt), RC = `%s`" % (pair["reverse"], len(pair["reverse"]),
                                                          revcomp(pair["reverse"])))
        L.append("")
        L.append("| transcript | forward hits | RC(reverse) hits | amplicon(s) |")
        L.append("|---|---|---|---|")
        for rec in records:
            tid = rec["id"]
            fh = ", ".join(str(p) for p in audit["fwd_hits"][tid]) or "-"
            rh = ", ".join(str(p) for p in audit["rev_hits"][tid]) or "-"
            amps = audit["amplicons"].get(tid)
            amptxt = ", ".join("%d..%d (%d bp)" % (a["start"], a["end"], a["length"])
                               for a in amps) if amps else "-"
            L.append("| %s | %s | %s | %s |" % (tid, fh, rh, amptxt))
        L.append("")
        L.append("- transcripts_matched: **%s**" % row["transcripts_matched"])
        L.append("- amplicon_length: **%s**" % row["amplicon_length"])
        L.append("- cds_compatible: **%s**" % row["cds_compatible"])
        L.append("- status: **%s**" % row["status"])
        L.append("- reason: %s" % row["reason"])
    L.append("")

    # ---------------- summary table ----------------
    L.append("## 4. Summary (contents of `output/primer_audit.csv`)")
    L.append("")
    L.append("| " + " | ".join(CSV_COLUMNS) + " |")
    L.append("|" + "---|" * len(CSV_COLUMNS))
    for ev in evaluations:
        r = ev["csv_row"]
        L.append("| " + " | ".join(str(r[c]) for c in CSV_COLUMNS) + " |")
    L.append("")

    # ---------------- conclusions ----------------
    L.append("## 5. Conclusions")
    L.append("")
    n_pass = sum(1 for ev in evaluations if ev["csv_row"]["status"] == "pass")
    L.append("- %d of %d primer pairs fully pass the audit." % (n_pass, len(evaluations)))
    for ev in evaluations:
        pair, row = ev["pair"], ev["csv_row"]
        if row["status"] == "length_mismatch":
            L.append("- **%s**: primers bind the expected transcript %s uniquely, but the "
                     "predicted product (%s bp) does not match expected_product_bp (%s); "
                     "the pair should not be used as-is."
                     % (pair["pair_id"], pair["expected_transcript"],
                        row["amplicon_length"], pair["expected_product_bp"]))
        elif row["status"] == "no_binding":
            L.append("- **%s**: neither primer has an exact binding site in any supplied "
                     "transcript; the pair cannot amplify the supplied isoforms and should "
                     "be discarded." % pair["pair_id"])
        elif row["status"] != "pass":
            L.append("- **%s**: status %s (see reason field)." % (pair["pair_id"], row["status"]))
    if bad:
        ids = ", ".join(r["id"] for r in bad)
        L.append("- Metadata: transcript(s) %s carry CDS annotations that exceed the actual "
                 "sequence length. Either the FASTA sequences are truncated excerpts or the "
                 "CDS coordinates are wrong; this cannot be resolved from the supplied data "
                 "alone, so the inconsistency is reported rather than repaired." % ids)
    L.append("")
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    records = [validate_record(r) for r in parse_fasta(FASTA_PATH)]
    if not records:
        raise SystemExit("no records parsed from %s" % FASTA_PATH)
    pairs = load_primers(PRIMER_CSV)
    if not pairs:
        raise SystemExit("no primer pairs parsed from %s" % PRIMER_CSV)

    evaluations = []
    for pair in pairs:
        if pair["input_issue"]:
            audit = {"fwd_hits": {r["id"]: [] for r in records},
                     "rev_hits": {r["id"]: [] for r in records},
                     "amplicons": {}}
        else:
            audit = audit_pair(pair, records)
        evaluations.append(evaluate_pair(pair, audit, records))

    with AUDIT_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        w.writeheader()
        for ev in evaluations:
            w.writerow(ev["csv_row"])

    REPORT_MD.write_text(build_report(records, pairs, evaluations), encoding="utf-8")

    print("Transcripts:")
    for r in records:
        print("  %-14s len=%-4d consistent=%s issues=%s"
              % (r["id"], len(r["seq"]), r["consistent"],
                 r["seq_issues"] + r["cds_issues"]))
    print("Audit rows written to %s" % AUDIT_CSV)
    for ev in evaluations:
        row = ev["csv_row"]
        print("  %s matched=%-14s amp=%-4s cds=%-22s status=%s"
              % (row["pair_id"], row["transcripts_matched"],
                 row["amplicon_length"], row["cds_compatible"], row["status"]))
    print("Report written to %s" % REPORT_MD)


if __name__ == "__main__":
    main()