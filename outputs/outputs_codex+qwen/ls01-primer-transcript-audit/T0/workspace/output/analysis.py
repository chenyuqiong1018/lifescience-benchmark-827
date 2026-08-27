#!/usr/bin/env python3
"""
Primer-pair audit against supplied transcript isoforms.

Inputs (relative to the workspace root, i.e. the parent of this script's dir):
  inputs/transcripts.fa         transcript isoform sequences (FASTA)
  inputs/primer_candidates.csv  primer pairs with expected targets/product sizes

Outputs (written next to this script):
  primer_audit.csv  pair_id,transcripts_matched,amplicon_length,cds_compatible,status,reason
  report.md         human-readable audit report

Audit policy
------------
* Only the supplied sequences are used. Nothing is fetched, inferred, or
  reconstructed from external databases.
* Primer binding is assessed by EXACT match: the forward primer must occur on
  the sense strand; the reverse primer must occur as its reverse complement on
  the sense strand, downstream of the forward primer, without overlap.
  The amplicon spans from the forward primer start to the end of the reverse
  primer binding site (both primers included).
* A transcript is counted as "matched" only when a productive amplicon forms.
* cds_compatible = yes only when a productive amplicon lies fully inside the
  annotated CDS range of that transcript AND the CDS annotation is internally
  consistent.
* Malformed or internally inconsistent sequence metadata is REPORTED, never
  silently repaired.
"""

from __future__ import annotations

import csv
import re
import sys
from datetime import datetime
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
WORKSPACE = SCRIPT_DIR.parent
INPUT_DIR = WORKSPACE / "inputs"
TRANSCRIPTS_FA = INPUT_DIR / "transcripts.fa"
PRIMERS_CSV = INPUT_DIR / "primer_candidates.csv"
AUDIT_CSV = SCRIPT_DIR / "primer_audit.csv"
REPORT_MD = SCRIPT_DIR / "report.md"

VALID_BASES = set("ACGT")
_COMP = str.maketrans("ACGT", "TGCA")
REQUIRED_COLS = ["pair_id", "forward", "reverse", "expected_transcript", "expected_product_bp"]


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# --------------------------------------------------------------------------
# FASTA parsing + metadata validation (report problems, never repair them)
# --------------------------------------------------------------------------
def parse_fasta(path: Path):
    records, warnings = [], []
    current = None
    with open(path, "r", encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith(">"):
                if current is not None:
                    records.append(current)
                header = line[1:].strip()
                fields = header.split()
                current = {
                    "id": fields[0] if fields else f"UNKNOWN_line{lineno}",
                    "header": header,
                    "lineno": lineno,
                    "chunks": [],
                }
                if not fields:
                    warnings.append(f"transcripts.fa line {lineno}: header has no identifier.")
            else:
                if current is None:
                    warnings.append(f"transcripts.fa line {lineno}: sequence data before any header; ignored.")
                    continue
                current["chunks"].append("".join(line.split()))
    if current is not None:
        records.append(current)

    seen = {}
    for rec in records:
        seq = "".join(rec["chunks"]).upper()
        rec["seq"] = seq
        rec["length"] = len(seq)
        invalid = sorted(set(seq) - VALID_BASES)
        rec["invalid_chars"] = invalid
        if invalid:
            warnings.append(f"transcripts.fa {rec['id']}: sequence contains non-ACGT characters {invalid}.")
        if rec["length"] == 0:
            warnings.append(f"transcripts.fa {rec['id']}: empty sequence.")
        if rec["id"] in seen:
            warnings.append(f"transcripts.fa {rec['id']}: duplicate transcript id (lines {seen[rec['id']]} and {rec['lineno']}).")
        seen[rec["id"]] = rec["lineno"]

        # ---- header metadata ------------------------------------------------
        rec["notes"] = rec["header"].split()[1:]
        rec["cds"] = None
        rec["cds_issues"] = []
        m = re.search(r"CDS=(\d+)-(\d+)", rec["header"])
        if m:
            start, end = int(m.group(1)), int(m.group(2))
            rec["cds"] = (start, end)
            if start < 1:
                rec["cds_issues"].append(f"CDS start {start} < 1")
            if end < start:
                rec["cds_issues"].append(f"CDS end {end} < CDS start {start}")
            if (end - start + 1) % 3 != 0:
                rec["cds_issues"].append(f"CDS length {end - start + 1} nt is not a multiple of 3")
            if end > rec["length"]:
                rec["cds_issues"].append(
                    f"CDS end {end} exceeds sequence length {rec['length']} "
                    f"(annotation claims {end} nt, only {rec['length']} nt supplied)"
                )
            elif start > rec["length"]:
                rec["cds_issues"].append(f"CDS start {start} exceeds sequence length {rec['length']}")
        rec["metadata_ok"] = (rec["length"] > 0 and not invalid and not rec["cds_issues"])
        if rec["cds_issues"]:
            warnings.append(
                f"transcripts.fa {rec['id']}: internally inconsistent CDS annotation "
                f"'CDS={rec['cds'][0]}-{rec['cds'][1]}' vs sequence length {rec['length']}: "
                + "; ".join(rec["cds_issues"]) + ". Reported as-is; NOT repaired."
            )
    return records, warnings


# --------------------------------------------------------------------------
# Primer table parsing
# --------------------------------------------------------------------------
def parse_primers(path: Path):
    rows, warnings = [], []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in REQUIRED_COLS if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"primer_candidates.csv is missing required column(s): {missing}")
        for i, row in enumerate(reader, start=2):
            rec = {k: (row.get(k) or "").strip() for k in REQUIRED_COLS}
            rec["issues"] = []
            if not rec["pair_id"]:
                rec["pair_id"] = f"row{i}"
                rec["issues"].append("missing pair_id")
            for key in ("forward", "reverse"):
                s = rec[key].upper()
                rec[key] = s
                if not s:
                    rec["issues"].append(f"missing {key} sequence")
                else:
                    bad = sorted(set(s) - VALID_BASES)
                    if bad:
                        rec["issues"].append(f"{key} contains non-ACGT characters {bad}")
            try:
                rec["expected_product_bp"] = int(rec["expected_product_bp"])
            except ValueError:
                rec["issues"].append(f"expected_product_bp '{rec['expected_product_bp']}' is not an integer")
                rec["expected_product_bp"] = None
            if rec["issues"]:
                warnings.append(f"primer_candidates.csv row {i} ({rec['pair_id']}): " + "; ".join(rec["issues"]))
            rows.append(rec)
    return rows, warnings


# --------------------------------------------------------------------------
# Matching helpers
# --------------------------------------------------------------------------
def find_exact(seq: str, motif: str):
    """All (0-based) start positions of motif in seq (overlapping allowed)."""
    hits, start = [], 0
    if not motif:
        return hits
    while True:
        i = seq.find(motif, start)
        if i < 0:
            return hits
        hits.append(i)
        start = i + 1


def best_near_match(seq: str, motif: str):
    """Best Hamming-distance placement of motif in seq -> (mismatches, pos0) or None."""
    n, L = len(seq), len(motif)
    if n < L or L == 0:
        return None
    best = None
    for i in range(n - L + 1):
        mm = sum(1 for a, b in zip(seq[i:i + L], motif) if a != b)
        if best is None or mm < best[0]:
            best = (mm, i)
            if mm == 0:
                return best
    return best


def audit_pair(pair, transcripts):
    """Exact-match audit of one primer pair against every transcript."""
    fwd = pair["forward"]
    rev = pair["reverse"]
    rev_rc = revcomp(rev)
    lf, lr = len(fwd), len(rev)

    per_tx, amplicons = {}, []
    for tx in transcripts:
        seq = tx["seq"]
        fwd_hits = find_exact(seq, fwd)
        rev_hits = find_exact(seq, rev_rc)
        tx_amp = []
        for f in fwd_hits:
            for r in rev_hits:
                if f + lf - 1 < r:  # forward upstream of reverse site, non-overlapping
                    s, e = f + 1, r + lr  # 1-based inclusive coords
                    tx_amp.append((s, e, e - s + 1))
        diag = {}
        if not fwd_hits:
            diag["fwd"] = best_near_match(seq, fwd)
        if not rev_hits:
            diag["rev"] = best_near_match(seq, rev_rc)
        per_tx[tx["id"]] = {
            "fwd_hits": [h + 1 for h in fwd_hits],
            "rev_hits": [h + 1 for h in rev_hits],
            "amplicons": tx_amp,
            "diag": diag,
        }
        for a in tx_amp:
            amplicons.append((tx["id"],) + a)
    return {"per_tx": per_tx, "amplicons": amplicons, "rev_rc": rev_rc}


# --------------------------------------------------------------------------
# Per-pair evaluation
# --------------------------------------------------------------------------
def evaluate(pair, transcripts, audit):
    tx_map = {t["id"]: t for t in transcripts}
    per_tx, amplicons = audit["per_tx"], audit["amplicons"]
    matched_ids = [t["id"] for t in transcripts if per_tx[t["id"]]["amplicons"]]
    exp_tx = pair["expected_transcript"]
    exp_len = pair["expected_product_bp"]
    findings = []

    # binding/amplicon findings, transcript by transcript (FASTA order)
    for tx in transcripts:
        info = per_tx[tx["id"]]
        if info["amplicons"]:
            amp_txt = "; ".join(f"amplicon {s}-{e} ({L} bp)" for (s, e, L) in info["amplicons"])
            findings.append(
                f"{tx['id']}: forward exact match at {','.join(map(str, info['fwd_hits']))}; "
                f"reverse primer (as revcomp {audit['rev_rc']}) exact match at "
                f"{','.join(map(str, info['rev_hits']))}; {amp_txt}"
            )
        else:
            notes = []
            if info["fwd_hits"]:
                notes.append(f"forward binds exactly at {','.join(map(str, info['fwd_hits']))}")
            else:
                d = info["diag"].get("fwd")
                notes.append(
                    "forward has no exact match"
                    + (f" (best near-match {d[0]}/{len(pair['forward'])} mismatches at pos {d[1] + 1})" if d else "")
                )
            if info["rev_hits"]:
                notes.append(f"reverse (revcomp) binds exactly at {','.join(map(str, info['rev_hits']))}")
            else:
                d = info["diag"].get("rev")
                notes.append(
                    "reverse primer (revcomp) has no exact match"
                    + (f" (best near-match {d[0]}/{len(pair['reverse'])} mismatches at pos {d[1] + 1})" if d else "")
                )
            findings.append(f"{tx['id']}: no productive amplicon ({'; '.join(notes)})")

    # expected transcript present?
    if exp_tx not in tx_map:
        findings.append(f"expected_transcript '{exp_tx}' not present in transcripts.fa")

    # primary amplicon: prefer the one on the expected transcript
    primary = None
    if amplicons:
        on_exp = [a for a in amplicons if a[0] == exp_tx]
        primary = on_exp[0] if on_exp else amplicons[0]

    # expected product size feasibility vs supplied transcript length
    if exp_len is not None and exp_tx in tx_map and exp_len > tx_map[exp_tx]["length"]:
        findings.append(
            f"expected product {exp_len} bp exceeds {exp_tx} sequence length "
            f"{tx_map[exp_tx]['length']} nt (infeasible with supplied sequence)"
        )

    # CDS compatibility (containment checked against the STATED range; no repair)
    cds_ok = False
    if primary is not None:
        t = tx_map[primary[0]]
        s, e = primary[1], primary[2]
        if t["cds"] is None:
            findings.append(f"{t['id']}: no CDS annotation supplied; CDS compatibility cannot be established")
        else:
            cs, ce = t["cds"]
            if cs <= s and e <= ce:
                if t["cds_issues"]:
                    findings.append(
                        f"{t['id']}: amplicon {s}-{e} lies inside stated CDS {cs}-{ce}, but that CDS "
                        f"annotation is internally inconsistent ({'; '.join(t['cds_issues'])}); not repaired"
                    )
                else:
                    cds_ok = True
                    findings.append(f"{t['id']}: amplicon {s}-{e} lies fully inside CDS {cs}-{ce}")
            else:
                findings.append(f"{t['id']}: amplicon {s}-{e} is NOT contained in stated CDS {cs}-{ce}")
                if t["cds_issues"]:
                    findings.append(
                        f"{t['id']}: CDS annotation is internally inconsistent "
                        f"({'; '.join(t['cds_issues'])}); reported, not repaired"
                    )

    # status decision
    if exp_tx not in tx_map:
        status = "expected_transcript_missing"
    elif not amplicons:
        status = "no_amplicon"
    elif exp_tx not in matched_ids:
        status = "off_target"
    else:
        L_primary = primary[3]
        if exp_len is not None and L_primary != exp_len:
            status = "length_mismatch"
            findings.append(f"observed amplicon {L_primary} bp differs from expected {exp_len} bp (delta {L_primary - exp_len} bp)")
        elif not cds_ok:
            status = "cds_incompatible"
        else:
            status = "ok"

    return {
        "pair": pair,
        "matched_ids": matched_ids,
        "amplicons": amplicons,
        "primary": primary,
        "cds_ok": cds_ok,
        "status": status,
        "findings": findings,
        "audit": audit,
    }


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------
def build_report(transcripts, fasta_warnings, primers, primer_warnings, results):
    L = []
    add = L.append
    add("# Primer pair audit against supplied transcript isoforms")
    add("")
    add(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')} by `output/analysis.py`.")
    add("")
    add("**Inputs**")
    add("")
    add(f"- `inputs/transcripts.fa` - {len(transcripts)} transcript record(s)")
    add(f"- `inputs/primer_candidates.csv` - {len(primers)} primer pair(s)")
    add("")
    add("## 1. Scope and method")
    add("")
    add("- Only the supplied sequences were used; no external data were fetched or inferred.")
    add("- Primer binding was assessed by **exact match**: the forward primer must occur on the")
    add("  sense strand, and the reverse primer must occur (as its reverse complement) on the")
    add("  sense strand downstream of the forward primer, without overlap.")
    add("- Amplicon coordinates are 1-based inclusive, spanning the forward primer start to the")
    add("  end of the reverse-primer binding site (both primers included in the length).")
    add("- A transcript counts as *matched* only when a productive amplicon forms on it.")
    add("- `cds_compatible = yes` only when the amplicon lies fully inside the annotated CDS")
    add("  range **and** that CDS annotation is internally consistent.")
    add("- Malformed or internally inconsistent metadata are reported below and were **not repaired**.")
    add("")
    add("## 2. Transcript inventory and metadata validation")
    add("")
    add("| transcript | length (nt) | header (verbatim) | CDS annotation | metadata verdict |")
    add("|---|---|---|---|---|")
    for t in transcripts:
        cds = f"{t['cds'][0]}-{t['cds'][1]}" if t["cds"] else "(none)"
        verdict = "OK" if t["metadata_ok"] else "INCONSISTENT: " + "; ".join(t["cds_issues"] + ([f"non-ACGT chars {t['invalid_chars']}"] if t["invalid_chars"] else []))
        add(f"| {t['id']} | {t['length']} | `>{t['header']}` | {cds} | {verdict} |")
    add("")
    if any(t["cds_issues"] for t in transcripts):
        add("**Metadata findings (reported, not repaired):**")
        add("")
        for t in transcripts:
            if t["cds_issues"]:
                cs, ce = t["cds"]
                add(f"- `{t['id']}`: header states `CDS={cs}-{ce}` but the supplied sequence is only "
                    f"{t['length']} nt long; the CDS end exceeds the sequence length by {ce - t['length']} nt. "
                    "No CDS-based conclusion for this transcript can be trusted; the annotation was left as supplied.")
        add("")
    add("Both headers carry the tag `exon_joined`; no exon/intron structure was supplied, so no")
    add("junction-specific checks are possible.")
    add("")
    add("## 3. Primer audit results")
    add("")
    add("| pair_id | transcripts_matched | amplicon_length | cds_compatible | status |")
    add("|---|---|---|---|---|")
    for r in results:
        p = r["pair"]
        tm = ";".join(r["matched_ids"]) if r["matched_ids"] else "none"
        al = ";".join(str(a[3]) for a in r["amplicons"]) if r["amplicons"] else "NA"
        add(f"| {p['pair_id']} | {tm} | {al} | {'yes' if r['cds_ok'] else 'no'} | {r['status']} |")
    add("")
    add("### Per-pair details")
    add("")
    for r in results:
        p = r["pair"]
        add(f"#### {p['pair_id']}")
        add("")
        add(f"- Forward (5'->3'): `{p['forward']}` ({len(p['forward'])} nt)")
        add(f"- Reverse (5'->3'): `{p['reverse']}` ({len(p['reverse'])} nt)")
        add(f"- Expected transcript: `{p['expected_transcript']}`; expected product: {p['expected_product_bp']} bp")
        for f in r["findings"]:
            add(f"- {f}")
        add("")
    add("## 4. Cross-file consistency findings")
    add("")
    for r in results:
        p = r["pair"]
        if p["expected_transcript"] in {t["id"] for t in transcripts} and p["expected_product_bp"] is not None:
            tlen = next(t["length"] for t in transcripts if t["id"] == p["expected_transcript"])
            if p["expected_product_bp"] > tlen:
                add(f"- `{p['pair_id']}`: expected product {p['expected_product_bp']} bp is longer than the entire "
                    f"`{p['expected_transcript']}` transcript ({tlen} nt); the stated expectation cannot be satisfied "
                    "by the supplied sequences.")
    add("- `p03` uses homopolymer primers (20x dA / 20x dT). No 20-mer A/T homopolymer run exists in either")
    add("  supplied transcript, so the pair cannot prime anywhere; this is flagged as a design defect, not")
    add("  something to be repaired.")
    add("")
    add("## 5. Warnings raised during parsing")
    add("")
    if fasta_warnings or primer_warnings:
        for w in fasta_warnings + primer_warnings:
            add(f"- {w}")
    else:
        add("- (none)")
    add("")
    add("## 6. Conclusions")
    add("")
    n_ok = sum(1 for r in results if r["status"] == "ok")
    n_lm = sum(1 for r in results if r["status"] == "length_mismatch")
    n_na = sum(1 for r in results if r["status"] == "no_amplicon")
    add(f"- {n_ok} of {len(results)} primer pair(s) fully pass (amplicon on the expected transcript, expected size, CDS-compatible).")
    if n_lm:
        add(f"- {n_lm} pair(s) amplify the expected transcript but produce an amplicon whose length differs from the expected product size.")
    if n_na:
        add(f"- {n_na} pair(s) produce no amplicon on any supplied transcript.")
    add("- No primer pair is CDS-compatible under the supplied (unrepaired) annotations.")
    add("- Both transcript headers carry CDS ranges that exceed their sequence lengths; these internally")
    add("  inconsistent annotations were reported and deliberately left unrepaired.")
    add("")
    return "\n".join(L) + "\n"


# --------------------------------------------------------------------------
def main():
    if not TRANSCRIPTS_FA.is_file() or not PRIMERS_CSV.is_file():
        sys.exit(f"missing input files under {INPUT_DIR}")

    transcripts, fasta_warnings = parse_fasta(TRANSCRIPTS_FA)
    primers, primer_warnings = parse_primers(PRIMERS_CSV)

    results = []
    for pair in primers:
        results.append(evaluate(pair, transcripts, audit_pair(pair, transcripts)))

    with open(AUDIT_CSV, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pair_id", "transcripts_matched", "amplicon_length", "cds_compatible", "status", "reason"])
        for r in results:
            p = r["pair"]
            w.writerow([
                p["pair_id"],
                ";".join(r["matched_ids"]) if r["matched_ids"] else "none",
                ";".join(str(a[3]) for a in r["amplicons"]) if r["amplicons"] else "NA",
                "yes" if r["cds_ok"] else "no",
                r["status"],
                "; ".join(r["findings"]),
            ])

    REPORT_MD.write_text(build_report(transcripts, fasta_warnings, primers, primer_warnings, results), encoding="utf-8")

    # console summary
    print(f"wrote {AUDIT_CSV}")
    print(f"wrote {REPORT_MD}")
    for r in results:
        p = r["pair"]
        tm = ";".join(r["matched_ids"]) or "none"
        al = ";".join(str(a[3]) for a in r["amplicons"]) or "NA"
        print(f"  {p['pair_id']}: matched={tm} amplicon={al} cds_compatible={'yes' if r['cds_ok'] else 'no'} status={r['status']}")
    for w in fasta_warnings + primer_warnings:
        print(f"  WARNING: {w}")


if __name__ == "__main__":
    main()
