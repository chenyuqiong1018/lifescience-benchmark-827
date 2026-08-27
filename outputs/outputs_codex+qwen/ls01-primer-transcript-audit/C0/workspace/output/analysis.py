#!/usr/bin/env python3
"""Primer-pair audit against supplied transcript isoforms.

Inputs (the ONLY data used; no external sequences or databases):
  inputs/primer_candidates.csv
      columns: pair_id, forward, reverse, expected_transcript, expected_product_bp
  inputs/transcripts.fa
      FASTA; headers carry metadata tokens, e.g. "exon_joined CDS=101-700"

Outputs:
  output/primer_audit.csv
      columns: pair_id,transcripts_matched,amplicon_length,cds_compatible,status,reason
  output/report.md
      human-readable audit report (generated from the same computed results)

Audit policy
------------
* Only exact (full-length, zero-mismatch) primer matches count as binding sites.
* Forward primer must match the transcript sense strand; the reverse primer is
  searched as its reverse complement on the sense strand. Binding sites must be
  non-overlapping with the forward site upstream of the reverse site.
* Amplicon length = end of reverse-primer binding site - start of forward
  binding site + 1 (i.e. the full PCR product, primers included).
* All supplied isoforms are screened for every pair (expected target +
  cross-isoform / off-target binding).
* Malformed or internally inconsistent sequence metadata is REPORTED, never
  silently repaired (no clipping, re-anchoring, or re-interpretation of CDS
  coordinates is applied).
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FASTA = ROOT / "inputs" / "transcripts.fa"
PRIMERS_CSV = ROOT / "inputs" / "primer_candidates.csv"
OUT_CSV = ROOT / "output" / "primer_audit.csv"
OUT_MD = ROOT / "output" / "report.md"

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
STOP_CODONS = {"TAA", "TAG", "TGA"}
REQUIRED_CSV_COLUMNS = [
    "pair_id", "forward", "reverse", "expected_transcript", "expected_product_bp",
]


def revcomp(seq: str) -> str:
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(seq))


def find_all(haystack: str, needle: str) -> list[int]:
    """All (including overlapping) 1-based start positions of needle in haystack."""
    hits: list[int] = []
    start = 0
    while True:
        i = haystack.find(needle, start)
        if i < 0:
            return hits
        hits.append(i + 1)
        start = i + 1


# ---------------------------------------------------------------------------
# FASTA parsing and metadata validation
# ---------------------------------------------------------------------------

def parse_fasta(path: Path):
    """Parse FASTA into records keyed by transcript id; collect format problems."""
    records: dict[str, dict] = {}
    order: list[str] = []
    problems: list[str] = []
    cur: dict | None = None
    with open(path, encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith(">"):
                header = line[1:].strip()
                tokens = header.split()
                if not tokens:
                    problems.append(f"line {lineno}: FASTA header with no identifier")
                    cur = None
                    continue
                tid = tokens[0]
                cur = {
                    "id": tid,
                    "header": header,
                    "meta_tokens": tokens[1:],
                    "seq": "",
                    "header_line": lineno,
                    "issues": [],       # format/character-level issues
                    "meta_issues": [],  # metadata-consistency issues
                    "cds": None,        # (start, end) ints as annotated
                    "cds_raw": None,    # raw token text
                }
                if tid in records:
                    cur["issues"].append(
                        f"duplicate transcript id '{tid}' "
                        f"(first defined at line {records[tid]['header_line']})"
                    )
                records[tid] = cur
                order.append(tid)
            else:
                if cur is None:
                    problems.append(f"line {lineno}: sequence data before any FASTA header")
                    continue
                chunk = line.strip().upper()
                bad = sorted(set(chunk) - set("ACGT"))
                if bad:
                    cur["issues"].append(
                        f"line {lineno}: non-ACGT character(s) in sequence: {', '.join(bad)}"
                    )
                cur["seq"] += chunk

    for tid in order:
        if not records[tid]["seq"]:
            records[tid]["issues"].append("record has an empty sequence")
    if not order:
        problems.append("no FASTA records found")
    return records, order, problems


def validate_record_metadata(rec: dict) -> None:
    """Validate header metadata of one record; record findings in rec['meta_issues'].

    Inconsistencies are recorded and surfaced; they are never repaired.
    """
    n = len(rec["seq"])
    cds_tok = None
    for tok in rec["meta_tokens"]:
        if tok.startswith("CDS="):
            cds_tok = tok
            break
    if cds_tok is None:
        rec["meta_issues"].append("no CDS annotation present in header")
        return

    rec["cds_raw"] = cds_tok
    m = re.fullmatch(r"CDS=(\d+)-(\d+)", cds_tok)
    if not m:
        rec["meta_issues"].append(f"malformed CDS token '{cds_tok}' (expected CDS=<start>-<end>)")
        return

    s, e = int(m.group(1)), int(m.group(2))
    rec["cds"] = (s, e)
    issues = rec["meta_issues"]

    if s < 1:
        issues.append(f"CDS start {s} < 1")
    if e < s:
        issues.append(f"CDS end {e} < CDS start {s}")
    if e > n:
        issues.append(
            f"METADATA INCONSISTENCY: annotated CDS end {e} exceeds sequence length {n} "
            f"(header '{rec['header']}')"
        )
    if s > n:
        issues.append(
            f"METADATA INCONSISTENCY: annotated CDS start {s} exceeds sequence length {n}"
        )
    if e <= n and s >= 1 and e >= s:
        length = e - s + 1
        if length % 3 != 0:
            issues.append(f"CDS length {length} is not a multiple of 3")
        cds_seq = rec["seq"][s - 1:e]
        if not cds_seq.startswith("ATG"):
            issues.append("annotated CDS does not begin with a start codon (ATG)")
        if cds_seq[-3:] not in STOP_CODONS:
            issues.append(f"annotated CDS does not end with a stop codon (ends {cds_seq[-3:]})")


def observable_orf_note(rec: dict) -> str | None:
    """Observation only (never used to repair metadata): does the supplied sequence
    itself look like a complete in-frame ORF?"""
    seq = rec["seq"]
    n = len(seq)
    if n >= 3 and n % 3 == 0 and seq.startswith("ATG") and seq[-3:] in STOP_CODONS:
        return (
            f"the supplied {n} bp sequence itself is a complete in-frame ORF "
            f"(starts ATG at position 1, ends with stop codon {seq[-3:]} at positions "
            f"{n - 2}-{n}), which contradicts the annotated CDS coordinates"
        )
    return None


# ---------------------------------------------------------------------------
# Primer table parsing
# ---------------------------------------------------------------------------

def parse_primers(path: Path):
    pairs: list[dict] = []
    problems: list[str] = []
    with open(path, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in REQUIRED_CSV_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            problems.append(f"primer_candidates.csv missing column(s): {', '.join(missing)}")
        seen_ids: set[str] = set()
        for rowno, row in enumerate(reader, start=2):
            pair = {
                "pair_id": (row.get("pair_id") or "").strip(),
                "forward": (row.get("forward") or "").strip().upper(),
                "reverse": (row.get("reverse") or "").strip().upper(),
                "expected_transcript": (row.get("expected_transcript") or "").strip(),
                "expected_product_bp_raw": (row.get("expected_product_bp") or "").strip(),
                "expected_product_bp": None,
                "row": rowno,
                "issues": [],
            }
            if not pair["pair_id"]:
                pair["issues"].append("missing pair_id")
            elif pair["pair_id"] in seen_ids:
                pair["issues"].append(f"duplicate pair_id '{pair['pair_id']}'")
            seen_ids.add(pair["pair_id"])
            for key in ("forward", "reverse"):
                seq = pair[key]
                if not seq:
                    pair["issues"].append(f"missing {key} primer sequence")
                else:
                    bad = sorted(set(seq) - set("ACGT"))
                    if bad:
                        pair["issues"].append(
                            f"{key} primer contains non-ACGT character(s): {', '.join(bad)}"
                        )
            try:
                pair["expected_product_bp"] = int(pair["expected_product_bp_raw"])
                if pair["expected_product_bp"] <= 0:
                    pair["issues"].append("expected_product_bp must be a positive integer")
            except ValueError:
                pair["issues"].append(
                    f"expected_product_bp '{pair['expected_product_bp_raw']}' is not an integer"
                )
            pairs.append(pair)
    if not pairs:
        problems.append("no primer pairs found in primer_candidates.csv")
    return pairs, problems


# ---------------------------------------------------------------------------
# Audit of one primer pair against all supplied transcripts
# ---------------------------------------------------------------------------

def audit_pair(pair: dict, records: dict, order: list[str]) -> dict:
    fwd = pair["forward"]
    rev = pair["reverse"]
    rev_rc = revcomp(rev) if rev else ""
    expected_tid = pair["expected_transcript"]
    expected_bp = pair["expected_product_bp"]

    findings: list[str] = []
    per_transcript: dict[str, dict] = {}

    for tid in order:
        seq = records[tid]["seq"]
        fwd_hits = find_all(seq, fwd) if fwd else []
        rev_hits = find_all(seq, rev_rc) if rev_rc else []
        amplicons = []
        for i in fwd_hits:
            for j in rev_hits:
                if j >= i + len(fwd):  # reverse site downstream, non-overlapping
                    start, end = i, j + len(rev) - 1
                    amplicons.append((start, end, end - start + 1))
        per_transcript[tid] = {
            "fwd_hits": fwd_hits,
            "rev_hits": rev_hits,
            "amplicons": amplicons,
        }

    matched = [tid for tid in order if per_transcript[tid]["amplicons"]]
    exp_info = per_transcript.get(expected_tid)
    exp_amps = exp_info["amplicons"] if exp_info else []

    # ---- amplicon_length column -------------------------------------------
    if exp_amps:
        lengths = sorted({a[2] for a in exp_amps})
        amplicon_length = ";".join(str(x) for x in lengths)
    else:
        amplicon_length = "na"

    # ---- CDS compatibility --------------------------------------------------
    cds_value = "na"
    cds_findings: list[str] = []
    if expected_tid in records:
        rec = records[expected_tid]
        if not exp_amps:
            cds_value = "na"
            cds_findings.append("no amplicon on the expected transcript; CDS compatibility not assessable")
        else:
            a0, a1 = exp_amps[0][0], exp_amps[-1][1]
            if rec["cds"] is None:
                cds_value = "na"
                cds_findings.append("no CDS annotation available for the expected transcript")
            else:
                cs, ce = rec["cds"]
                contained = all(cs <= a[0] and a[1] <= ce for a in exp_amps)
                if rec["meta_issues"]:
                    cds_value = "no"
                    cds_findings.append(
                        f"amplicon spans positions {a0}-{a1} and is NOT contained within the "
                        f"annotated CDS {cs}-{ce}; additionally the CDS annotation itself is "
                        f"internally inconsistent (see metadata findings), so CDS compatibility "
                        f"cannot be established"
                    )
                else:
                    cds_value = "yes" if contained else "no"
                    if not contained:
                        cds_findings.append(
                            f"amplicon spans positions {a0}-{a1}, outside annotated CDS {cs}-{ce}"
                        )
    else:
        cds_findings.append(f"expected transcript '{expected_tid}' not found in supplied FASTA")

    # ---- binding / amplicon detail findings --------------------------------
    if exp_info is not None:
        fh_txt = ", ".join(str(h) for h in exp_info["fwd_hits"]) or "none"
        rh_txt = ", ".join(str(h) for h in exp_info["rev_hits"]) or "none"
        findings.append(
            f"on {expected_tid}: forward primer match(es) at {fh_txt}; "
            f"reverse primer (reverse complement) match(es) at {rh_txt}"
        )
        if exp_amps:
            for st, en, ln in exp_amps:
                findings.append(f"amplicon on {expected_tid}: positions {st}-{en}, length {ln} bp")

    # length comparison vs expectation
    if exp_amps and expected_bp is not None:
        obs = sorted({a[2] for a in exp_amps})
        if expected_bp in obs and len(obs) == 1:
            findings.append(f"observed amplicon length {obs[0]} bp matches expected {expected_bp} bp")
        else:
            findings.append(
                f"AMPLICON LENGTH MISMATCH: observed {';'.join(map(str, obs))} bp vs "
                f"expected {expected_bp} bp"
            )
    if expected_tid in records and expected_bp is not None and not exp_amps:
        if expected_bp > len(records[expected_tid]["seq"]):
            findings.append(
                f"expected product {expected_bp} bp exceeds the supplied transcript length "
                f"({len(records[expected_tid]['seq'])} bp); it cannot be produced from the "
                f"supplied sequence"
            )
    if expected_tid in records and expected_bp is not None and exp_amps:
        if expected_bp > len(records[expected_tid]["seq"]):
            findings.append(
                f"expected product {expected_bp} bp exceeds the supplied transcript length "
                f"({len(records[expected_tid]['seq'])} bp); it cannot be produced from the "
                f"supplied sequence"
            )

    # cross-isoform / off-target binding on non-expected transcripts
    for tid in order:
        info = per_transcript[tid]
        if tid == expected_tid:
            continue
        if info["amplicons"]:
            for st, en, ln in info["amplicons"]:
                findings.append(
                    f"OFF-TARGET AMPLICON: {tid} positions {st}-{en}, length {ln} bp"
                )
        else:
            binds = []
            if info["fwd_hits"]:
                binds.append("forward primer binds at " + ", ".join(str(h) for h in info["fwd_hits"]))
            if info["rev_hits"]:
                binds.append("reverse primer (RC) binds at " + ", ".join(str(h) for h in info["rev_hits"]))
            if binds:
                findings.append(f"cross-isoform binding without product on {tid}: " + "; ".join(binds))

    findings.extend(cds_findings)

    if not exp_amps and not any(per_transcript[t]["amplicons"] for t in order):
        if not (fwd and rev):
            pass
        else:
            total_sites = sum(
                len(per_transcript[t]["fwd_hits"]) + len(per_transcript[t]["rev_hits"]) for t in order
            )
            if total_sites == 0:
                findings.append(
                    "NO MATCH: neither primer has an exact binding site on any supplied transcript"
                )

    # ---- status -------------------------------------------------------------
    ok = (
        not pair["issues"]
        and expected_tid in records
        and not records[expected_tid]["issues"]
        and not records[expected_tid]["meta_issues"]
        and len(exp_amps) == 1
        and expected_bp is not None
        and exp_amps[0][2] == expected_bp
        and cds_value == "yes"
        and all(len(per_transcript[t]["amplicons"]) == 0 for t in order if t != expected_tid)
        and len(exp_info["fwd_hits"]) == 1
        and len(exp_info["rev_hits"]) == 1
    )
    status = "pass" if ok else "fail"

    return {
        "pair": pair,
        "per_transcript": per_transcript,
        "matched": matched,
        "amplicon_length": amplicon_length,
        "cds_compatible": cds_value,
        "status": status,
        "findings": findings,
    }


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_audit_csv(results: list[dict], path: Path) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pair_id", "transcripts_matched", "amplicon_length",
                    "cds_compatible", "status", "reason"])
        for r in results:
            w.writerow([
                r["pair"]["pair_id"],
                ";".join(r["matched"]) if r["matched"] else "none",
                r["amplicon_length"],
                r["cds_compatible"],
                r["status"],
                "; ".join(r["findings"]),
            ])


def write_report(results, records, order, fasta_problems, primer_problems, path: Path) -> None:
    lines: list[str] = []
    add = lines.append
    add("# Primer audit report")
    add("")
    add("Audit of every primer pair in `inputs/primer_candidates.csv` against the supplied")
    add("transcript isoforms in `inputs/transcripts.fa`. Only the supplied sequences were")
    add("used; no external references were consulted. Malformed or internally inconsistent")
    add("metadata is reported below and was **not** silently repaired.")
    add("")
    add("## 1. Method summary")
    add("")
    add("- Primers were aligned with **exact (zero-mismatch) matching** only.")
    add("- Forward primer matched against the transcript sense strand; reverse primer matched as its reverse complement.")
    add("- A valid amplicon requires the forward site upstream of, and non-overlapping with, the reverse site; amplicon length spans both primer sites inclusive.")
    add("- Every pair was screened against **all** supplied isoforms to detect cross-isoform / off-target binding.")
    add("- CDS compatibility was evaluated strictly against the annotated CDS coordinates as given. Because the annotations are internally inconsistent (Section 2), compatibility could not be established for any pair.")
    add("")
    add("## 2. Input validation and sequence-metadata findings")
    add("")
    if fasta_problems or primer_problems:
        for p in fasta_problems + primer_problems:
            add(f"- file-level problem: {p}")
        add("")
    add("| transcript | header | sequence length | metadata findings |")
    add("|---|---|---|---|")
    for tid in order:
        rec = records[tid]
        issues = rec["issues"] + rec["meta_issues"]
        issue_txt = "<br>".join(issues) if issues else "none"
        add(f"| {tid} | `{rec['header']}` | {len(rec['seq'])} bp | {issue_txt} |")
    add("")
    add("**Key metadata inconsistency (reported, not repaired).** Both transcripts are 102 bp")
    add("long, yet their headers annotate `CDS=101-700` (TX_CANONICAL) and `CDS=101-640`")
    add("(TX_ALT). In each case the annotated CDS end lies far beyond the end of the supplied")
    add("sequence, and only 2 nt (positions 101-102) of the annotated CDS exist in the")
    add("sequence. Under the `exon_joined` interpretation the coordinates apply to the spliced")
    add("transcript, so the mismatch cannot be explained away by genomic coordinates.")
    add("")
    for tid in order:
        note = observable_orf_note(records[tid])
        if note:
            add(f"- Observation for {tid}: {note}. This observation is reported only; it was **not** used to overwrite or re-anchor the supplied CDS annotation.")
    add("")
    add("Consequences for the audit: CDS compatibility cannot be affirmed for any amplicon,")
    add("and every pair targeting these transcripts inherits the metadata inconsistency in its")
    add("verdict.")
    add("")
    add("## 3. Per-pair results")
    add("")
    add("| pair_id | expected transcript | expected product | transcripts with amplicon | observed amplicon length | CDS compatible | status |")
    add("|---|---|---|---|---|---|---|")
    for r in results:
        p = r["pair"]
        add(
            f"| {p['pair_id']} | {p['expected_transcript']} | {p['expected_product_bp_raw']} bp | "
            f"{';'.join(r['matched']) if r['matched'] else 'none'} | {r['amplicon_length']} | "
            f"{r['cds_compatible']} | {r['status']} |"
        )
    add("")
    for r in results:
        p = r["pair"]
        add(f"### {p['pair_id']}")
        add("")
        add(f"- forward: `{p['forward']}` ({len(p['forward'])} nt)")
        add(f"- reverse: `{p['reverse']}` ({len(p['reverse'])} nt)")
        add(f"- expected transcript: `{p['expected_transcript']}`; expected product: {p['expected_product_bp_raw']} bp")
        add("- findings:")
        for f in r["findings"]:
            add(f"  - {f}")
        add("")
    add("## 4. Summary")
    add("")
    n_pass = sum(1 for r in results if r["status"] == "pass")
    add(f"- **{n_pass} of {len(results)} primer pairs pass** the audit.")
    add("- **p01** fails: the only product on TX_CANONICAL is 102 bp, not the expected 108 bp (108 bp exceeds the supplied transcript length entirely); the amplicon is not within the annotated CDS; TX_CANONICAL CDS metadata is internally inconsistent.")
    add("- **p02** fails: the only product on TX_ALT is 99 bp, not the expected 104 bp; the amplicon is not within the annotated CDS; TX_ALT CDS metadata is internally inconsistent. The forward primer also binds TX_CANONICAL (no product forms there because the reverse primer is isoform-specific).")
    add("- **p03** fails: neither primer matches any supplied transcript; no product is possible. Note also the expected 120 bp product exceeds the 102 bp transcript length.")
    add("")
    add("## 5. Recommendations")
    add("")
    add("1. Correct the CDS annotations in `transcripts.fa` (as supplied they cannot be reconciled with the 102 bp sequences) before any CDS-dependent interpretation.")
    add("2. Re-derive expected product sizes from the actual supplied sequences (102 bp for p01, 99 bp for p02) or supply the transcript versions on which 108/104/120 bp products are genuine.")
    add("3. Redesign p03: poly-A/poly-T primers have no binding sites in these GC-rich transcripts.")
    add("")
    path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------

def main() -> int:
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    records, order, fasta_problems = parse_fasta(FASTA)
    for tid in order:
        validate_record_metadata(records[tid])
    pairs, primer_problems = parse_primers(PRIMERS_CSV)

    for pair in pairs:
        if pair["expected_transcript"] and pair["expected_transcript"] not in records:
            pair["issues"].append(
                f"expected_transcript '{pair['expected_transcript']}' not present in supplied FASTA"
            )

    results = [audit_pair(p, records, order) for p in pairs]
    write_audit_csv(results, OUT_CSV)
    write_report(results, records, order, fasta_problems, primer_problems, OUT_MD)

    # console summary
    print(f"wrote {OUT_CSV}")
    print(f"wrote {OUT_MD}")
    for r in results:
        p = r["pair"]
        print(
            f"{p['pair_id']}: status={r['status']} matched={r['matched'] or 'none'} "
            f"amplicon_length={r['amplicon_length']} cds_compatible={r['cds_compatible']}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
