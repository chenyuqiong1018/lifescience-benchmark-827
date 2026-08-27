#!/usr/bin/env python3
"""Infer the genome build (hg18 / hg19 / hg38 / T2T) of a chr20 VCF.

Method
------
The call is based on reproducible allele/coordinate checks, not on
chromosome naming:

1. Verify integrity of the three supplied chr20 reference FASTAs against
   the SHA-256 values in inputs/references/reference_manifest.json.
2. Parse every variant record of inputs/vcf.infer.build.q1.vcf.gz
   (CHROM, POS, REF; 1-based VCF coordinates).
3. For each supplied build, fetch the reference base(s) at the declared
   POS and compare with the VCF REF allele (case-insensitive, full
   string compare so indels are checked too).
4. Count matches / mismatches / unverifiable (ref base 'N' or POS beyond
   contig length) per build. The build whose REF alleles reproduce the
   VCF essentially perfectly while all others fail is the call.
5. T2T/hs1 has no reference file in inputs/; per the supplied README it
   is explicitly excluded as a candidate (not counted as a mismatch).

Outputs: output/build_call.json and output/report.md.
Standard library only; deterministic given the inputs.
"""

import gzip
import hashlib
import json
import os

WS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VCF = os.path.join(WS, "inputs", "vcf.infer.build.q1.vcf.gz")
REFDIR = os.path.join(WS, "inputs", "references")
MANIFEST = os.path.join(REFDIR, "reference_manifest.json")
OUTDIR = os.path.join(WS, "output")

BUILDS = ["hg18", "hg19", "hg38"]          # T2T: no reference file supplied
COMP = str.maketrans("", "")               # identity (keep uppercase handling manual)


def sha256_of(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_fasta_seq(path: str):
    """Return (header, uppercase sequence) for a single-record FASTA.gz."""
    header = None
    parts = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    raise ValueError(f"multiple records in {path}")
                header = line[1:].split()[0]
            else:
                parts.append(line)
    if header is None:
        raise ValueError(f"no FASTA record in {path}")
    return header, "".join(parts).upper()


def parse_vcf(path: str):
    """Yield (chrom, pos, ref) for every variant record."""
    records = []
    chroms = set()
    header_lines = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                header_lines += 1
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                raise ValueError("malformed VCF line: " + line[:80])
            chrom, pos, ref = cols[0], int(cols[1]), cols[3]
            chroms.add(chrom)
            records.append((chrom, pos, ref))
    return records, sorted(chroms), header_lines


def check_build(seq: str, records):
    """Compare every VCF REF allele against seq at its 1-based POS."""
    n_match = n_mismatch = n_unverifiable = 0
    first_mismatches = []
    for chrom, pos, ref in records:
        start = pos - 1
        end = start + len(ref)
        if start < 0 or end > len(seq):
            n_unverifiable += 1
            continue
        refseg = seq[start:end]
        if "N" in refseg or "N" in ref.upper():
            n_unverifiable += 1
            continue
        if refseg == ref.upper():
            n_match += 1
        else:
            n_mismatch += 1
            if len(first_mismatches) < 5:
                first_mismatches.append(
                    {"chrom": chrom, "pos": pos, "vcf_ref": ref, "ref_genome_base": refseg}
                )
    return n_match, n_mismatch, n_unverifiable, first_mismatches


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    manifest = json.load(open(MANIFEST, "r", encoding="utf-8"))

    # 1. reference integrity + load
    refs = {}
    integrity = {}
    for entry in manifest["files"]:
        build = next(b for b in BUILDS if entry["file"].startswith(b))
        path = os.path.join(REFDIR, entry["file"])
        digest = sha256_of(path)
        ok = digest == entry["sha256"]
        integrity[build] = {"sha256_ok": ok, "source": entry["source"]}
        if not ok:
            raise SystemExit(f"reference integrity check failed for {entry['file']}")
        header, seq = load_fasta_seq(path)
        refs[build] = {"header": header, "seq": seq, "length": len(seq)}

    # 2. parse VCF
    records, chroms, header_lines = parse_vcf(VCF)
    if chroms != ["20"]:
        raise SystemExit(f"unexpected contigs in VCF: {chroms}")
    n_variants = len(records)
    min_pos = min(r[1] for r in records)
    max_pos = max(r[1] for r in records)
    n_snps = sum(1 for r in records if len(r[2]) == 1 and r[2] in "ACGT")
    n_indels = n_variants - n_snps

    # 3. REF allele verification per build
    stats = {}
    for build in BUILDS:
        m, mm, unv, examples = check_build(refs[build]["seq"], records)
        stats[build] = {
            "fasta_header": refs[build]["header"],
            "contig_length_bp": refs[build]["length"],
            "n_ref_matches": m,
            "n_ref_mismatches": mm,
            "n_unverifiable": unv,
            "match_rate": round(m / n_variants, 6),
            "first_mismatch_examples": examples,
        }

    # 4. decide: unique build with (near-)perfect REF reproduction
    ordered = sorted(BUILDS, key=lambda b: (stats[b]["n_ref_matches"],
                                            -stats[b]["n_ref_mismatches"]), reverse=True)
    best, second = ordered[0], ordered[1]
    best_s, second_s = stats[best], stats[second]
    dominant = (best_s["n_ref_mismatches"] == 0
                and best_s["n_unverifiable"] == 0
                and best_s["n_ref_matches"] == n_variants
                and second_s["match_rate"] < 0.90)
    build_call = best if dominant else best
    confidence = "high" if dominant else (
        "medium" if best_s["match_rate"] > 0.99 and best_s["match_rate"] > second_s["match_rate"] + 0.05
        else "low")

    t2t_note = ("T2T/hs1 excluded: no T2T chr20 reference is supplied in inputs/references "
                "(per its README, an unavailable candidate is excluded, not counted as a mismatch). "
                "A T2T coordinate call would additionally be incompatible with an exact, "
                "genome-wide REF match against one of the supplied linear-build references.")

    evidence = [
        f"VCF contains {n_variants} variant records on contig '{chroms[0]}' "
        f"(POS {min_pos}..{max_pos}; {n_snps} SNVs, {n_indels} indels/other); "
        "chromosome naming alone was not used as evidence.",
        f"All three supplied chr20 references passed SHA-256 verification against "
        f"reference_manifest.json (hg18 {refs['hg18']['length']} bp, hg19 {refs['hg19']['length']} bp, "
        f"hg38 {refs['hg38']['length']} bp).",
    ]
    for b in BUILDS:
        s = stats[b]
        evidence.append(
            f"{b}: {s['n_ref_matches']}/{n_variants} VCF REF alleles match the reference at the "
            f"declared 1-based coordinates ({s['match_rate'] * 100:.4f}%); "
            f"{s['n_ref_mismatches']} mismatches, {s['n_unverifiable']} unverifiable.")
    if dominant:
        evidence.append(
            f"{build_call} reproduces 100% of REF alleles while every alternative build mismatches "
            f"the majority of sites; the coordinate system is therefore identified unambiguously.")
    evidence.append(t2t_note)

    result = {
        "build": build_call,
        "confidence": confidence,
        "n_variants_checked": n_variants,
        "n_ref_matches": best_s["n_ref_matches"],
        "n_ref_mismatches": best_s["n_ref_mismatches"],
        "evidence": evidence,
        "per_build": stats,
        "reference_integrity": integrity,
        "vcf": {
            "file": "inputs/vcf.infer.build.q1.vcf.gz",
            "contigs": chroms,
            "min_pos": min_pos,
            "max_pos": max_pos,
            "n_snv": n_snps,
            "n_indel_or_other": n_indels,
        },
        "t2t": {"candidate": True, "reference_available": False,
                "status": "explicitly_excluded", "note": t2t_note},
    }
    with open(os.path.join(OUTDIR, "build_call.json"), "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=2)
        fh.write("\n")

    # 5. report
    lines = [
        "# Genome build inference report",
        "",
        f"**Input:** `inputs/vcf.infer.build.q1.vcf.gz` ({n_variants} chr20 variant records, "
        f"POS {min_pos}\u2013{max_pos}; {n_snps} SNVs, {n_indels} indels/other).",
        "",
        f"## Call: **{build_call}** (confidence: {confidence})",
        "",
        "The VCF REF allele of every record was compared against each supplied chr20 "
        "reference at the declared 1-based VCF coordinate (full-string comparison, so "
        "indels are verified as well). Chromosome naming was recorded but deliberately "
        "not used as evidence.",
        "",
        "| Build | Contig length (bp) | REF matches | Mismatches | Unverifiable | Match rate |",
        "|---|---|---|---|---|---|",
    ]
    for b in BUILDS:
        s = stats[b]
        lines.append(f"| {b} | {s['contig_length_bp']:,} | {s['n_ref_matches']:,} | "
                     f"{s['n_ref_mismatches']:,} | {s['n_unverifiable']:,} | "
                     f"{s['match_rate'] * 100:.4f}% |")
    lines += [
        f"| T2T (hs1) | \u2014 | \u2014 | \u2014 | \u2014 | excluded (no reference supplied) |",
        "",
        "## Interpretation",
        "",
    ]
    if dominant:
        lines += [
            f"- All {n_variants} REF alleles match **{build_call}** exactly; no other build "
            f"comes close ({second} matches only {second_s['match_rate'] * 100:.2f}% of sites).",
            "- A VCF aligned to a different build would show widespread REF mismatches at "
            "these same coordinates, so the coordinate system is unambiguous.",
        ]
    else:
        lines.append("- Match rates did not show a single dominant build; review per-build stats.")
    lines += [
        "",
        "Example REF mismatches under the wrong builds (VCF REF vs reference base at the same coordinate):",
        "",
    ]
    for b in BUILDS:
        if b == build_call:
            continue
        ex = stats[b]["first_mismatch_examples"][:3]
        shown = "; ".join(f"pos {e['pos']}: VCF={e['vcf_ref']} vs {b}={e['ref_genome_base']}" for e in ex)
        lines.append(f"- **{b}**: {shown}")
    lines += [
        "",
        "## T2T/hs1",
        "",
        "No T2T chr20 FASTA is present in `inputs/references` (its README records that the "
        "chromosome-only endpoint was unavailable). Following the supplied instructions, T2T is "
        "explicitly excluded as a candidate rather than counted as a mismatch. The exact "
        f"genome-wide REF agreement with {build_call} further argues against T2T coordinates.",
        "",
        "## Reproducibility",
        "",
        "- Reference integrity: SHA-256 of all three FASTAs matches `reference_manifest.json`.",
        "- Re-run: `python output/analysis.py` regenerates `output/build_call.json` and this report.",
        "- Only files under `./inputs` were read; only `./output` was written.",
        "",
    ]
    with open(os.path.join(OUTDIR, "report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))

    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "first_mismatch_examples"}
                      for k, v in stats.items()}, indent=2))
    print("CALL:", build_call, "| confidence:", confidence,
          f"| checked={n_variants} matches={best_s['n_ref_matches']} mismatches={best_s['n_ref_mismatches']}")


if __name__ == "__main__":
    main()
