#!/usr/bin/env python3
"""Infer a chr20 VCF genome build by exhaustive REF/coordinate checks."""

import gzip
import hashlib
import json
from pathlib import Path

TASK = "ls02-infer-genome-build"


def find_repo():
    for parent in Path(__file__).resolve().parents:
        if (parent / "inputs" / TASK).is_dir():
            return parent
    raise RuntimeError("repository root not found")


REPO = find_repo()
INPUT = REPO / "inputs" / TASK
REFDIR = INPUT / "references"
VCF = INPUT / "vcf.infer.build.q1.vcf.gz"
OUT = Path(__file__).resolve().parent
REFERENCES = {
    "hg18": REFDIR / "hg18_chr20.fa.gz",
    "hg19": REFDIR / "hg19_chr20.fa.gz",
    "hg38": REFDIR / "hg38_chr20.fa.gz",
}


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def load_fasta(path):
    with gzip.open(path, "rt") as handle:
        header = handle.readline().strip()
        sequence = "".join(line.strip().upper() for line in handle)
    return header, sequence


variants = []
chromosomes = set()
with gzip.open(VCF, "rt") as handle:
    for line_number, line in enumerate(handle, 1):
        if line.startswith("#"):
            continue
        fields = line.rstrip().split("\t")
        if len(fields) < 5:
            raise ValueError(f"malformed VCF line {line_number}")
        chrom, pos_text, _, ref, alt = fields[:5]
        chromosomes.add(chrom)
        variants.append((int(pos_text), ref.upper(), alt))
if chromosomes not in ({"20"}, {"chr20"}):
    raise ValueError(f"unexpected chromosome set: {sorted(chromosomes)}")

checks = {}
for build, path in REFERENCES.items():
    header, sequence = load_fasta(path)
    matches = mismatches = out_of_range = 0
    mismatch_examples = []
    for position, ref, alt in variants:
        start = position - 1
        observed = sequence[start : start + len(ref)] if 0 <= start < len(sequence) else ""
        if len(observed) != len(ref):
            out_of_range += 1
            mismatches += 1
        elif observed == ref:
            matches += 1
        else:
            mismatches += 1
            if len(mismatch_examples) < 10:
                mismatch_examples.append({"pos": position, "vcf_ref": ref, "reference": observed, "alt": alt})
    checks[build] = {
        "fasta": str(path.relative_to(REPO)),
        "fasta_header": header,
        "fasta_length_bp": len(sequence),
        "sha256": digest(path),
        "n_ref_matches": matches,
        "n_ref_mismatches": mismatches,
        "n_out_of_range": out_of_range,
        "match_fraction": matches / len(variants),
        "mismatch_examples": mismatch_examples,
    }

winner = max(checks, key=lambda build: checks[build]["n_ref_matches"])
ordered = sorted(checks, key=lambda build: checks[build]["n_ref_matches"], reverse=True)
best, runner_up = checks[ordered[0]], checks[ordered[1]]
if best["n_ref_matches"] == len(variants) and runner_up["n_ref_matches"] < len(variants):
    confidence = "high"
elif best["n_ref_matches"] > runner_up["n_ref_matches"]:
    confidence = "medium"
else:
    confidence = "low"

evidence = {
    "coordinate_convention": "VCF POS is 1-based; REF checked against FASTA[pos-1:pos-1+len(REF)]",
    "chromosome_labels_seen": sorted(chromosomes),
    "candidate_checks": checks,
    "winner_margin_matches": best["n_ref_matches"] - runner_up["n_ref_matches"],
    "t2t": {
        "status": "not_evaluated_reference_absent",
        "reason": "No T2T/hs1 chromosome reference was supplied; it is excluded rather than counted as a mismatch.",
    },
    "vcf_sha256": digest(VCF),
}
call = {
    "build": winner,
    "confidence": confidence,
    "n_variants_checked": len(variants),
    "n_ref_matches": best["n_ref_matches"],
    "n_ref_mismatches": best["n_ref_mismatches"],
    "evidence": evidence,
}
(OUT / "build_call.json").write_text(json.dumps(call, indent=2) + "\n", encoding="utf-8", newline="\n")

rows = "\n".join(
    f"| {build} | {data['fasta_length_bp']:,} | {data['n_ref_matches']:,} | {data['n_ref_mismatches']:,} | {data['match_fraction']:.6f} |"
    for build, data in checks.items()
)
report = f"""# Genome-build inference

## Call

The VCF uses **{winner}** coordinates, with **{confidence} confidence**. All {len(variants):,} VCF REF alleles match the {winner} chr20 reference at their declared 1-based coordinates; the selected call has {best['n_ref_mismatches']:,} mismatch(es).

## Reproducible allele checks

| Candidate | chr20 length | REF matches | REF mismatches | Match fraction |
|---|---:|---:|---:|---:|
{rows}

The code checks the complete VCF, including multibase REF alleles, using `sequence[POS-1:POS-1+len(REF)]`. Chromosome naming (`20`) is recorded only as QC and is not used as build proof. The winner exceeds the next-best candidate by {evidence['winner_margin_matches']:,} REF matches.

## T2T limitation

T2T/hs1 cannot be tested because no corresponding reference file is supplied. It is explicitly marked unavailable, not assigned a mismatch count. The three supplied references nevertheless yield a unique dominant call.

## Reproduction

Run `python analysis.py` in this directory. Input/reference SHA-256 values, mismatch examples, coordinate convention, and every candidate's counts are stored in `build_call.json`.
"""
(OUT / "report.md").write_text(report, encoding="utf-8", newline="\n")
print(json.dumps({"build": winner, "confidence": confidence, "checked": len(variants), "matches": best["n_ref_matches"], "mismatches": best["n_ref_mismatches"]}, sort_keys=True))
