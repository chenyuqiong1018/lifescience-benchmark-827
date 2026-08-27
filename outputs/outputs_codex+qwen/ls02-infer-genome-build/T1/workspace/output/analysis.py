#!/usr/bin/env python3
"""Infer the genome build (hg18 / hg19 / hg38 / T2T) of a chr20 VCF.

Method
------
For every variant record in the VCF we fetch the reference base(s) at the
declared 1-based POS from each supplied chromosome-20 FASTA and compare them
with the VCF REF allele (case-insensitive).  A VCF whose coordinates belong
to a given build must have REF == reference at (virtually) every site, so
the build whose reference reproduces the REF alleles is the correct one.

Chromosome naming ("20" vs "chr20") is recorded but deliberately NOT used
as evidence for the build call.  T2T/hs1 has no supplied reference file in
this workspace, so it is explicitly excluded and never counted as a
mismatch (per inputs/references/README.md).

Only the Python standard library is used, so the analysis is fully
reproducible:  python output/analysis.py
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
VCF_PATH = ROOT / "inputs" / "vcf.infer.build.q1.vcf.gz"
REF_DIR = ROOT / "inputs" / "references"
MANIFEST_PATH = REF_DIR / "reference_manifest.json"
OUT_JSON = ROOT / "output" / "build_call.json"

BUILDS = ["hg18", "hg19", "hg38"]
REF_FILES = {b: REF_DIR / f"{b}_chr20.fa.gz" for b in BUILDS}
VALID_REF = re.compile(r"^[ACGTNacgtn]+$")
ACCEPTED_CHROM = {"20", "chr20"}


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_fasta(path: Path) -> tuple[str, list[str]]:
    """Return (uppercase sequence, list of record names) from a gz FASTA."""
    seq_parts: list[str] = []
    names: list[str] = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                names.append(line[1:].split()[0])
            else:
                seq_parts.append(line.strip())
    return "".join(seq_parts).upper(), names


def iter_vcf_records(path: Path):
    """Yield (chrom, pos, ref, alt, raw_line) for each data line."""
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            yield line.rstrip("\n")


def main() -> int:
    # ------------------------------------------------------------------
    # 1. Reference integrity + loading
    # ------------------------------------------------------------------
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    manifest_sha = {
        e["file"]: e["sha256"] for e in manifest.get("files", [])
    }

    ref_integrity = {}
    sequences: dict[str, str] = {}
    headers: dict[str, list[str]] = {}
    for build, path in REF_FILES.items():
        if not path.exists():
            ref_integrity[build] = "MISSING_FILE"
            continue
        digest = sha256_of(path)
        expected = manifest_sha.get(path.name)
        ref_integrity[build] = {
            "sha256": digest,
            "manifest_sha256": expected,
            "verified": bool(expected and digest == expected),
        }
        seq, names = load_fasta(path)
        sequences[build] = seq
        headers[build] = names

    # ------------------------------------------------------------------
    # 2. Scan the VCF once and check every REF allele against every build
    # ------------------------------------------------------------------
    stats = {
        b: Counter(
            matches=0,
            mismatches=0,
            uninformative=0,   # N in REF or reference slice
            out_of_range=0,    # POS/REF extends beyond the contig
        )
        for b in BUILDS
    }
    chrom_names: Counter = Counter()
    n_records = 0
    n_skipped = 0
    min_pos: int | None = None
    max_pos: int | None = None
    mismatch_examples: dict[str, list[dict]] = {b: [] for b in BUILDS}
    unique_support: Counter = Counter()  # variants matching exactly one build

    for line in iter_vcf_records(VCF_PATH):
        fields = line.split("\t")
        if len(fields) < 5:
            n_skipped += 1
            continue
        chrom, pos_s, _id, ref, alt = fields[:5]
        try:
            pos = int(pos_s)
        except ValueError:
            n_skipped += 1
            continue
        if chrom not in ACCEPTED_CHROM:
            n_skipped += 1
            continue
        if not VALID_REF.match(ref) or pos < 1:
            n_skipped += 1
            continue

        n_records += 1
        chrom_names[chrom] += 1
        min_pos = pos if min_pos is None else min(min_pos, pos)
        max_pos = pos if max_pos is None else max(max_pos, pos)

        matching_builds: list[str] = []
        for build in BUILDS:
            seq = sequences.get(build)
            if seq is None:
                continue
            start0 = pos - 1
            end0 = start0 + len(ref)
            if end0 > len(seq):
                stats[build]["out_of_range"] += 1
                continue
            ref_slice = seq[start0:end0]
            if "N" in ref_slice or "N" in ref.upper():
                stats[build]["uninformative"] += 1
                continue
            if ref.upper() == ref_slice:
                stats[build]["matches"] += 1
                matching_builds.append(build)
            else:
                stats[build]["mismatches"] += 1
                if len(mismatch_examples[build]) < 5:
                    mismatch_examples[build].append(
                        {
                            "pos": pos,
                            "vcf_ref": ref,
                            "reference_base": ref_slice,
                        }
                    )
        if len(matching_builds) == 1:
            unique_support[matching_builds[0]] += 1

    # ------------------------------------------------------------------
    # 3. Score builds and make the call
    # ------------------------------------------------------------------
    per_build = {}
    for b in BUILDS:
        c = stats[b]
        checked = c["matches"] + c["mismatches"]
        rate = c["matches"] / checked if checked else 0.0
        per_build[b] = {
            "contig_length": len(sequences.get(b, "")),
            "fasta_header": headers.get(b, [])[:1],
            "ref_matches": c["matches"],
            "ref_mismatches": c["mismatches"],
            "uninformative": c["uninformative"],
            "out_of_range": c["out_of_range"],
            "variants_checked": checked,
            "ref_match_rate": round(rate, 6),
            "unique_support_variants": unique_support[b],
            "first_mismatch_examples": mismatch_examples[b],
        }

    ranked = sorted(
        (b for b in BUILDS if per_build[b]["variants_checked"] > 0),
        key=lambda b: (
            per_build[b]["ref_match_rate"],
            per_build[b]["ref_matches"],
        ),
        reverse=True,
    )
    best = ranked[0] if ranked else None
    best_rate = per_build[best]["ref_match_rate"] if best else 0.0
    runner_rate = per_build[ranked[1]]["ref_match_rate"] if len(ranked) > 1 else 0.0
    margin = best_rate - runner_rate

    if best is None:
        build, confidence = "unknown", "none"
    elif best_rate >= 0.99 and margin >= 0.05:
        build, confidence = best, "high"
    elif best_rate >= 0.95 and margin >= 0.02:
        build, confidence = best, "medium"
    else:
        build, confidence = "ambiguous", "low"

    # T2T: no reference supplied -> explicitly excluded, never a mismatch.
    t2t_status = {
        "candidate": "T2T (hs1)",
        "status": "excluded_no_reference_supplied",
        "note": (
            "No T2T chr20 FASTA is present under inputs/references; per the "
            "references README it is excluded from scoring rather than "
            "treated as a mismatch."
        ),
    }

    called = per_build.get(build, {})
    result = {
        "build": build,
        "confidence": confidence,
        "n_variants_checked": n_records,
        "n_ref_matches": called.get("ref_matches", 0),
        "n_ref_mismatches": called.get("ref_mismatches", 0),
        "evidence": {
            "method": (
                "Each VCF REF allele was fetched at its declared 1-based POS "
                "from every supplied chr20 reference and compared "
                "case-insensitively; the build whose reference reproduces "
                "the REF alleles is called. Chromosome naming was recorded "
                "but not used as evidence."
            ),
            "vcf": {
                "file": "inputs/vcf.infer.build.q1.vcf.gz",
                "records_analyzed": n_records,
                "records_skipped": n_skipped,
                "chrom_field_values": dict(chrom_names),
                "pos_range": [min_pos, max_pos],
            },
            "per_build": per_build,
            "build_ranking": [
                {
                    "build": b,
                    "ref_match_rate": per_build[b]["ref_match_rate"],
                    "ref_matches": per_build[b]["ref_matches"],
                    "ref_mismatches": per_build[b]["ref_mismatches"],
                }
                for b in ranked
            ],
            "best_vs_runner_up_margin": round(margin, 6),
            "reference_integrity": ref_integrity,
            "t2t": t2t_status,
        },
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    # Console summary ---------------------------------------------------
    print(f"VCF records analyzed : {n_records} (skipped {n_skipped})")
    print(f"chrom field values   : {dict(chrom_names)}")
    print(f"POS range            : {min_pos}..{max_pos}")
    for b in BUILDS:
        p = per_build[b]
        print(
            f"{b}: len={p['contig_length']:,}  matches={p['ref_matches']:,}  "
            f"mismatches={p['ref_mismatches']:,}  uninformative={p['uninformative']:,}  "
            f"out_of_range={p['out_of_range']:,}  rate={p['ref_match_rate']:.4%}  "
            f"unique_support={p['unique_support_variants']:,}"
        )
    print(f"CALL: build={build} confidence={confidence}")
    print(f"Wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
