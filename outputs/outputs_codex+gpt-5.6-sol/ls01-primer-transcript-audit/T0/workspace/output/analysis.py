from __future__ import annotations

import csv
import re
from pathlib import Path


OUT = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[6]
DATA = ROOT / "inputs" / "ls01-primer-transcript-audit"
DNA = frozenset("ACGT")
HEADER = re.compile(r"^(\S+)\s+exon_joined\s+CDS=(\d+)-(\d+)$")


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def motif_starts(seq: str, motif: str) -> list[int]:
    return [i for i in range(len(seq)) if seq.startswith(motif, i)]


def read_fasta(path: Path) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    raw_header: str | None = None
    pieces: list[str] = []

    def finish() -> None:
        if raw_header is None:
            return
        sequence = "".join(pieces).upper()
        if set(sequence) - DNA:
            raise ValueError(f"non-DNA symbol in {raw_header}")
        parsed = HEADER.fullmatch(raw_header)
        if not parsed:
            name = raw_header.split()[0]
            records[name] = {
                "sequence": sequence,
                "metadata_error": f"malformed header: {raw_header}",
            }
            return
        name, start_text, end_text = parsed.groups()
        start, end = int(start_text), int(end_text)
        error = None
        if not (1 <= start <= end <= len(sequence)):
            error = f"CDS={start}-{end} exceeds transcript length {len(sequence)}"
        records[name] = {"sequence": sequence, "metadata_error": error}

    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith(">"):
            finish()
            raw_header = line[1:]
            pieces = []
        elif raw_header is None:
            raise ValueError("FASTA sequence precedes header")
        else:
            pieces.append(line)
    finish()
    return records


def main() -> None:
    with (DATA / "primer_candidates.csv").open(newline="", encoding="utf-8") as handle:
        primers = list(csv.DictReader(handle))
    transcripts = read_fasta(DATA / "transcripts.fa")
    rows: list[dict[str, str]] = []

    for primer in primers:
        forward, reverse = primer["forward"].upper(), primer["reverse"].upper()
        if len(forward) != 20 or len(reverse) != 20 or set(forward + reverse) - DNA:
            raise ValueError(f"invalid primer sequence in {primer['pair_id']}")
        reverse_site = reverse_complement(reverse)
        matches: list[tuple[str, int]] = []
        for name, record in transcripts.items():
            sequence = str(record["sequence"])
            for left in motif_starts(sequence, forward):
                for right in motif_starts(sequence, reverse_site):
                    if right >= left + len(forward):
                        matches.append((name, right + len(reverse_site) - left))

        names = sorted({name for name, _ in matches})
        lengths = sorted({length for _, length in matches})
        expected_name = primer["expected_transcript"]
        expected_length = int(primer["expected_product_bp"])
        reasons: list[str] = []
        if not matches:
            reasons.append("no valid primer-pair amplicon")
        if names and names != [expected_name]:
            reasons.append(f"matched {','.join(names)} rather than {expected_name}")
        if len(lengths) > 1:
            reasons.append(f"ambiguous amplicon lengths: {','.join(map(str, lengths))}")
        elif lengths and lengths[0] != expected_length:
            reasons.append(f"observed {lengths[0]} bp, expected {expected_length} bp")

        relevant = names or ([expected_name] if expected_name in transcripts else [])
        metadata_errors = [
            f"{name}: {transcripts[name]['metadata_error']}"
            for name in relevant
            if transcripts[name]["metadata_error"]
        ]
        reasons.extend(metadata_errors)
        rows.append(
            {
                "pair_id": primer["pair_id"],
                "transcripts_matched": ";".join(names),
                "amplicon_length": ";".join(map(str, lengths)),
                "cds_compatible": "NA" if not matches else ("false" if metadata_errors else "true"),
                "status": "fail" if reasons else "pass",
                "reason": "; ".join(reasons),
            }
        )

    fields = ["pair_id", "transcripts_matched", "amplicon_length", "cds_compatible", "status", "reason"]
    with (OUT / "primer_audit.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
