from __future__ import annotations

import csv
import re
from pathlib import Path


OUT = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[6]
DATA = ROOT / "inputs" / "ls01-primer-transcript-audit"
HEADER_RE = re.compile(r"^(\S+)\s+exon_joined\s+CDS=(\d+)-(\d+)$")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def parse_fasta(path: Path) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    header: str | None = None
    sequence_parts: list[str] = []

    def store() -> None:
        if header is None:
            return
        match = HEADER_RE.fullmatch(header)
        if match is None:
            records[header.split()[0]] = {
                "sequence": "".join(sequence_parts),
                "cds": None,
                "metadata_error": f"malformed FASTA header: {header}",
            }
            return
        name, start_text, end_text = match.groups()
        sequence = "".join(sequence_parts)
        start, end = int(start_text), int(end_text)
        error = None
        if not (1 <= start <= end <= len(sequence)):
            error = (
                f"CDS={start}-{end} is outside transcript length {len(sequence)}"
            )
        records[name] = {
            "sequence": sequence,
            "cds": (start, end),
            "metadata_error": error,
        }

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line.startswith(">"):
            store()
            header = line[1:]
            sequence_parts = []
        elif header is None:
            raise ValueError("sequence data encountered before FASTA header")
        else:
            sequence_parts.append(line.upper())
    store()
    return records


def all_starts(sequence: str, motif: str) -> list[int]:
    return [index for index in range(len(sequence)) if sequence.startswith(motif, index)]


def main() -> None:
    with (DATA / "primer_candidates.csv").open(newline="", encoding="utf-8") as handle:
        primers = list(csv.DictReader(handle))
    transcripts = parse_fasta(DATA / "transcripts.fa")
    rows: list[dict[str, str]] = []

    for primer in primers:
        pair_id = primer["pair_id"]
        forward = primer["forward"].upper()
        reverse_site = reverse_complement(primer["reverse"].upper())
        matches: list[tuple[str, int]] = []
        reasons: list[str] = []

        for name, record in transcripts.items():
            sequence = str(record["sequence"])
            for forward_start in all_starts(sequence, forward):
                for reverse_start in all_starts(sequence, reverse_site):
                    if reverse_start >= forward_start + len(forward):
                        length = reverse_start + len(reverse_site) - forward_start
                        matches.append((name, length))

        matched_names = sorted({name for name, _ in matches})
        lengths = sorted({length for _, length in matches})
        expected_name = primer["expected_transcript"]
        expected_length = int(primer["expected_product_bp"])

        if not matches:
            reasons.append("no valid forward/reverse primer amplicon found")
        if matched_names and matched_names != [expected_name]:
            reasons.append(
                f"matched transcript(s) {','.join(matched_names)} differ from expected {expected_name}"
            )
        if len(lengths) > 1:
            reasons.append("multiple distinct amplicon lengths found")
        elif lengths and lengths[0] != expected_length:
            reasons.append(
                f"observed amplicon {lengths[0]} bp differs from expected {expected_length} bp"
            )

        relevant_names = matched_names or ([expected_name] if expected_name in transcripts else [])
        metadata_errors = [
            f"{name}: {transcripts[name]['metadata_error']}"
            for name in relevant_names
            if transcripts[name]["metadata_error"]
        ]
        reasons.extend(metadata_errors)

        cds_compatible = "NA"
        if matches:
            cds_compatible = "false" if metadata_errors else "true"
        status = "pass" if not reasons else "fail"
        rows.append(
            {
                "pair_id": pair_id,
                "transcripts_matched": ";".join(matched_names),
                "amplicon_length": ";".join(map(str, lengths)),
                "cds_compatible": cds_compatible,
                "status": status,
                "reason": "; ".join(reasons),
            }
        )

    fields = [
        "pair_id",
        "transcripts_matched",
        "amplicon_length",
        "cds_compatible",
        "status",
        "reason",
    ]
    with (OUT / "primer_audit.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
