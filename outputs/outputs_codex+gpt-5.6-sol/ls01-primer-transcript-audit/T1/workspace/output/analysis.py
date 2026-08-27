from __future__ import annotations

import csv
import re
from pathlib import Path


OUT = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[6]
DATA = ROOT / "inputs" / "ls01-primer-transcript-audit"
DNA = set("ACGT")
HEADER_RE = re.compile(r"^(\S+)\s+exon_joined\s+CDS=(\d+)-(\d+)$")


def rc(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def starts(sequence: str, motif: str) -> list[int]:
    return [i for i in range(len(sequence)) if sequence.startswith(motif, i)]


def fasta_records(path: Path) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    header: str | None = None
    chunks: list[str] = []

    def emit() -> None:
        if header is None:
            return
        sequence = "".join(chunks).upper()
        if set(sequence) - DNA:
            raise ValueError(f"invalid nucleotide in {header}")
        match = HEADER_RE.fullmatch(header)
        if not match:
            records[header.split()[0]] = {
                "sequence": sequence,
                "metadata_error": f"malformed header: {header}",
            }
            return
        name, start_raw, end_raw = match.groups()
        cds_start, cds_end = int(start_raw), int(end_raw)
        error = None
        if not (1 <= cds_start <= cds_end <= len(sequence)):
            error = f"CDS={cds_start}-{cds_end} outside sequence length {len(sequence)}"
        records[name] = {"sequence": sequence, "metadata_error": error}

    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith(">"):
            emit()
            header, chunks = line[1:], []
        elif header is None:
            raise ValueError("FASTA sequence before header")
        else:
            chunks.append(line)
    emit()
    return records


def main() -> None:
    with (DATA / "primer_candidates.csv").open(newline="", encoding="utf-8") as handle:
        primers = list(csv.DictReader(handle))
    transcripts = fasta_records(DATA / "transcripts.fa")
    output: list[dict[str, str]] = []

    for primer in primers:
        forward = primer["forward"].upper()
        reverse = primer["reverse"].upper()
        if len(forward) != 20 or len(reverse) != 20 or set(forward + reverse) - DNA:
            raise ValueError(f"malformed primer pair {primer['pair_id']}")
        reverse_site = rc(reverse)
        amplicons: list[tuple[str, int]] = []
        for name, record in transcripts.items():
            sequence = str(record["sequence"])
            for left in starts(sequence, forward):
                for right in starts(sequence, reverse_site):
                    if right >= left + len(forward):
                        amplicons.append((name, right + len(reverse_site) - left))

        names = sorted({name for name, _ in amplicons})
        lengths = sorted({length for _, length in amplicons})
        expected_name = primer["expected_transcript"]
        expected_bp = int(primer["expected_product_bp"])
        reason: list[str] = []
        if not amplicons:
            reason.append("no valid forward/reverse amplicon")
        if names and names != [expected_name]:
            reason.append(f"matched {','.join(names)} instead of {expected_name}")
        if len(lengths) > 1:
            reason.append(f"multiple product lengths {','.join(map(str, lengths))}")
        elif lengths and lengths[0] != expected_bp:
            reason.append(f"observed {lengths[0]} bp versus expected {expected_bp} bp")

        relevant = names or ([expected_name] if expected_name in transcripts else [])
        metadata_errors = [
            f"{name}: {transcripts[name]['metadata_error']}"
            for name in relevant
            if transcripts[name]["metadata_error"]
        ]
        reason.extend(metadata_errors)
        output.append(
            {
                "pair_id": primer["pair_id"],
                "transcripts_matched": ";".join(names),
                "amplicon_length": ";".join(map(str, lengths)),
                "cds_compatible": "NA" if not amplicons else ("false" if metadata_errors else "true"),
                "status": "fail" if reason else "pass",
                "reason": "; ".join(reason),
            }
        )

    fields = ["pair_id", "transcripts_matched", "amplicon_length", "cds_compatible", "status", "reason"]
    with (OUT / "primer_audit.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()
