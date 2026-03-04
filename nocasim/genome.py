from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


@dataclass
class GenomeRecord:
    name: str
    sequence: str
    length: int


def load_fasta(path: Path) -> dict[str, GenomeRecord]:
    records = {}
    for rec in SeqIO.parse(str(path), "fasta"):
        seq = str(rec.seq).upper()
        records[rec.id] = GenomeRecord(name=rec.id, sequence=seq, length=len(seq))
    return records


def extract_vp1(genome: GenomeRecord, vp1_start: int, vp1_end: int) -> GenomeRecord:
    subseq = genome.sequence[vp1_start:vp1_end]
    return GenomeRecord(
        name=f"{genome.name}_VP1",
        sequence=subseq,
        length=len(subseq),
    )


def gc_content(sequence: str) -> float:
    if not sequence:
        return 0.0
    upper = sequence.upper()
    gc = sum(1 for c in upper if c in "GC")
    return gc / len(upper)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]
