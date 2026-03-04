from pathlib import Path

from nocasim.genome import (
    GenomeRecord,
    gc_content,
    load_fasta,
    reverse_complement,
    extract_vp1,
)

FIXTURES = Path(__file__).parent / "fixtures"


def test_load_fasta():
    records = load_fasta(FIXTURES / "test_reference.fasta")
    assert len(records) == 1
    name = next(iter(records))
    assert records[name].length > 0
    assert set(records[name].sequence) <= {"A", "C", "G", "T"}


def test_gc_content():
    assert gc_content("GGCC") == 1.0
    assert gc_content("AATT") == 0.0
    assert gc_content("ATGC") == 0.5
    assert gc_content("") == 0.0


def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("") == ""


def test_extract_vp1():
    genome = GenomeRecord(
        name="test", sequence="A" * 100 + "G" * 50 + "T" * 100, length=250
    )
    vp1 = extract_vp1(genome, 100, 150)
    assert vp1.sequence == "G" * 50
    assert vp1.length == 50
    assert vp1.name == "test_VP1"
