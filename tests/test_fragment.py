from pathlib import Path

import numpy as np

from nocasim.config import SimConfig
from nocasim.fragment import sample_fragments
from nocasim.genome import load_fasta

FIXTURES = Path(__file__).parent / "fixtures"


def _make_config() -> SimConfig:
    return SimConfig(
        references_dir=FIXTURES,
        art_modern_bin=Path("/usr/bin/false"),
        outdir=Path("/tmp/test_out"),
    )


def test_fragment_count():
    config = _make_config()
    records = load_fasta(FIXTURES / "test_reference.fasta")
    genome = next(iter(records.values()))
    rng = np.random.default_rng(42)

    frags = sample_fragments(genome, 100, config, "viral", 0, genome.length, rng)
    assert len(frags) == 100


def test_fragment_lengths_within_bounds():
    config = _make_config()
    records = load_fasta(FIXTURES / "test_reference.fasta")
    genome = next(iter(records.values()))
    rng = np.random.default_rng(42)

    frags = sample_fragments(genome, 500, config, "viral", 0, genome.length, rng)
    lengths = [f.end - f.start for f in frags]
    for length in lengths:
        assert config.fragment_min <= length <= min(config.fragment_max, genome.length)


def test_fragment_strands():
    config = _make_config()
    records = load_fasta(FIXTURES / "test_reference.fasta")
    genome = next(iter(records.values()))
    rng = np.random.default_rng(42)

    frags = sample_fragments(genome, 1000, config, "viral", 0, genome.length, rng)
    strands = [f.strand for f in frags]
    plus_count = strands.count("+")
    assert 400 < plus_count < 600, (
        f"Expected ~50/50 strand split, got {plus_count}/1000"
    )


def test_fragment_vp1_overlap():
    config = _make_config()
    records = load_fasta(FIXTURES / "test_reference.fasta")
    genome = next(iter(records.values()))
    rng = np.random.default_rng(42)

    frags = sample_fragments(genome, 100, config, "viral", 500, 1000, rng)
    overlapping = [f for f in frags if f.overlaps_vp1]
    assert len(overlapping) > 0
    for f in overlapping:
        assert f.start < 1000 and f.end > 500


def test_fragment_ids_unique():
    config = _make_config()
    records = load_fasta(FIXTURES / "test_reference.fasta")
    genome = next(iter(records.values()))
    rng = np.random.default_rng(42)

    frags = sample_fragments(genome, 100, config, "viral", 0, genome.length, rng)
    ids = [f.id for f in frags]
    assert len(ids) == len(set(ids))
