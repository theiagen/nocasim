import numpy as np

from nocasim.background import (
    generate_background_fragments,
    _synthetic_sequence,
)
from nocasim.config import SimConfig
from pathlib import Path


def _config(sample_type: str = "stool") -> SimConfig:
    return SimConfig(
        references_dir=Path("."),
        art_modern_bin=Path("."),
        outdir=Path("."),
        sample_type=sample_type,
    )


def test_synthetic_sequence_length():
    rng = np.random.default_rng(42)
    seq = _synthetic_sequence(500, 0.45, rng)
    assert len(seq) == 500


def test_synthetic_sequence_gc_approx():
    rng = np.random.default_rng(42)
    seq = _synthetic_sequence(1000, 0.50, rng)
    gc = sum(1 for c in seq if c in "GC") / len(seq)
    assert 0.45 < gc < 0.55


def test_generate_background_count():
    config = _config()
    rng = np.random.default_rng(42)
    total = 0
    for batch in generate_background_fragments(1000, config, None, None, rng):
        total += len(batch)
        for frag in batch:
            assert not frag.overlaps_vp1
            assert frag.source in ("human_bg", "microbiome_bg")
    assert total == 1000


def test_background_gc_distribution():
    config = _config()
    rng = np.random.default_rng(42)
    all_gc = []
    for batch in generate_background_fragments(2000, config, None, None, rng):
        all_gc.extend(f.gc for f in batch)
    mean_gc = sum(all_gc) / len(all_gc)
    assert 0.35 < mean_gc < 0.55


def test_wastewater_background_sources():
    config = _config(sample_type="wastewater")
    rng = np.random.default_rng(42)
    sources = set()
    total = 0
    for batch in generate_background_fragments(1000, config, None, None, rng):
        total += len(batch)
        for frag in batch:
            sources.add(frag.source)
            assert not frag.overlaps_vp1
    assert total == 1000
    assert "wastewater_bg" in sources
    assert "human_bg" in sources
    assert "microbiome_bg" in sources


def test_wastewater_mix_fractions():
    config = _config(sample_type="wastewater")
    rng = np.random.default_rng(42)
    counts = {"human_bg": 0, "microbiome_bg": 0, "wastewater_bg": 0}
    for batch in generate_background_fragments(10_000, config, None, None, rng):
        for frag in batch:
            counts[frag.source] += 1
    total = sum(counts.values())
    assert counts["wastewater_bg"] / total > 0.5
    assert counts["human_bg"] / total < 0.2
    assert counts["microbiome_bg"] / total > 0.2


def test_stool_has_no_wastewater():
    config = _config(sample_type="stool")
    rng = np.random.default_rng(42)
    for batch in generate_background_fragments(1000, config, None, None, rng):
        for frag in batch:
            assert frag.source != "wastewater_bg"
