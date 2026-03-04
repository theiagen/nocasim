import numpy as np

from nocasim.capture import gc_capture_prob, capture_probability, run_capture
from nocasim.config import SimConfig
from nocasim.fragment import Fragment
from pathlib import Path


def _config() -> SimConfig:
    return SimConfig(
        references_dir=Path("."), art_modern_bin=Path("."), outdir=Path("."),
    )


def _frag(overlaps_vp1: bool, gc: float = 0.47, strand: str = "+") -> Fragment:
    return Fragment(
        id="test_frag", sequence="ACGT" * 50, start=0, end=200,
        strand=strand, gc=gc, source="viral", overlaps_vp1=overlaps_vp1,
    )


def test_gc_optimal_gives_high_prob():
    p = gc_capture_prob(0.47)
    assert p > 0.99


def test_gc_extreme_gives_low_prob():
    p = gc_capture_prob(0.10)
    assert p < 0.5


def test_non_overlapping_fragment_zero_prob():
    frag = _frag(overlaps_vp1=False)
    assert capture_probability(frag) == 0.0


def test_overlapping_fragment_nonzero_prob():
    frag = _frag(overlaps_vp1=True)
    assert capture_probability(frag) > 0.0


def test_run_capture_on_target_rate():
    config = _config()
    rng = np.random.default_rng(42)

    viral = [_frag(overlaps_vp1=True, gc=0.47) for _ in range(1000)]
    bg = [_frag(overlaps_vp1=False, gc=0.45) for _ in range(10000)]

    on, off = run_capture(viral, bg, config, rng)

    assert len(on) > 0
    assert len(off) > 0

    total = len(on) + len(off)
    on_rate = len(on) / total
    assert 0.30 < on_rate < 0.50, f"On-target rate {on_rate} outside expected range"
