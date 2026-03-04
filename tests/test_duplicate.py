import numpy as np

from nocasim.duplicate import assign_copy_counts
from nocasim.fragment import Fragment


def _frag(idx: int) -> Fragment:
    return Fragment(
        id=f"frag_{idx}", sequence="ACGT" * 50, start=0, end=200,
        strand="+", gc=0.5, source="viral", overlaps_vp1=True,
    )


def test_copy_counts_at_least_one():
    frags = [_frag(i) for i in range(100)]
    rng = np.random.default_rng(42)
    result = assign_copy_counts(frags, 0.40, rng)
    assert all(count >= 1 for _, count in result)


def test_mean_copies_near_expected():
    frags = [_frag(i) for i in range(10000)]
    rng = np.random.default_rng(42)
    result = assign_copy_counts(frags, 0.40, rng)
    mean_copies = sum(c for _, c in result) / len(result)
    expected = 1.0 / (1.0 - 0.40)
    assert abs(mean_copies - expected) < 0.1, f"Mean {mean_copies} vs expected {expected}"


def test_empty_input():
    result = assign_copy_counts([], 0.40, np.random.default_rng(42))
    assert result == []


def test_zero_dup_rate():
    frags = [_frag(i) for i in range(100)]
    rng = np.random.default_rng(42)
    result = assign_copy_counts(frags, 0.0, rng)
    assert all(count == 1 for _, count in result)
