import numpy as np

from nocasim.fragment import Fragment


def assign_copy_counts(
    fragments: list[Fragment],
    dup_rate: float,
    rng: np.random.Generator,
) -> list[tuple[Fragment, int]]:
    if not fragments:
        return []

    p = 1.0 - dup_rate
    counts = rng.geometric(p, size=len(fragments))
    return [(frag, int(count)) for frag, count in zip(fragments, counts)]
