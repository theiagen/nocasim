import math

import numpy as np

from nocasim.config import SimConfig
from nocasim.fragment import Fragment

GC_OPT = 0.47
GC_BIAS_K = 6.0


def gc_capture_prob(
    gc: float, gc_opt: float = GC_OPT, gc_bias_k: float = GC_BIAS_K
) -> float:
    return math.exp(-gc_bias_k * (gc - gc_opt) ** 2)


def capture_probability(frag: Fragment) -> float:
    if not frag.overlaps_vp1:
        return 0.0
    return gc_capture_prob(frag.gc)


def run_capture(
    viral_frags: list[Fragment],
    bg_frags: list[Fragment],
    config: SimConfig,
    rng: np.random.Generator,
) -> tuple[list[Fragment], list[Fragment]]:
    on_target = []
    for frag in viral_frags:
        p = capture_probability(frag)
        if rng.random() < p:
            on_target.append(frag)

    if not on_target:
        return on_target, []

    n_on = len(on_target)
    n_off = int(n_on * config.off_target_rate / (1.0 - config.off_target_rate))
    n_off = min(n_off, len(bg_frags))

    if n_off > 0 and bg_frags:
        indices = rng.choice(len(bg_frags), size=n_off, replace=len(bg_frags) < n_off)
        off_target = [bg_frags[i] for i in indices]
    else:
        off_target = []

    return on_target, off_target
