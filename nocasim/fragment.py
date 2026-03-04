from dataclasses import dataclass

import numpy as np
from scipy.stats import truncnorm

from nocasim.config import SimConfig
from nocasim.genome import GenomeRecord, gc_content, reverse_complement


@dataclass
class Fragment:
    id: str
    sequence: str
    start: int
    end: int
    strand: str
    gc: float
    source: str
    overlaps_vp1: bool


def _truncnorm_rvs(
    mean: int, sd: int, low: int, high: int, size: int, rng: np.random.Generator
) -> np.ndarray:
    a = (low - mean) / sd
    b = (high - mean) / sd
    return truncnorm.rvs(a, b, loc=mean, scale=sd, size=size, random_state=rng).astype(int)


def sample_fragments(
    genome: GenomeRecord,
    n: int,
    config: SimConfig,
    source: str,
    vp1_start: int,
    vp1_end: int,
    rng: np.random.Generator,
) -> list[Fragment]:
    if n == 0:
        return []

    lengths = _truncnorm_rvs(
        config.fragment_mean, config.fragment_sd,
        config.fragment_min, config.fragment_max,
        size=n, rng=rng,
    )
    lengths = np.clip(lengths, config.fragment_min, min(config.fragment_max, genome.length))

    max_starts = genome.length - lengths
    max_starts = np.maximum(max_starts, 0)
    starts = (rng.random(n) * (max_starts + 1)).astype(int)
    ends = starts + lengths

    strands = rng.choice(["+", "-"], size=n)

    fragments = []
    for i in range(n):
        s, e = int(starts[i]), int(ends[i])
        seq = genome.sequence[s:e]
        strand = strands[i]
        if strand == "-":
            seq = reverse_complement(seq)

        overlaps = s < vp1_end and e > vp1_start

        frag = Fragment(
            id=f"frag_{i:04d}_{source}_{s}_{e}_{strand}",
            sequence=seq,
            start=s,
            end=e,
            strand=strand,
            gc=gc_content(seq),
            source=source,
            overlaps_vp1=overlaps,
        )
        fragments.append(frag)

    return fragments
