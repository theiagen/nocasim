from collections.abc import Iterator

import numpy as np

from nocasim.config import SimConfig
from nocasim.fragment import Fragment, sample_fragments
from nocasim.genome import GenomeRecord, gc_content

HUMAN_GC_MEAN = 0.41
HUMAN_GC_SD = 0.03
MICROBIOME_GC_MEAN = 0.52
MICROBIOME_GC_SD = 0.05
WASTEWATER_GC_MEAN = 0.48
WASTEWATER_GC_SD = 0.07
BATCH_SIZE = 50_000

SAMPLE_TYPE_MIX = {
    "stool": {"human": 0.80, "microbiome": 0.20, "wastewater": 0.00},
    "wastewater": {"human": 0.10, "microbiome": 0.30, "wastewater": 0.60},
}


def _synthetic_sequence(length: int, gc_target: float, rng: np.random.Generator) -> str:
    gc_target = max(0.01, min(0.99, gc_target))
    n_gc = int(length * gc_target)
    n_at = length - n_gc
    bases = (
        ["G"] * (n_gc // 2)
        + ["C"] * (n_gc - n_gc // 2)
        + ["A"] * (n_at // 2)
        + ["T"] * (n_at - n_at // 2)
    )
    rng.shuffle(bases)
    return "".join(bases)


def _sample_gc(
    n: int, gc_mean: float, gc_sd: float, rng: np.random.Generator
) -> np.ndarray:
    gc_values = rng.normal(gc_mean, gc_sd, size=n)
    return np.clip(gc_values, 0.30, 0.70)


def generate_background_fragments(
    n: int,
    config: SimConfig,
    human_ref: dict[str, GenomeRecord] | None,
    microbiome_ref: dict[str, GenomeRecord] | None,
    rng: np.random.Generator,
    wastewater_ref: dict[str, GenomeRecord] | None = None,
) -> Iterator[list[Fragment]]:
    mix = SAMPLE_TYPE_MIX.get(config.sample_type, SAMPLE_TYPE_MIX["stool"])
    n_human = int(n * mix["human"])
    n_wastewater = int(n * mix["wastewater"])
    n_microbiome = n - n_human - n_wastewater

    yield from _generate_source_fragments(
        n_human,
        "human_bg",
        human_ref,
        HUMAN_GC_MEAN,
        HUMAN_GC_SD,
        config,
        rng,
    )
    yield from _generate_source_fragments(
        n_microbiome,
        "microbiome_bg",
        microbiome_ref,
        MICROBIOME_GC_MEAN,
        MICROBIOME_GC_SD,
        config,
        rng,
    )
    if n_wastewater > 0:
        yield from _generate_source_fragments(
            n_wastewater,
            "wastewater_bg",
            wastewater_ref,
            WASTEWATER_GC_MEAN,
            WASTEWATER_GC_SD,
            config,
            rng,
        )


def _generate_source_fragments(
    n: int,
    source: str,
    ref: dict[str, GenomeRecord] | None,
    gc_mean: float,
    gc_sd: float,
    config: SimConfig,
    rng: np.random.Generator,
) -> Iterator[list[Fragment]]:
    generated = 0
    while generated < n:
        batch_n = min(BATCH_SIZE, n - generated)

        if ref:
            ref_names = list(ref.keys())
            chosen_name = rng.choice(ref_names)
            genome = ref[chosen_name]
            batch = sample_fragments(
                genome,
                batch_n,
                config,
                source,
                vp1_start=-1,
                vp1_end=-1,
                rng=rng,
            )
        else:
            gc_values = _sample_gc(batch_n, gc_mean, gc_sd, rng)
            from scipy.stats import truncnorm

            a = (config.fragment_min - config.fragment_mean) / config.fragment_sd
            b = (config.fragment_max - config.fragment_mean) / config.fragment_sd
            lengths = truncnorm.rvs(
                a,
                b,
                loc=config.fragment_mean,
                scale=config.fragment_sd,
                size=batch_n,
                random_state=rng,
            ).astype(int)

            batch = []
            for i in range(batch_n):
                idx = generated + i
                length = int(lengths[i])
                gc_val = float(gc_values[i])
                seq = _synthetic_sequence(length, gc_val, rng)
                frag = Fragment(
                    id=f"frag_{idx:06d}_{source}_0_{length}_+",
                    sequence=seq,
                    start=0,
                    end=length,
                    strand="+",
                    gc=gc_content(seq),
                    source=source,
                    overlaps_vp1=False,
                )
                batch.append(frag)

        generated += batch_n
        yield batch
