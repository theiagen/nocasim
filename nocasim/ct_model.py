from nocasim.config import SimConfig

BASE_FRACTION = 0.03
CALIBRATION_CT = 25.0
EFFICIENCY = 3.32


def ct_to_viral_fraction(ct: float, base_fraction: float = BASE_FRACTION,
                         calibration_ct: float = CALIBRATION_CT,
                         efficiency: float = EFFICIENCY) -> float:
    return base_fraction * 2 ** ((calibration_ct - ct) / efficiency)


def ct_to_expected_vp1_depth(ct: float, total_reads: int, config: SimConfig,
                             vp1_length: int = 1700) -> float:
    vf = ct_to_viral_fraction(ct)
    on_target_rate = 1.0 - config.off_target_rate
    viral_on_target = total_reads * vf * on_target_rate
    mean_frag_len = config.fragment_mean
    bases_on_vp1 = viral_on_target * mean_frag_len
    return bases_on_vp1 / vp1_length


def predict_completeness(ct: float, total_reads: int = 1_000_000,
                         config: SimConfig | None = None) -> str:
    if config is None:
        from nocasim.config import SimConfig
        from pathlib import Path
        config = SimConfig(
            references_dir=Path("."), art_modern_bin=Path("."), outdir=Path("."),
        )
    depth = ct_to_expected_vp1_depth(ct, total_reads, config)
    if depth >= 20:
        return "complete"
    elif depth >= 5:
        return "low_coverage"
    return "incomplete"


def viral_fragment_count(ct: float, target_vp1_depth: float,
                         config: SimConfig, vp1_length: int = 1700) -> int:
    vf = ct_to_viral_fraction(ct)
    on_target_rate = 1.0 - config.off_target_rate
    mean_frag_len = config.fragment_mean
    bases_needed = target_vp1_depth * vp1_length
    frags_on_target = bases_needed / mean_frag_len
    total_viral = frags_on_target / on_target_rate
    return max(int(total_viral / vf), 100)
