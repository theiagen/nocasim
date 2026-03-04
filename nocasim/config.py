from dataclasses import dataclass
from pathlib import Path


@dataclass
class SimConfig:
    references_dir: Path
    art_modern_bin: Path
    outdir: Path

    probes_fasta: Path | None = None
    human_bg_fasta: Path | None = None
    microbiome_bg_fasta: Path | None = None
    wastewater_bg_fasta: Path | None = None

    sample_type: str = "stool"

    read_len: int = 150
    dup_rate: float = 0.40
    off_target_rate: float = 0.592

    fragment_mean: int = 380
    fragment_sd: int = 80
    fragment_min: int = 200
    fragment_max: int = 600

    seed: int = 42

    vp1_start: int = 5_100
    vp1_end: int = 6_800
