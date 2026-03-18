from pathlib import Path

import numpy as np
import pytest

from nocasim.config import SimConfig
from nocasim.genome import load_fasta
from nocasim.mixture import (
    LineageSpec,
    MixtureSpec,
    generate_mixture_fragments,
    parse_mixture,
    resolve_mixture,
    resolve_preset,
    PRESETS,
)

FIXTURES = Path(__file__).parent / "fixtures"


def _make_config() -> SimConfig:
    return SimConfig(
        references_dir=FIXTURES,
        art_modern_bin=Path("/usr/bin/false"),
        outdir=Path("/tmp/test_out"),
    )


def _make_genome_records() -> dict:
    records = load_fasta(FIXTURES / "test_reference.fasta")
    genome = next(iter(records.values()))
    return {"GII.4": genome, "GII.17": genome}


def test_parse_bare_genotype():
    result = parse_mixture("GII.4")
    assert result == MixtureSpec(lineages=[LineageSpec("GII.4", 1.0)])


def test_parse_two_lineages():
    result = parse_mixture("GII.17:0.75,GII.4:0.25")
    assert len(result.lineages) == 2
    assert result.lineages[0] == LineageSpec("GII.17", 0.75)
    assert result.lineages[1] == LineageSpec("GII.4", 0.25)


def test_parse_three_lineages_normalizes():
    # Sum is 1.005, within tolerance 0.01, should normalize to exactly 1.0
    result = parse_mixture("GII.17:0.74,GII.4:0.16,GI.1:0.105")
    total = sum(spec.abundance for spec in result.lineages)
    assert abs(total - 1.0) < 1e-9


def test_parse_bad_sum_rejects():
    with pytest.raises(ValueError, match="sum"):
        parse_mixture("GII.17:0.5,GII.4:0.3")


def test_parse_duplicate_genotype():
    with pytest.raises(ValueError, match="uplicate"):
        parse_mixture("GII.4:0.5,GII.4:0.5")


def test_parse_negative_abundance():
    with pytest.raises(ValueError):
        parse_mixture("GII.4:-0.5,GII.17:1.5")


def test_parse_empty_string():
    with pytest.raises(ValueError):
        parse_mixture("")


def test_resolve_preset_us_2024():
    result = resolve_preset("us-2024")
    assert len(result.lineages) == 6
    assert result.lineages[0].genotype == "GII.17"
    total = sum(spec.abundance for spec in result.lineages)
    assert abs(total - 1.0) < 1e-9


def test_resolve_preset_all_defined():
    for name in PRESETS:
        result = resolve_preset(name)
        total = sum(spec.abundance for spec in result.lineages)
        assert abs(total - 1.0) < 1e-9


def test_resolve_unknown_preset():
    with pytest.raises(ValueError, match="Unknown preset"):
        resolve_preset("nonexistent")


def test_resolve_mixture_tsv_colon_wins():
    result = resolve_mixture(
        genotype_col="GII.4:0.6,GI.1:0.4",
        cli_mixture="GII.17:1.0",
        cli_preset=None,
    )
    assert result.lineages[0].genotype == "GII.4"


def test_resolve_mixture_cli_mixture_fallback():
    result = resolve_mixture(
        genotype_col="GII.4",
        cli_mixture="GII.17:0.7,GII.4:0.3",
        cli_preset=None,
    )
    assert len(result.lineages) == 2
    assert result.lineages[0].genotype == "GII.17"


def test_resolve_mixture_cli_preset_fallback():
    result = resolve_mixture(
        genotype_col="GII.4",
        cli_mixture=None,
        cli_preset="us-2024",
    )
    assert len(result.lineages) == 6


def test_resolve_mixture_bare_genotype():
    result = resolve_mixture(
        genotype_col="GII.4",
        cli_mixture=None,
        cli_preset=None,
    )
    assert result == MixtureSpec(lineages=[LineageSpec("GII.4", 1.0)])


def test_generate_splits_proportionally():
    config = _make_config()
    genomes = _make_genome_records()
    rng = np.random.default_rng(42)
    mixture = MixtureSpec(
        lineages=[
            LineageSpec("GII.4", 0.7),
            LineageSpec("GII.17", 0.3),
        ]
    )

    result = generate_mixture_fragments(mixture, 1000, genomes, config, rng)
    assert "GII.4" in result
    assert "GII.17" in result
    n_gii4 = len(result["GII.4"][2])
    n_gii17 = len(result["GII.17"][2])
    assert n_gii4 + n_gii17 == 1000
    assert n_gii4 == 700
    assert n_gii17 == 300


def test_generate_single_lineage():
    config = _make_config()
    genomes = _make_genome_records()
    rng = np.random.default_rng(42)
    mixture = MixtureSpec(lineages=[LineageSpec("GII.4", 1.0)])

    result = generate_mixture_fragments(mixture, 500, genomes, config, rng)
    assert len(result) == 1
    assert len(result["GII.4"][2]) == 500


def test_generate_tiny_abundance_skipped():
    config = _make_config()
    genomes = _make_genome_records()
    rng = np.random.default_rng(42)
    mixture = MixtureSpec(
        lineages=[
            LineageSpec("GII.4", 0.999),
            LineageSpec("GII.17", 0.001),
        ]
    )

    result = generate_mixture_fragments(mixture, 100, genomes, config, rng)
    assert "GII.4" in result
    # GII.17 gets 0 fragments, should be skipped
    if "GII.17" in result:
        assert len(result["GII.17"][2]) == 0


def test_generate_returns_vp1_coords():
    config = _make_config()
    genomes = _make_genome_records()
    rng = np.random.default_rng(42)
    mixture = MixtureSpec(lineages=[LineageSpec("GII.4", 1.0)])

    result = generate_mixture_fragments(mixture, 100, genomes, config, rng)
    vp1_start, vp1_end, frags = result["GII.4"]
    assert isinstance(vp1_start, int)
    assert isinstance(vp1_end, int)
    assert vp1_start < vp1_end


def test_generate_fragment_source_is_genotype():
    config = _make_config()
    genomes = _make_genome_records()
    rng = np.random.default_rng(42)
    mixture = MixtureSpec(
        lineages=[
            LineageSpec("GII.4", 0.5),
            LineageSpec("GII.17", 0.5),
        ]
    )

    result = generate_mixture_fragments(mixture, 100, genomes, config, rng)
    for frag in result["GII.4"][2]:
        assert frag.source == "GII.4"
    for frag in result["GII.17"][2]:
        assert frag.source == "GII.17"
