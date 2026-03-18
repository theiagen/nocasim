import pytest

from nocasim.mixture import LineageSpec, MixtureSpec, parse_mixture


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
