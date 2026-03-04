from pathlib import Path
from unittest.mock import patch, MagicMock

from nocasim.art_runner import verify_art_modern, write_pbsim3_tsv
from nocasim.fragment import Fragment


def test_write_pbsim3_tsv_format(tmp_path):
    frags = [
        (
            Fragment(
                id="frag_0001_viral_0_350_+",
                sequence="ACGT" * 50,
                start=0,
                end=350,
                strand="+",
                gc=0.5,
                source="viral",
                overlaps_vp1=True,
            ),
            2,
        ),
        (
            Fragment(
                id="frag_0002_viral_200_580_-",
                sequence="TGCA" * 50,
                start=200,
                end=580,
                strand="-",
                gc=0.5,
                source="viral",
                overlaps_vp1=True,
            ),
            1,
        ),
    ]

    tsv_path = tmp_path / "test.tsv"
    write_pbsim3_tsv(frags, tsv_path)

    lines = tsv_path.read_text().strip().split("\n")
    assert lines[0].startswith("#")
    assert len(lines) == 3

    fields = lines[1].split("\t")
    assert fields[0] == "frag_0001_viral_0_350_+"
    assert float(fields[1]) == 2.0
    assert float(fields[2]) == 0.0

    fields = lines[2].split("\t")
    assert fields[0] == "frag_0002_viral_200_580_-"
    assert float(fields[1]) == 0.0
    assert float(fields[2]) == 1.0


def test_verify_art_modern_not_found():
    import pytest

    with pytest.raises(FileNotFoundError):
        verify_art_modern(Path("/nonexistent/art_modern"))


@patch("nocasim.art_runner.subprocess.run")
def test_verify_art_modern_success(mock_run):
    mock_run.return_value = MagicMock(
        stdout="art_modern v1.0.0", stderr="", returncode=0
    )
    version = verify_art_modern(Path("/usr/bin/art_modern"))
    assert "art_modern" in version
