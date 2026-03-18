from pathlib import Path

from click.testing import CliRunner

from nocasim.cli import main


FIXTURES = Path(__file__).parent / "fixtures"


def test_simulate_with_preset_validates_genotypes(tmp_path):
    """Preset references genotypes not in test fixtures -> should error."""
    sample_sheet = tmp_path / "samples.tsv"
    sample_sheet.write_text("sample_id\tgenotype\tct_value\nsample_001\tGII.4\t25.0\n")

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "simulate",
            "--sample-sheet",
            str(sample_sheet),
            "--references",
            str(FIXTURES),
            "--preset",
            "us-2024",
            "--art-modern",
            "/usr/bin/false",
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    # us-2024 preset requires GII.17, GII.4, GII.2, etc. — only test_reference.fasta exists
    assert result.exit_code != 0
    assert (
        "not found" in result.output.lower()
        or "not found" in str(result.exception).lower()
    )


def test_simulate_mixture_mutual_exclusion(tmp_path):
    """Passing both --mixture and --preset should error."""
    sample_sheet = tmp_path / "samples.tsv"
    sample_sheet.write_text("sample_id\tgenotype\tct_value\nsample_001\tGII.4\t25.0\n")

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "simulate",
            "--sample-sheet",
            str(sample_sheet),
            "--references",
            str(FIXTURES),
            "--mixture",
            "GII.4:1.0",
            "--preset",
            "outbreak",
            "--art-modern",
            "/usr/bin/false",
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code != 0
