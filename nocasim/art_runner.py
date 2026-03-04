import subprocess
from pathlib import Path

from nocasim.config import SimConfig
from nocasim.fragment import Fragment


def verify_art_modern(bin_path: Path) -> str:
    try:
        result = subprocess.run(
            [str(bin_path), "--version"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        output = result.stdout.strip() or result.stderr.strip()
        if not output:
            raise RuntimeError(f"art_modern at {bin_path} produced no version output")
        return output
    except FileNotFoundError:
        raise FileNotFoundError(f"art_modern not found at {bin_path}")


def write_pbsim3_tsv(
    fragments: list[tuple[Fragment, int]],
    path: Path,
) -> None:
    with open(path, "w") as f:
        f.write("#ID\tCOV_POS\tCOV_NEG\tSEQ\n")
        for frag, copy_count in fragments:
            if frag.strand == "+":
                cov_pos, cov_neg = float(copy_count), 0.0
            else:
                cov_pos, cov_neg = 0.0, float(copy_count)
            f.write(f"{frag.id}\t{cov_pos}\t{cov_neg}\t{frag.sequence}\n")


def run_art_modern(
    tsv_path: Path,
    sample_id: str,
    outdir: Path,
    config: SimConfig,
) -> Path:
    output_prefix = outdir / sample_id
    cmd = [
        str(config.art_modern_bin),
        "--mode",
        "template",
        "--lc",
        "pe",
        "--i-file",
        str(tsv_path),
        "--i-type",
        "pbsim3_transcripts",
        "--i-parser",
        "stream",
        "--o-fastq",
        str(output_prefix),
        "--builtin_qual_file",
        "HiSeq2500_150bp",
        "--read_len",
        str(config.read_len),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"art_modern exited with code {result.returncode}:\n{result.stderr}"
        )

    return output_prefix


def deinterleave_fastq(interleaved: Path, r1_path: Path, r2_path: Path) -> None:
    with open(interleaved) as fin, open(r1_path, "w") as f1, open(r2_path, "w") as f2:
        while True:
            lines_r1 = [fin.readline() for _ in range(4)]
            if not lines_r1[0]:
                break
            lines_r2 = [fin.readline() for _ in range(4)]
            if not lines_r2[0]:
                break
            f1.writelines(lines_r1)
            f2.writelines(lines_r2)
