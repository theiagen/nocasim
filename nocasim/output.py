import gzip
import json
import shutil
from pathlib import Path

from nocasim.art_runner import deinterleave_fastq


def compress_fastq(src: Path, dst: Path) -> None:
    with open(src, "rb") as f_in, gzip.open(dst, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def organise_outputs(
    sample_id: str,
    art_output_prefix: Path,
    manifest: dict,
    outdir: Path,
) -> None:
    sample_dir = outdir / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    interleaved = Path(str(art_output_prefix))
    if interleaved.exists() and interleaved.stat().st_size > 0:
        r1_tmp = interleaved.parent / f"{sample_id}_R1.fq"
        r2_tmp = interleaved.parent / f"{sample_id}_R2.fq"
        deinterleave_fastq(interleaved, r1_tmp, r2_tmp)
        compress_fastq(r1_tmp, sample_dir / f"{sample_id}_R1.fastq.gz")
        compress_fastq(r2_tmp, sample_dir / f"{sample_id}_R2.fastq.gz")
    else:
        for suffix in ("_R1.fq", "_1.fq"):
            r1_raw = Path(str(art_output_prefix) + suffix)
            if r1_raw.exists():
                compress_fastq(r1_raw, sample_dir / f"{sample_id}_R1.fastq.gz")
                break
        for suffix in ("_R2.fq", "_2.fq"):
            r2_raw = Path(str(art_output_prefix) + suffix)
            if r2_raw.exists():
                compress_fastq(r2_raw, sample_dir / f"{sample_id}_R2.fastq.gz")
                break

    manifest_path = sample_dir / f"{sample_id}_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
