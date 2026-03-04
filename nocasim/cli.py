import csv
import tempfile
from pathlib import Path

import click
import numpy as np

from nocasim.art_runner import run_art_modern, verify_art_modern, write_pbsim3_tsv
from nocasim.background import generate_background_fragments
from nocasim.capture import run_capture
from nocasim.config import SimConfig
from nocasim.ct_model import ct_to_viral_fraction
from nocasim.duplicate import assign_copy_counts
from nocasim.fragment import sample_fragments
from nocasim.genome import load_fasta
from nocasim.output import organise_outputs
from nocasim.truth import compute_vp1_coverage_from_fragments, write_summary


@click.group()
def main():
    """nocasim - Norovirus Capture Sequencing Simulator."""


def _simulate_sample(
    sample_id: str,
    genotype: str,
    ct: float,
    genome_records: dict,
    config: SimConfig,
    rng: np.random.Generator,
    tmpdir: Path,
    keep_intermediates: bool = False,
    total_library: int = 500_000,
    human_ref: dict | None = None,
    microbiome_ref: dict | None = None,
    wastewater_ref: dict | None = None,
) -> dict:
    if genotype not in genome_records:
        raise click.ClickException(
            f"Genotype '{genotype}' not found in references. "
            f"Available: {list(genome_records.keys())}"
        )
    genome = genome_records[genotype]

    vp1_start = config.vp1_start if genome.length > config.vp1_end else 0
    vp1_end = config.vp1_end if genome.length > config.vp1_end else genome.length
    vp1_length = vp1_end - vp1_start

    vf = ct_to_viral_fraction(ct)
    n_viral = max(int(total_library * vf), 100)
    n_bg = total_library - n_viral
    click.echo(
        f"  [{sample_id}] Ct={ct}, viral_fraction={vf:.4f}, "
        f"n_viral={n_viral}, n_bg={n_bg}"
    )

    viral_frags = sample_fragments(
        genome,
        n_viral,
        config,
        "viral",
        vp1_start,
        vp1_end,
        rng,
    )
    bg_frags = []
    for batch in generate_background_fragments(
        n_bg, config, human_ref, microbiome_ref, rng, wastewater_ref
    ):
        bg_frags.extend(batch)

    on_target, off_target = run_capture(viral_frags, bg_frags, config, rng)
    click.echo(
        f"  [{sample_id}] Captured: {len(on_target)} on-target, "
        f"{len(off_target)} off-target"
    )

    on_copies = assign_copy_counts(on_target, config.dup_rate, rng)
    off_copies = assign_copy_counts(off_target, config.dup_rate, rng)

    all_copies = on_copies + off_copies
    if not all_copies:
        click.echo(f"  [{sample_id}] WARNING: no fragments after capture")
        return {
            "sample_id": sample_id,
            "genotype": genotype,
            "ct_value": ct,
            "vp1_mean_depth": 0.0,
            "vp1_completeness_20x": 0.0,
            "completeness_call": "incomplete",
        }

    tsv_path = tmpdir / f"{sample_id}_fragments.tsv"
    write_pbsim3_tsv(all_copies, tsv_path)

    art_prefix = run_art_modern(tsv_path, sample_id, tmpdir, config)

    coverage = compute_vp1_coverage_from_fragments(on_copies, vp1_start, vp1_end)

    total_on = sum(c for _, c in on_copies)
    total_off = sum(c for _, c in off_copies)
    total_reads = total_on + total_off
    achieved_on_target = total_on / total_reads if total_reads > 0 else 0.0
    achieved_dup = (
        1.0 - len(all_copies) / sum(c for _, c in all_copies) if all_copies else 0.0
    )
    frag_lengths = [f.end - f.start for f, _ in all_copies]
    mean_frag = sum(frag_lengths) / len(frag_lengths) if frag_lengths else 0.0

    stats = {
        "sample_id": sample_id,
        "genotype": genotype,
        "ct_value": ct,
        "seed": config.seed,
        "viral_fraction_input": round(vf, 4),
        "vp1_length_bp": coverage.get("vp1_length_bp", vp1_length),
        "total_read_pairs": total_reads,
        "on_target_read_pairs": total_on,
        "off_target_read_pairs": total_off,
        "on_target_rate_achieved": round(achieved_on_target, 3),
        "vp1_mean_depth": coverage.get("vp1_mean_depth", 0.0),
        "vp1_completeness_20x": coverage.get("vp1_completeness_20x", 0.0),
        "completeness_call": coverage.get("completeness_call", "incomplete"),
        "duplicate_rate_achieved": round(achieved_dup, 3),
        "fragment_mean_achieved": round(mean_frag, 1),
    }

    organise_outputs(sample_id, art_prefix, stats, config.outdir)
    click.echo(
        f"  [{sample_id}] depth={stats['vp1_mean_depth']}, "
        f"call={stats['completeness_call']}"
    )

    if not keep_intermediates and tsv_path.exists():
        tsv_path.unlink()

    return stats


@main.command()
@click.option(
    "--sample-sheet", required=True, type=click.Path(exists=True, path_type=Path)
)
@click.option(
    "--references", required=True, type=click.Path(exists=True, path_type=Path)
)
@click.option("--probes", type=click.Path(exists=True, path_type=Path), default=None)
@click.option("--human-bg", type=click.Path(exists=True, path_type=Path), default=None)
@click.option(
    "--microbiome-bg", type=click.Path(exists=True, path_type=Path), default=None
)
@click.option(
    "--wastewater-bg", type=click.Path(exists=True, path_type=Path), default=None
)
@click.option(
    "--sample-type",
    default="stool",
    type=click.Choice(["stool", "wastewater"]),
    help="Sample type controls background composition mix",
)
@click.option("--art-modern", required=True, type=click.Path(path_type=Path))
@click.option("--read-len", default=150, type=int)
@click.option("--dup-rate", default=0.40, type=float)
@click.option("--off-target", default=0.592, type=float)
@click.option("--outdir", required=True, type=click.Path(path_type=Path))
@click.option("--seed", default=42, type=int)
@click.option("--keep-intermediates", is_flag=True, default=False)
@click.option(
    "--total-fragments",
    default=500_000,
    type=int,
    help="Total pre-capture library size per sample",
)
def simulate(
    sample_sheet,
    references,
    probes,
    human_bg,
    microbiome_bg,
    wastewater_bg,
    sample_type,
    art_modern,
    read_len,
    dup_rate,
    off_target,
    outdir,
    seed,
    keep_intermediates,
    total_fragments,
):
    """Simulate a batch of samples from a sample sheet."""
    config = SimConfig(
        references_dir=references,
        art_modern_bin=art_modern,
        outdir=outdir,
        probes_fasta=probes,
        human_bg_fasta=human_bg,
        microbiome_bg_fasta=microbiome_bg,
        wastewater_bg_fasta=wastewater_bg,
        sample_type=sample_type,
        read_len=read_len,
        dup_rate=dup_rate,
        off_target_rate=off_target,
        seed=seed,
    )

    version = verify_art_modern(config.art_modern_bin)
    click.echo(f"art_modern version: {version}")

    genome_records = {}
    for fasta_path in references.glob("*.fasta"):
        records = load_fasta(fasta_path)
        for name, rec in records.items():
            genome_records[fasta_path.stem] = rec

    human_ref = load_fasta(human_bg) if human_bg else None
    microbiome_ref = load_fasta(microbiome_bg) if microbiome_bg else None
    wastewater_ref = load_fasta(wastewater_bg) if wastewater_bg else None
    if human_ref:
        click.echo(f"Loaded human background: {len(human_ref)} sequences")
    if microbiome_ref:
        click.echo(f"Loaded microbiome background: {len(microbiome_ref)} sequences")
    if wastewater_ref:
        click.echo(f"Loaded wastewater background: {len(wastewater_ref)} sequences")

    samples = []
    with open(sample_sheet) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            samples.append(
                {
                    "sample_id": row["sample_id"],
                    "genotype": row["genotype"],
                    "ct_value": float(row["ct_value"]),
                }
            )

    click.echo(f"Loaded {len(samples)} samples, {len(genome_records)} references")
    outdir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(seed)
    all_stats = []

    with tempfile.TemporaryDirectory() as tmpdir:
        for sample in samples:
            stats = _simulate_sample(
                sample_id=sample["sample_id"],
                genotype=sample["genotype"],
                ct=sample["ct_value"],
                genome_records=genome_records,
                config=config,
                rng=rng,
                tmpdir=Path(tmpdir),
                keep_intermediates=keep_intermediates,
                total_library=total_fragments,
                human_ref=human_ref,
                microbiome_ref=microbiome_ref,
                wastewater_ref=wastewater_ref,
            )
            all_stats.append(stats)

    summary_path = outdir / "summary.tsv"
    write_summary(all_stats, summary_path)
    click.echo(f"Summary written to {summary_path}")


@main.command()
@click.option(
    "--reference", required=True, type=click.Path(exists=True, path_type=Path)
)
@click.option("--ct", required=True, type=float)
@click.option("--outdir", required=True, type=click.Path(path_type=Path))
@click.option("--art-modern", required=True, type=click.Path(path_type=Path))
@click.option("--read-len", default=150, type=int)
@click.option("--dup-rate", default=0.40, type=float)
@click.option("--off-target", default=0.592, type=float)
@click.option("--seed", default=42, type=int)
@click.option(
    "--n-fragments",
    default=10_000,
    type=int,
    help="Number of viral fragments to generate",
)
def single(
    reference, ct, outdir, art_modern, read_len, dup_rate, off_target, seed, n_fragments
):
    """Simulate a single sample for testing."""
    config = SimConfig(
        references_dir=reference.parent,
        art_modern_bin=art_modern,
        outdir=outdir,
        read_len=read_len,
        dup_rate=dup_rate,
        off_target_rate=off_target,
        seed=seed,
    )

    version = verify_art_modern(config.art_modern_bin)
    click.echo(f"art_modern version: {version}")

    rng = np.random.default_rng(seed)

    records = load_fasta(reference)
    genome = next(iter(records.values()))
    click.echo(f"Loaded reference: {genome.name} ({genome.length} bp)")

    vp1_start = config.vp1_start if genome.length > config.vp1_end else 0
    vp1_end = config.vp1_end if genome.length > config.vp1_end else genome.length

    viral_frags = sample_fragments(
        genome,
        n_fragments,
        config,
        "viral",
        vp1_start,
        vp1_end,
        rng,
    )
    click.echo(f"Sampled {len(viral_frags)} viral fragments")

    on_target, off_target_frags = run_capture(viral_frags, [], config, rng)
    click.echo(f"Captured: {len(on_target)} on-target")

    on_copies = assign_copy_counts(on_target, config.dup_rate, rng)

    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        tsv_path = Path(tmpdir) / "fragments.tsv"
        write_pbsim3_tsv(on_copies, tsv_path)

        sample_id = reference.stem
        art_prefix = run_art_modern(tsv_path, sample_id, Path(tmpdir), config)

        coverage = compute_vp1_coverage_from_fragments(
            on_copies,
            vp1_start,
            vp1_end,
        )

        manifest = {
            "sample_id": sample_id,
            "reference": str(reference),
            "ct_value": ct,
            "seed": seed,
            "n_fragments": n_fragments,
            "genome_length": genome.length,
            **coverage,
        }
        organise_outputs(sample_id, art_prefix, manifest, outdir)

    click.echo(f"Output written to {outdir / sample_id}/")


@main.command("download-probes")
@click.option("--outdir", required=True, type=click.Path(path_type=Path))
def download_probes(outdir):
    """Download probe sequences from PMC12216758 Supplementary File 2."""
    import urllib.request

    outdir.mkdir(parents=True, exist_ok=True)
    url = (
        "https://static-content.springer.com/esm/"
        "art%3A10.1038%2Fs41598-025-03398-6/MediaObjects/"
        "41598_2025_3398_MOESM2_ESM.txt"
    )
    dest = outdir / "hunov_probes.txt"

    click.echo(f"Downloading probes from {url}")
    try:
        urllib.request.urlretrieve(url, dest)
        click.echo(f"Saved to {dest}")
    except Exception as e:
        raise click.ClickException(f"Download failed: {e}")
