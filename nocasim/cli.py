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
from nocasim.mixture import (
    PRESETS,
    MixtureSpec,
    generate_mixture_fragments,
    resolve_mixture,
)
from nocasim.output import organise_outputs
from nocasim.truth import (
    compute_mixture_coverage,
    compute_vp1_coverage_from_fragments,
    write_summary,
)


@click.group()
def main():
    """nocasim - Norovirus Capture Sequencing Simulator."""


def _simulate_sample(
    sample_id: str,
    mixture: MixtureSpec,
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
    genotype_str = mixture.to_spec_string()

    vf = ct_to_viral_fraction(ct)
    n_viral = max(int(total_library * vf), 100)
    n_bg = total_library - n_viral
    click.echo(
        f"  [{sample_id}] Ct={ct}, viral_fraction={vf:.4f}, "
        f"n_viral={n_viral}, n_bg={n_bg}"
    )

    lineage_data = generate_mixture_fragments(
        mixture, n_viral, genome_records, config, rng
    )

    viral_frags = []
    for _genotype, (_vs, _ve, frags) in lineage_data.items():
        viral_frags.extend(frags)

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
            "genotype": genotype_str,
            "ct_value": ct,
            "vp1_mean_depth": 0.0,
            "vp1_completeness_20x": 0.0,
            "completeness_call": "incomplete",
            "lineage_detail": "",
        }

    tsv_path = tmpdir / f"{sample_id}_fragments.tsv"
    write_pbsim3_tsv(all_copies, tsv_path)

    art_prefix = run_art_modern(tsv_path, sample_id, tmpdir, config)

    abundances = {lineage.genotype: lineage.abundance for lineage in mixture.lineages}
    lineage_on_copies = {
        g: (vs, ve, [(f, c) for f, c in on_copies if f.source == g])
        for g, (vs, ve, _) in lineage_data.items()
    }
    mix_coverage = compute_mixture_coverage(lineage_on_copies, abundances)
    coverage = mix_coverage["aggregate"]

    if len(lineage_data) > 1:
        detail_parts = []
        for g, lstats in mix_coverage["per_lineage"].items():
            detail_parts.append(
                f"{g}:{lstats['vp1_mean_depth']}x/{lstats['vp1_completeness_20x']}"
            )
        lineage_detail = ";".join(detail_parts)
    else:
        lineage_detail = ""

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
        "genotype": genotype_str,
        "ct_value": ct,
        "seed": config.seed,
        "viral_fraction_input": round(vf, 4),
        "vp1_length_bp": coverage.get("vp1_length_bp", 0),
        "total_read_pairs": total_reads,
        "on_target_read_pairs": total_on,
        "off_target_read_pairs": total_off,
        "on_target_rate_achieved": round(achieved_on_target, 3),
        "vp1_mean_depth": coverage.get("vp1_mean_depth", 0.0),
        "vp1_completeness_20x": coverage.get("vp1_completeness_20x", 0.0),
        "completeness_call": coverage.get("completeness_call", "incomplete"),
        "duplicate_rate_achieved": round(achieved_dup, 3),
        "fragment_mean_achieved": round(mean_frag, 1),
        "lineage_detail": lineage_detail,
    }

    manifest = dict(stats)
    if len(mixture.lineages) > 1:
        manifest["mixture"] = [
            {"genotype": lineage.genotype, "abundance": lineage.abundance}
            for lineage in mixture.lineages
        ]
        manifest["per_lineage"] = {
            g: {
                "abundance": abundances[g],
                "n_fragments": len(lineage_data[g][2]),
                **lstats,
            }
            for g, lstats in mix_coverage["per_lineage"].items()
        }
        manifest["aggregate"] = coverage

    organise_outputs(sample_id, art_prefix, manifest, config.outdir)
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
@click.option(
    "--mixture",
    "mixture_spec",
    default=None,
    type=str,
    help='Mixture spec for all samples, e.g. "GII.17:0.75,GII.4:0.25"',
)
@click.option(
    "--preset",
    default=None,
    type=click.Choice(list(PRESETS.keys())),
    help="Built-in mixture preset name",
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
    mixture_spec,
    preset,
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

    if mixture_spec and preset:
        raise click.ClickException("--mixture and --preset are mutually exclusive")

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

    for sample in samples:
        mix = resolve_mixture(sample["genotype"], mixture_spec, preset)
        for lineage in mix.lineages:
            if lineage.genotype not in genome_records:
                raise click.ClickException(
                    f"Genotype '{lineage.genotype}' (sample {sample['sample_id']}) "
                    f"not found in references. "
                    f"Available: {sorted(genome_records.keys())}"
                )

    rng = np.random.default_rng(seed)
    all_stats = []

    with tempfile.TemporaryDirectory() as tmpdir:
        for sample in samples:
            mix = resolve_mixture(sample["genotype"], mixture_spec, preset)
            stats = _simulate_sample(
                sample_id=sample["sample_id"],
                mixture=mix,
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
