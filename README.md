# nocasim

Norovirus Capture Sequencing Simulator. Generates realistic paired-end FASTQ
from simulated hybrid capture enrichment of Norovirus VP1, calibrated from
empirical data in Bhamidipati et al. (2025) *Sci Rep* 15:20526.

## Background

Human Norovirus genotype is defined by the VP1 region (ORF2, ~1,700 bp).
Hybrid capture enrichment followed by paired-end sequencing on NovaSeq is the
current gold standard for lineage surveillance. This simulator produces
controlled FASTQ datasets — with known genotypes, Ct values, and ground truth
coverage — so that full analysis pipelines can be benchmarked before deployment.

**Empirical calibration source:**
> Bhamidipati et al. (2025) "Complete genomic characterization of RSV and HuNoV
> using probe-based capture enrichment." *Sci Rep* 15:20526 (PMC12216758)

## How It Works

The simulation pipeline has six stages that mirror the wet-lab workflow:

**1. Ct-to-viral-fraction model.** Each sample's Ct value is converted to a
viral fraction using `vf = 0.03 * 2^((25 - ct) / 3.32)`. The exponential
relationship follows standard qPCR theory — each 3.32 Ct difference
corresponds to a 10-fold change in template concentration (`2^3.32 ≈ 10`).
The calibration parameters (base fraction 0.03 at Ct 25) were tuned to
produce coverage outputs consistent with the empirical results in Table 1
of Bhamidipati et al. (2025); they are not directly reported in the paper.
A Ct of 20 gives ~30% viral reads; a Ct of 35 gives ~0.03%. The total
pre-capture library (default 500,000 fragments) is split into `n_viral` and
`n_background` based on this fraction.

**2. Fragment generation.** Viral fragments are sampled from the reference
genome using a truncated normal length distribution (mean=380 bp, sd=80 bp,
range 200-600 bp) matching heat fragmentation at 94 C for 10 minutes. Each
fragment records its genomic coordinates, strand, GC content, and whether it
overlaps VP1. Background fragments come from human DNA (GC ~0.41), gut
microbiome (GC ~0.52), or wastewater environmental sources (GC ~0.48), mixed
according to sample type (stool: 80/20 human/microbiome; wastewater:
10/30/60 human/microbiome/environmental).

**3. Capture enrichment.** Hybrid capture is modelled per-fragment. Viral
fragments that overlap VP1 are retained with probability
`P = exp(-6.0 * (gc - 0.47)^2)`, a Gaussian centered on the optimal probe GC
of 0.47. Fragments outside VP1 and all background fragments pass through at
the configured off-target rate (default 59.2%, from 1 - 0.408 on-target rate
in Table 1 of Bhamidipati et al.). This produces on-target and off-target
fragment pools.

**4. PCR duplicate model.** Each captured fragment is assigned a copy count
drawn from a geometric distribution with `p = 1 - dup_rate`. At the default
40% duplicate rate, most fragments appear once, but some are amplified 2-5x,
matching the 12-20 PCR cycles used in post-capture amplification.

**5. Read simulation.** The duplicated fragment pool is written as a PBSIM3
TSV (sequence, start, end per fragment) and passed to
[art_modern](https://github.com/YU-Zhejian/art_modern), which generates
paired-end 2x150 bp reads with Illumina NovaSeq-calibrated base quality
scores and error profiles. The output is deinterleaved into R1/R2 FASTQ files
and gzip-compressed.

**6. Truth and QC.** Ground truth VP1 coverage is computed directly from
fragment coordinates: mean depth and the fraction of VP1 bases at >= 20x.
Each sample receives a completeness call — `complete` (>= 90% of VP1 at
>= 20x), `low_coverage` (>= 90% breadth but < 20x mean), or `incomplete`
(< 90% breadth) — following the criteria in the paper.

## Requirements

- Python >= 3.11
- [art_modern](https://github.com/YU-Zhejian/art_modern) (external binary, see below)
- Python dependencies: `click`, `biopython`, `numpy`, `scipy`

## Installation

```bash
git clone <repo>
cd nocasim
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[dev]"
```

### Installing art_modern (Ubuntu 24.04)

Download the `.deb` from the
[art_modern releases page](https://github.com/YU-Zhejian/art_modern/releases),
then install it:

```bash
sudo dpkg -i art-modern_*.deb
sudo apt-get install -f -y
```

Verify the install:

```bash
art_modern --version
```

### Docker

```bash
docker build -t nocasim .
```

Run with mounted output directory:

```bash
docker run --rm -v $(pwd)/results:/output nocasim single \
  --reference /opt/nocasim/data/references/GII.4.fasta \
  --ct 28.0 --outdir /output --art-modern art_modern
```

For batch mode, mount your sample sheet and references:

```bash
docker run --rm \
  -v $(pwd)/samples.tsv:/data/samples.tsv \
  -v $(pwd)/results:/output \
  nocasim simulate \
    --sample-sheet /data/samples.tsv \
    --references /opt/nocasim/data/references/ \
    --art-modern art_modern \
    --outdir /output
```

## Quick Start

### Single sample

```bash
nocasim single \
  --reference data/references/GII.4.fasta \
  --ct 28.0 \
  --outdir results/ \
  --art-modern art_modern
```

### Batch mode

```bash
nocasim simulate \
  --sample-sheet samples.tsv \
  --references data/references/ \
  --art-modern art_modern \
  --outdir results/
```

### Wastewater sample type

Use `--sample-type wastewater` to simulate wastewater surveillance samples.
The background composition shifts from stool (80% human / 20% microbiome) to
wastewater (10% human / 30% microbiome / 60% environmental):

```bash
nocasim simulate \
  --sample-sheet samples.tsv \
  --references data/references/ \
  --art-modern art_modern \
  --sample-type wastewater \
  --outdir results/
```

Optionally provide a wastewater background FASTA for realistic environmental
sequences instead of synthetic ones:

```bash
nocasim simulate \
  --sample-sheet samples.tsv \
  --references data/references/ \
  --art-modern art_modern \
  --sample-type wastewater \
  --wastewater-bg data/background/wastewater_metagenome.fasta \
  --outdir results/
```

### Fetch probe sequences from the paper

```bash
nocasim download-probes --outdir data/probes/
```

This downloads Supplementary File 2 from PMC12216758 and saves it as
`data/probes/hunov_probes.txt` (BED-like probe interval file).

## Sample Sheet Format

The sample sheet is a TSV file with three columns:

```
sample_id    genotype    ct_value
sample_001   GII.4       24.5
sample_002   GII.17      30.0
sample_003   GI.1        31.7
sample_004   GII.2       26.0
sample_005   GII.4       35.2
```

The `genotype` value must match the filename stem of a FASTA file in the
`--references` directory. For example, `GII.4` maps to
`data/references/GII.4.fasta`.

An example sample sheet is included at `data/samples.tsv`.

## Reference Files

The repository includes 35 reference sequences (9 GI + 26 GII genotypes)
from the [CDC Calicivirus Typing Tool](https://calicivirustypingtool.cdc.gov/becerance.cgi).
One FASTA file per genotype in `data/references/`, where the filename stem
matches the genotype string in the sample sheet:

```
data/references/
├── GI.1.fasta    (EU085529, 3081 bp)
├── GI.2.fasta    (AF435807, 2354 bp)
├── ...
├── GII.4.fasta   (X76716,   3881 bp)
├── GII.17.fasta  (KJ156329, 3723 bp)
├── ...
└── GII.27.fasta  (MK733205, 7308 bp)
```

References range from partial VP1 (~1,000 bp) to full genome (~7,700 bp).
Full-genome references use the configured `vp1_start`/`vp1_end` coordinates
to extract the target region (defaults: 5,100-6,800 bp). Shorter sequences
are treated as the target region in their entirety.

## Output Structure

```
results/
├── summary.tsv
└── sample_001/
    ├── sample_001_R1.fastq.gz
    ├── sample_001_R2.fastq.gz
    └── sample_001_manifest.json
```

`summary.tsv` contains one row per sample with columns: `sample_id`,
`genotype`, `ct_value`, `vp1_depth`, `vp1_completeness_20x`,
`completeness_call`.

`sample_001_manifest.json` contains per-sample ground truth metrics including
achieved on-target rate, duplicate rate, mean VP1 depth, and completeness call.

Completeness calls follow the paper's criteria:
- `complete`: >= 20x coverage across >= 90% of VP1
- `low_coverage`: >= 90% of VP1 breadth covered but mean depth < 20x
- `incomplete`: < 90% of VP1 covered

## Parameters

| Parameter | Default | Source |
|---|---|---|
| `--read-len` | 150 | NovaSeq 2x150 bp |
| `--dup-rate` | 0.40 | 12-20 post-capture PCR cycles |
| `--off-target` | 0.592 | 1 - 0.408 on-target rate (Table 1) |
| `--total-fragments` | 500000 | Pre-capture library size |
| `--sample-type` | stool | Background composition model (`stool` or `wastewater`) |
| `--seed` | 42 | RNG seed for reproducibility |

### Background composition by sample type

| Source | Stool | Wastewater |
|---|---|---|
| Human | 80% | 10% |
| Microbiome | 20% | 30% |
| Environmental | 0% | 60% |

Background FASTA files (`--human-bg`, `--microbiome-bg`, `--wastewater-bg`) are
optional. When omitted, synthetic sequences with realistic GC content profiles
are generated as placeholders.

Fragment length distribution is a truncated normal (mean=380 bp, sd=80 bp,
range 200-600 bp) matching heat fragmentation at 94 degrees C for 10 minutes
as described in the paper. This differs from sonication (log-normal) — do not
substitute a log-normal distribution.

## Running Tests

```bash
pytest tests/ -v
```

## What This Simulator Does Not Model

- ORF1 / polymerase region (recombination makes it irrelevant for genotyping)
- Intra-host variant frequencies or quasispecies dynamics
- RNA extraction efficiency
- Reverse transcription efficiency
