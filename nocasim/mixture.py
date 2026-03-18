from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class LineageSpec:
    genotype: str
    abundance: float


@dataclass(frozen=True)
class MixtureSpec:
    lineages: list[LineageSpec]

    def to_spec_string(self) -> str:
        if len(self.lineages) == 1 and self.lineages[0].abundance == 1.0:
            return self.lineages[0].genotype
        return ",".join(f"{spec.genotype}:{spec.abundance}" for spec in self.lineages)


_ABUNDANCE_SUM_TOLERANCE = 0.01


def parse_mixture(spec: str) -> MixtureSpec:
    spec = spec.strip()
    if not spec:
        raise ValueError("Empty mixture spec")

    if ":" not in spec:
        return MixtureSpec(lineages=[LineageSpec(spec, 1.0)])

    parts = spec.split(",")
    lineages: list[LineageSpec] = []
    seen: set[str] = set()

    for part in parts:
        genotype, _, abundance_str = part.strip().partition(":")
        genotype = genotype.strip()
        abundance = float(abundance_str.strip())

        if abundance <= 0.0 or abundance > 1.0:
            raise ValueError(
                f"Abundance for '{genotype}' must be > 0 and <= 1, got {abundance}"
            )
        if genotype in seen:
            raise ValueError(f"Duplicate genotype: '{genotype}'")

        seen.add(genotype)
        lineages.append(LineageSpec(genotype, abundance))

    total = sum(spec.abundance for spec in lineages)
    if abs(total - 1.0) > _ABUNDANCE_SUM_TOLERANCE:
        raise ValueError(
            f"Abundances sum to {total:.4f}, expected 1.0 "
            f"(tolerance {_ABUNDANCE_SUM_TOLERANCE})"
        )

    if abs(total - 1.0) > 1e-9:
        lineages = [
            LineageSpec(spec.genotype, spec.abundance / total) for spec in lineages
        ]

    return MixtureSpec(lineages=lineages)


# CaliciNet 2024-25 (PMC12205451), minor genotypes estimated
PRESETS: dict[str, MixtureSpec] = {
    "us-2024": MixtureSpec(
        lineages=[
            LineageSpec("GII.17", 0.75),
            LineageSpec("GII.4", 0.11),
            LineageSpec("GII.2", 0.05),
            LineageSpec("GI.1", 0.04),
            LineageSpec("GII.6", 0.03),
            LineageSpec("GI.3", 0.02),
        ]
    ),
    "diverse": MixtureSpec(
        lineages=[
            LineageSpec("GII.17", 0.25),
            LineageSpec("GII.4", 0.20),
            LineageSpec("GII.2", 0.15),
            LineageSpec("GI.1", 0.15),
            LineageSpec("GII.6", 0.13),
            LineageSpec("GI.7", 0.12),
        ]
    ),
    "gi-dominant": MixtureSpec(
        lineages=[
            LineageSpec("GI.1", 0.40),
            LineageSpec("GI.3", 0.25),
            LineageSpec("GI.7", 0.20),
            LineageSpec("GII.17", 0.15),
        ]
    ),
    "outbreak": MixtureSpec(
        lineages=[
            LineageSpec("GII.17", 0.90),
            LineageSpec("GII.4", 0.10),
        ]
    ),
}


def resolve_preset(name: str) -> MixtureSpec:
    if name not in PRESETS:
        available = ", ".join(sorted(PRESETS.keys()))
        raise ValueError(f"Unknown preset '{name}'. Available: {available}")
    return PRESETS[name]


def resolve_mixture(
    genotype_col: str,
    cli_mixture: str | None,
    cli_preset: str | None,
) -> MixtureSpec:
    if ":" in genotype_col:
        return parse_mixture(genotype_col)
    if cli_mixture is not None:
        return parse_mixture(cli_mixture)
    if cli_preset is not None:
        return resolve_preset(cli_preset)
    return parse_mixture(genotype_col)
