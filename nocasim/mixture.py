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
