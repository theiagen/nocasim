import csv
import json
from pathlib import Path

from nocasim.fragment import Fragment


def compute_vp1_coverage_from_fragments(
    fragments: list[tuple[Fragment, int]],
    vp1_start: int,
    vp1_end: int,
) -> dict:
    vp1_len = vp1_end - vp1_start
    depth = [0] * vp1_len

    for frag, copy_count in fragments:
        if not frag.overlaps_vp1:
            continue
        overlap_start = max(frag.start, vp1_start) - vp1_start
        overlap_end = min(frag.end, vp1_end) - vp1_start
        if overlap_start < overlap_end:
            for i in range(overlap_start, overlap_end):
                if 0 <= i < vp1_len:
                    depth[i] += copy_count

    return _depth_to_stats(depth, vp1_len)


def _depth_to_stats(depth: list[int], vp1_len: int) -> dict:
    covered_20x = sum(1 for d in depth if d >= 20)
    covered_any = sum(1 for d in depth if d > 0)
    mean_depth = sum(depth) / vp1_len if vp1_len > 0 else 0.0
    completeness_20x = covered_20x / vp1_len if vp1_len > 0 else 0.0
    breadth = covered_any / vp1_len if vp1_len > 0 else 0.0

    if completeness_20x >= 0.90:
        call = "complete"
    elif breadth >= 0.90:
        call = "low_coverage"
    else:
        call = "incomplete"

    return {
        "vp1_length_bp": vp1_len,
        "vp1_mean_depth": round(mean_depth, 1),
        "vp1_completeness_20x": round(completeness_20x, 3),
        "vp1_breadth": round(breadth, 3),
        "completeness_call": call,
    }


def compute_vp1_coverage(sam_path: Path, vp1_start: int, vp1_end: int) -> dict:
    vp1_len = vp1_end - vp1_start
    depth = [0] * vp1_len

    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            flag = int(fields[1])
            if flag & 4:
                continue
            pos = int(fields[3]) - 1
            seq = fields[9]
            read_end = pos + len(seq)

            overlap_start = max(pos, vp1_start) - vp1_start
            overlap_end = min(read_end, vp1_end) - vp1_start

            if overlap_start < overlap_end:
                for i in range(overlap_start, overlap_end):
                    if 0 <= i < vp1_len:
                        depth[i] += 1

    covered_20x = sum(1 for d in depth if d >= 20)
    covered_any = sum(1 for d in depth if d > 0)
    mean_depth = sum(depth) / vp1_len if vp1_len > 0 else 0.0
    completeness_20x = covered_20x / vp1_len if vp1_len > 0 else 0.0
    breadth = covered_any / vp1_len if vp1_len > 0 else 0.0

    if completeness_20x >= 0.90:
        call = "complete"
    elif breadth >= 0.90:
        call = "low_coverage"
    else:
        call = "incomplete"

    return {
        "vp1_length_bp": vp1_len,
        "vp1_mean_depth": round(mean_depth, 1),
        "vp1_completeness_20x": round(completeness_20x, 3),
        "vp1_breadth": round(breadth, 3),
        "completeness_call": call,
    }


def compute_mixture_coverage(
    lineage_data: dict[str, tuple[int, int, list[tuple[Fragment, int]]]],
    abundances: dict[str, float],
) -> dict:
    per_lineage = {}
    for genotype, (vp1_start, vp1_end, frags) in lineage_data.items():
        per_lineage[genotype] = compute_vp1_coverage_from_fragments(
            frags, vp1_start, vp1_end
        )

    aggregate = _aggregate_lineage_stats(per_lineage, abundances)
    return {"aggregate": aggregate, "per_lineage": per_lineage}


def _aggregate_lineage_stats(
    per_lineage: dict[str, dict],
    abundances: dict[str, float],
) -> dict:
    if len(per_lineage) == 1:
        return dict(next(iter(per_lineage.values())))

    mean_depth = sum(
        s["vp1_mean_depth"] * abundances[g] for g, s in per_lineage.items()
    )
    completeness = sum(
        s["vp1_completeness_20x"] * abundances[g] for g, s in per_lineage.items()
    )
    breadth = sum(s["vp1_breadth"] * abundances[g] for g, s in per_lineage.items())
    vp1_length = max(s["vp1_length_bp"] for s in per_lineage.values())

    if completeness >= 0.90:
        call = "complete"
    elif breadth >= 0.90:
        call = "low_coverage"
    else:
        call = "incomplete"

    return {
        "vp1_length_bp": vp1_length,
        "vp1_mean_depth": round(mean_depth, 1),
        "vp1_completeness_20x": round(completeness, 3),
        "vp1_breadth": round(breadth, 3),
        "completeness_call": call,
    }


def write_manifest(stats: dict, path: Path) -> None:
    with open(path, "w") as f:
        json.dump(stats, f, indent=2)


def write_summary(all_stats: list[dict], path: Path) -> None:
    if not all_stats:
        return

    fieldnames = [
        "sample_id",
        "genotype",
        "ct_value",
        "vp1_mean_depth",
        "vp1_completeness_20x",
        "completeness_call",
        "lineage_detail",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for stats in all_stats:
            writer.writerow(stats)
