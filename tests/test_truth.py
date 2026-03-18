from pathlib import Path

from nocasim.fragment import Fragment
from nocasim.truth import (
    compute_mixture_coverage,
    compute_vp1_coverage,
    write_manifest,
    write_summary,
)


def _write_sam(path: Path, reads: list[tuple[int, str]]) -> None:
    with open(path, "w") as f:
        f.write("@HD\tVN:1.6\n")
        f.write("@SQ\tSN:ref\tLN:2000\n")
        for i, (pos, seq) in enumerate(reads):
            f.write(
                f"read_{i}\t0\tref\t{pos}\t60\t{len(seq)}M\t*\t0\t0\t{seq}\t{'I' * len(seq)}\n"
            )


def test_compute_vp1_coverage_basic(tmp_path):
    sam = tmp_path / "test.sam"
    seq = "A" * 100
    reads = [(501, seq)] * 25
    _write_sam(sam, reads)

    result = compute_vp1_coverage(sam, 500, 600)
    assert result["vp1_length_bp"] == 100
    assert result["vp1_mean_depth"] > 0


def test_compute_vp1_full_coverage(tmp_path):
    sam = tmp_path / "test.sam"
    reads = []
    for pos in range(100, 300, 5):
        reads.append((pos + 1, "A" * 50))
    reads = reads * 5
    _write_sam(sam, reads)

    result = compute_vp1_coverage(sam, 100, 300)
    assert result["vp1_breadth"] > 0.5


def test_write_manifest(tmp_path):
    stats = {"sample_id": "s1", "genotype": "GII.4", "ct_value": 25.0}
    path = tmp_path / "manifest.json"
    write_manifest(stats, path)
    import json

    data = json.loads(path.read_text())
    assert data["sample_id"] == "s1"


def test_write_summary(tmp_path):
    stats_list = [
        {
            "sample_id": "s1",
            "genotype": "GII.4",
            "ct_value": 25.0,
            "vp1_mean_depth": 100.0,
            "vp1_completeness_20x": 0.95,
            "completeness_call": "complete",
        },
        {
            "sample_id": "s2",
            "genotype": "GI.1",
            "ct_value": 35.0,
            "vp1_mean_depth": 5.0,
            "vp1_completeness_20x": 0.30,
            "completeness_call": "incomplete",
        },
    ]
    path = tmp_path / "summary.tsv"
    write_summary(stats_list, path)
    lines = path.read_text().strip().split("\n")
    assert len(lines) == 3
    assert "sample_id" in lines[0]
    assert "s1" in lines[1]
    assert "s2" in lines[2]


def _make_fragment(start, end, overlaps_vp1=True):
    return Fragment(
        id=f"frag_{start}_{end}",
        sequence="A" * (end - start),
        start=start,
        end=end,
        strand="+",
        gc=0.47,
        source="test",
        overlaps_vp1=overlaps_vp1,
    )


def test_compute_mixture_coverage_two_lineages():
    abundances = {"GII.17": 0.75, "GII.4": 0.25}
    lineage_data = {
        "GII.17": (
            100,
            200,
            [(_make_fragment(100, 200), 10)],
        ),
        "GII.4": (
            100,
            200,
            [(_make_fragment(100, 200), 5)],
        ),
    }
    result = compute_mixture_coverage(lineage_data, abundances)
    assert "aggregate" in result
    assert "per_lineage" in result
    assert "GII.17" in result["per_lineage"]
    assert "GII.4" in result["per_lineage"]
    assert result["per_lineage"]["GII.17"]["vp1_mean_depth"] > 0
    assert result["per_lineage"]["GII.4"]["vp1_mean_depth"] > 0


def test_compute_mixture_coverage_weighted_aggregate():
    # GII.17 at 75% with depth 100, GII.4 at 25% with depth 20
    # Weighted aggregate depth should be 0.75*100 + 0.25*20 = 80
    abundances = {"GII.17": 0.75, "GII.4": 0.25}
    lineage_data = {
        "GII.17": (
            0,
            100,
            [(_make_fragment(0, 100), 100)],
        ),
        "GII.4": (
            0,
            100,
            [(_make_fragment(0, 100), 20)],
        ),
    }
    result = compute_mixture_coverage(lineage_data, abundances)
    agg = result["aggregate"]
    assert abs(agg["vp1_mean_depth"] - 80.0) < 0.1


def test_compute_mixture_coverage_single_lineage():
    abundances = {"GII.4": 1.0}
    lineage_data = {
        "GII.4": (
            100,
            200,
            [(_make_fragment(120, 180), 20)],
        ),
    }
    result = compute_mixture_coverage(lineage_data, abundances)
    assert len(result["per_lineage"]) == 1
    assert result["aggregate"]["vp1_mean_depth"] > 0


def test_write_summary_with_lineage_detail(tmp_path):
    stats_list = [
        {
            "sample_id": "s1",
            "genotype": "GII.4",
            "ct_value": 25.0,
            "vp1_mean_depth": 100.0,
            "vp1_completeness_20x": 0.95,
            "completeness_call": "complete",
            "lineage_detail": "",
        },
        {
            "sample_id": "s2",
            "genotype": "GII.17:0.75,GII.4:0.25",
            "ct_value": 28.0,
            "vp1_mean_depth": 80.0,
            "vp1_completeness_20x": 0.90,
            "completeness_call": "complete",
            "lineage_detail": "GII.17:60.0x/0.88;GII.4:20.0x/0.72",
        },
    ]
    path = tmp_path / "summary.tsv"
    write_summary(stats_list, path)
    lines = path.read_text().strip().split("\n")
    assert "lineage_detail" in lines[0]
    assert len(lines) == 3
