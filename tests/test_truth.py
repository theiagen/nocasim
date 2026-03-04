from pathlib import Path

from nocasim.truth import compute_vp1_coverage, write_manifest, write_summary


def _write_sam(path: Path, reads: list[tuple[int, str]]) -> None:
    with open(path, "w") as f:
        f.write("@HD\tVN:1.6\n")
        f.write("@SQ\tSN:ref\tLN:2000\n")
        for i, (pos, seq) in enumerate(reads):
            f.write(f"read_{i}\t0\tref\t{pos}\t60\t{len(seq)}M\t*\t0\t0\t{seq}\t{'I'*len(seq)}\n")


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
        {"sample_id": "s1", "genotype": "GII.4", "ct_value": 25.0,
         "vp1_mean_depth": 100.0, "vp1_completeness_20x": 0.95,
         "completeness_call": "complete"},
        {"sample_id": "s2", "genotype": "GI.1", "ct_value": 35.0,
         "vp1_mean_depth": 5.0, "vp1_completeness_20x": 0.30,
         "completeness_call": "incomplete"},
    ]
    path = tmp_path / "summary.tsv"
    write_summary(stats_list, path)
    lines = path.read_text().strip().split("\n")
    assert len(lines) == 3
    assert "sample_id" in lines[0]
    assert "s1" in lines[1]
    assert "s2" in lines[2]
