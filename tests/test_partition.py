"""Tests for :mod:`karyoscope.core.io.partition`."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.core.io import partition
from karyoscope.core.io.partition import PartitionError


def _write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def test_parse_bed_basic(tmp_path: Path) -> None:
    bed = _write(tmp_path / "a.bed", "chr1\t0\t10\tLINE\nchr1\t10\t20\tSINE\n# c\n\n")
    assert partition.parse_bed(bed) == [
        ("chr1", 0, 10, "LINE"),
        ("chr1", 10, 20, "SINE"),
    ]


def test_parse_bed_rejects_short_rows(tmp_path: Path) -> None:
    bed = _write(tmp_path / "a.bed", "chr1\t0\t10\n")
    with pytest.raises(PartitionError, match="at least 4"):
        partition.parse_bed(bed)


def test_parse_bed_rejects_bad_coords(tmp_path: Path) -> None:
    bed = _write(tmp_path / "a.bed", "chr1\tx\t10\tL\n")
    with pytest.raises(PartitionError, match="integers"):
        partition.parse_bed(bed)


def test_read_fai(tmp_path: Path) -> None:
    fai = _write(tmp_path / "g.fa.fai", "chr1\t100\t6\t60\t61\nchr2\t50\t120\t60\t61\n")
    assert partition.read_fai(fai) == {"chr1": 100, "chr2": 50}


def test_labels_in_preserves_first_seen_order() -> None:
    ivs = [("c", 0, 1, "b"), ("c", 1, 2, "a"), ("c", 2, 3, "b")]
    assert partition.labels_in(ivs) == ["b", "a"]


def test_compute_background_fills_gaps_and_handles_overlap() -> None:
    fai = {"chr1": 40, "chr2": 20}
    ivs = [
        ("chr1", 0, 10, "LINE"),
        ("chr1", 5, 15, "SINE"),  # overlaps LINE; union covers [0,15)
        ("chr1", 20, 30, "LINE"),
    ]
    bg = partition.compute_background_intervals(ivs, fai, "bg")
    assert bg == [
        ("chr1", 15, 20, "bg"),
        ("chr1", 30, 40, "bg"),
        ("chr2", 0, 20, "bg"),  # unannotated chromosome fully filled
    ]


def test_compute_background_full_tiling_yields_nothing() -> None:
    fai = {"chr1": 10}
    ivs = [("chr1", 0, 10, "x")]
    assert partition.compute_background_intervals(ivs, fai, "bg") == []


def test_flatten_by_priority_resolves_overlap() -> None:
    ivs = [
        ("chr1", 0, 10, "LINE"),
        ("chr1", 5, 15, "SINE"),
    ]
    # SINE outranks LINE, so the overlap [5,10) becomes SINE.
    assert partition.flatten_by_priority(ivs, ["SINE", "LINE"]) == [
        ("chr1", 0, 5, "LINE"),
        ("chr1", 5, 15, "SINE"),
    ]
    # Reversed priority: LINE wins the overlap.
    assert partition.flatten_by_priority(ivs, ["LINE", "SINE"]) == [
        ("chr1", 0, 10, "LINE"),
        ("chr1", 10, 15, "SINE"),
    ]


def test_flatten_merges_adjacent_same_label() -> None:
    ivs = [("chr1", 0, 5, "a"), ("chr1", 5, 10, "a")]
    assert partition.flatten_by_priority(ivs, ["a"]) == [("chr1", 0, 10, "a")]


def _genome(tmp_path: Path) -> Path:
    # chr1 = 40 bp, chr2 = 20 bp
    return _write(
        tmp_path / "g.fa",
        ">chr1\n" + "ACGTACGTAC" * 4 + "\n>chr2 desc\n" + "GGCCTTAAGG" * 2 + "\n",
    )


def test_slice_extends_by_k_minus_1_and_clamps(tmp_path: Path) -> None:
    fa = _genome(tmp_path)
    ivs = [("chr1", 0, 10, "L"), ("chr1", 35, 40, "L"), ("chr2", 0, 20, "S")]
    paths = partition.slice_features_to_fastas(fa, ivs, k=4, outdir=tmp_path / "out")

    assert set(paths) == {"L", "S"}
    l_lines = paths["L"].read_text().splitlines()
    # [0,10) extended by k-1=3 -> 13 bp
    assert l_lines[0] == ">chr1:0-10_L"
    assert l_lines[1] == ("ACGTACGTAC" * 4)[0:13]
    # [35,40) extension clamps at contig length 40 -> 5 bp
    assert l_lines[2] == ">chr1:35-40_L"
    assert len(l_lines[3]) == 5
    # chr2 whole contig
    s_lines = paths["S"].read_text().splitlines()
    assert len(s_lines[1]) == 20


def test_slice_skips_absent_contigs(tmp_path: Path) -> None:
    fa = _genome(tmp_path)
    ivs = [("chrZ", 0, 5, "L")]
    paths = partition.slice_features_to_fastas(fa, ivs, k=4, outdir=tmp_path / "out")
    assert paths == {}
