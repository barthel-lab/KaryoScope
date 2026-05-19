"""Tests for :mod:`karyoscope.core.io.hierarchy`."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.core.io.hierarchy import (
    Hierarchy,
    HierarchyError,
    HierarchyRow,
    parse_hierarchy,
)

_SAMPLE = """\
feature_set\tfeature\tparent
chromosome\tchr1\t.
chromosome\tchr2\t.
region\t1p\tchr1
region\t1q\tchr1
region\t2p\tchr2
"""


def _write(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "hierarchy.tsv"
    p.write_text(text)
    return p


# --- happy path --------------------------------------------------------


def test_parse_minimal(tmp_path: Path) -> None:
    h = parse_hierarchy(_write(tmp_path, _SAMPLE))
    assert len(h.rows) == 5
    assert h.rows[0] == HierarchyRow("chromosome", "chr1", None)
    assert h.rows[2] == HierarchyRow("region", "1p", "chr1")


def test_feature_sets_preserves_first_occurrence_order(tmp_path: Path) -> None:
    h = parse_hierarchy(_write(tmp_path, _SAMPLE))
    assert h.feature_sets() == ["chromosome", "region"]


def test_features_in(tmp_path: Path) -> None:
    h = parse_hierarchy(_write(tmp_path, _SAMPLE))
    assert h.features_in("chromosome") == ["chr1", "chr2"]
    assert h.features_in("region") == ["1p", "1q", "2p"]
    assert h.features_in("nope") == []


def test_count_by_feature_set(tmp_path: Path) -> None:
    h = parse_hierarchy(_write(tmp_path, _SAMPLE))
    assert h.count_by_feature_set() == {"chromosome": 2, "region": 3}


def test_roots_all(tmp_path: Path) -> None:
    h = parse_hierarchy(_write(tmp_path, _SAMPLE))
    roots = h.roots()
    assert {r.feature for r in roots} == {"chr1", "chr2"}


def test_roots_restricted_to_feature_set(tmp_path: Path) -> None:
    h = parse_hierarchy(_write(tmp_path, _SAMPLE))
    assert h.roots(feature_set="chromosome") == [
        HierarchyRow("chromosome", "chr1", None),
        HierarchyRow("chromosome", "chr2", None),
    ]
    assert h.roots(feature_set="region") == []


def test_parse_tolerates_blank_trailing_lines(tmp_path: Path) -> None:
    h = parse_hierarchy(_write(tmp_path, _SAMPLE + "\n\n"))
    assert len(h.rows) == 5


def test_parse_real_dummy_db(unpacked_dummy_db: Path) -> None:
    """The committed dummy db's hierarchy.tsv parses cleanly."""
    h = parse_hierarchy(unpacked_dummy_db / "hierarchy.tsv")
    assert "chromosome" in h.feature_sets()
    assert "region" in h.feature_sets()
    counts = h.count_by_feature_set()
    assert counts["chromosome"] == 2  # dummy db has chr1 and chr2
    assert counts["region"] == 3  # rA, rB, rC


# --- error paths -------------------------------------------------------


def test_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(HierarchyError, match="not found"):
        parse_hierarchy(tmp_path / "nope.tsv")


def test_empty_file_raises(tmp_path: Path) -> None:
    with pytest.raises(HierarchyError, match="empty"):
        parse_hierarchy(_write(tmp_path, ""))


def test_bad_header_raises(tmp_path: Path) -> None:
    bad = "wrong\theader\tcolumns\nchromosome\tchr1\t.\n"
    with pytest.raises(HierarchyError, match="unexpected header"):
        parse_hierarchy(_write(tmp_path, bad))


def test_wrong_column_count_raises(tmp_path: Path) -> None:
    bad = "feature_set\tfeature\tparent\nchromosome\tchr1\n"
    with pytest.raises(HierarchyError, match=r"3 tab-separated columns"):
        parse_hierarchy(_write(tmp_path, bad))


# --- Hierarchy direct construction ------------------------------------


def test_hierarchy_dataclass_round_trip() -> None:
    rows = [
        HierarchyRow("a", "x", None),
        HierarchyRow("a", "y", "x"),
        HierarchyRow("b", "z", None),
    ]
    h = Hierarchy(rows=rows)
    assert h.feature_sets() == ["a", "b"]
    assert h.features_in("a") == ["x", "y"]
    assert len(h.roots()) == 2
