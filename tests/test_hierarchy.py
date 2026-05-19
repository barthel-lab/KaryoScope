"""Tests for :mod:`karyoscope.core.io.hierarchy`."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.core.io.hierarchy import (
    REQUIRED_ROOT,
    Hierarchy,
    HierarchyError,
    HierarchyRow,
    parse_hierarchy,
    validate_hierarchy,
)

# A sample resembling the real production hierarchy in miniature.
_SAMPLE = """\
feature_set\tchild\tparent
chromosome\tautosome\tcategorized
chromosome\tchr1\tautosome
chromosome\tchr2\tautosome
region\tcentromeric\tcategorized
region\taSat\tcentromeric
region\trA\taSat
region\trB\taSat
"""


def _write(tmp_path: Path, text: str, name: str = "hierarchy.tsv") -> Path:
    p = tmp_path / name
    p.write_text(text)
    return p


# --- parser happy path ---------------------------------------------------


def test_parse_sample(tmp_path: Path) -> None:
    h = parse_hierarchy(_write(tmp_path, _SAMPLE))
    assert h.feature_sets() == ["chromosome", "region"]
    assert h.count_by_feature_set() == {"chromosome": 3, "region": 4}


def test_rows_in_returns_only_named_set(tmp_path: Path) -> None:
    h = parse_hierarchy(_write(tmp_path, _SAMPLE))
    chrom_rows = h.rows_in("chromosome")
    assert all(r.feature_set == "chromosome" for r in chrom_rows)
    assert len(chrom_rows) == 3


def test_parent_map(tmp_path: Path) -> None:
    h = parse_hierarchy(_write(tmp_path, _SAMPLE))
    p = h.parent_map("region")
    assert p == {
        "centromeric": "categorized",
        "aSat": "centromeric",
        "rA": "aSat",
        "rB": "aSat",
    }


def test_nodes_includes_root(tmp_path: Path) -> None:
    """``nodes()`` should include the root, which only appears as a parent."""
    h = parse_hierarchy(_write(tmp_path, _SAMPLE))
    assert REQUIRED_ROOT in h.nodes("chromosome")
    assert REQUIRED_ROOT in h.nodes("region")


def test_parse_tolerates_blank_lines(tmp_path: Path) -> None:
    h = parse_hierarchy(_write(tmp_path, _SAMPLE + "\n\n"))
    assert sum(h.count_by_feature_set().values()) == 7


def test_parse_dummy_db_file(unpacked_dummy_db: Path) -> None:
    """The committed dummy db's hierarchy.tsv parses cleanly."""
    h = parse_hierarchy(unpacked_dummy_db / "hierarchy.tsv")
    assert "chromosome" in h.feature_sets()
    assert "region" in h.feature_sets()
    # dummy db's region: rA, rB → aSat → centromeric → categorized;
    #                    rC     → HSat → centromeric → categorized
    p = h.parent_map("region")
    assert p["rA"] == "aSat"
    assert p["rB"] == "aSat"
    assert p["rC"] == "HSat"
    assert p["aSat"] == "centromeric"
    assert p["HSat"] == "centromeric"
    assert p["centromeric"] == "categorized"


# --- parser error paths -------------------------------------------------


def test_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(HierarchyError, match="not found"):
        parse_hierarchy(tmp_path / "nope.tsv")


def test_empty_file_raises(tmp_path: Path) -> None:
    with pytest.raises(HierarchyError, match="empty"):
        parse_hierarchy(_write(tmp_path, ""))


def test_bad_header_raises(tmp_path: Path) -> None:
    bad = "feature_set\tfeature\tparent\nchromosome\tchr1\tcategorized\n"
    with pytest.raises(HierarchyError, match="unexpected header"):
        parse_hierarchy(_write(tmp_path, bad))


def test_wrong_column_count_raises(tmp_path: Path) -> None:
    bad = "feature_set\tchild\tparent\nchromosome\tchr1\n"
    with pytest.raises(HierarchyError, match="3 tab-separated"):
        parse_hierarchy(_write(tmp_path, bad))


def test_blank_column_raises(tmp_path: Path) -> None:
    bad = "feature_set\tchild\tparent\nchromosome\t\tcategorized\n"
    with pytest.raises(HierarchyError, match="non-empty"):
        parse_hierarchy(_write(tmp_path, bad))


# --- validator: happy path ----------------------------------------------


def test_validate_clean_sample(tmp_path: Path) -> None:
    h = parse_hierarchy(_write(tmp_path, _SAMPLE))
    assert validate_hierarchy(h) == []


# --- validator: error detection -----------------------------------------


def test_validate_detects_multiple_parents() -> None:
    """A child appearing in two rows within one set is flagged."""
    h = Hierarchy(
        rows=[
            HierarchyRow("chromosome", "chr1", "categorized"),
            HierarchyRow("chromosome", "chr1", "autosome"),  # duplicate
        ]
    )
    issues = validate_hierarchy(h)
    assert any("multiple parent rows" in i for i in issues)


def test_validate_detects_wrong_root_name() -> None:
    """A root that isn't named 'categorized' is flagged."""
    h = Hierarchy(
        rows=[
            HierarchyRow("chromosome", "chr1", "wrongroot"),
            HierarchyRow("chromosome", "chr2", "wrongroot"),
        ]
    )
    issues = validate_hierarchy(h)
    assert any("must be 'categorized'" in i for i in issues)


def test_validate_detects_multiple_roots() -> None:
    """Two roots within one feature set is flagged."""
    h = Hierarchy(
        rows=[
            HierarchyRow("chromosome", "chr1", "categorized"),
            HierarchyRow("chromosome", "chr2", "other_root"),
        ]
    )
    issues = validate_hierarchy(h)
    assert any("multiple roots" in i for i in issues)


def test_validate_detects_cycle() -> None:
    """A cycle (a → b → a) is flagged."""
    h = Hierarchy(
        rows=[
            HierarchyRow("region", "a", "b"),
            HierarchyRow("region", "b", "a"),
        ]
    )
    issues = validate_hierarchy(h)
    # The cycle prevents finding any root, so we get "no root" rather
    # than a cycle message — either is fine, but at least one issue
    # must fire.
    assert issues


def test_validate_detects_features_tsv_inconsistency() -> None:
    """A feature in features.tsv that has no row in hierarchy is flagged."""
    h = Hierarchy(
        rows=[
            HierarchyRow("region", "rA", "categorized"),
        ]
    )
    feature_columns = {"region": {"rA", "rB"}}  # rB has no hierarchy node
    issues = validate_hierarchy(h, feature_columns=feature_columns)
    assert any("'rB'" in i and "no row" in i for i in issues)


def test_validate_separate_sets_are_independent() -> None:
    """Errors in one feature set don't suppress validation of another."""
    h = Hierarchy(
        rows=[
            # chromosome: well-formed
            HierarchyRow("chromosome", "chr1", "categorized"),
            # region: malformed (no categorized root)
            HierarchyRow("region", "rA", "wrongroot"),
        ]
    )
    issues = validate_hierarchy(h)
    assert any("region" in i and "wrongroot" in i for i in issues)
    # No spurious complaints about the chromosome set
    assert not any(i.startswith("feature set 'chromosome'") for i in issues)


def test_validate_real_production_hierarchy(tmp_path: Path) -> None:
    """The real production CHM13 hierarchy validates clean.

    This is the hierarchy file Rhyker shared; if our parser/validator
    breaks it, we've broken the contract.
    """
    real_hierarchy = """\
feature_set\tchild\tparent
acrocentric\tacrocentric\tcategorized
acrocentric\tarray\tacrocentric
acrocentric\tDJ\tarray
acrocentric\tPJ\tarray
acrocentric\trDNA\tarray
acrocentric\tSST1\tacrocentric
acrocentric\tPHR\tacrocentric
acrocentric\tnonacrocentric\tcategorized
chromosome\tautosome\tcategorized
chromosome\tacrocentric\tautosome
chromosome\tchr13\tacrocentric
chromosome\tchr14\tacrocentric
chromosome\tchr1\tautosome
chromosome\tsex\tcategorized
chromosome\tchrX\tsex
chromosome\tchrY\tsex
gene\texon\tcategorized
gene\tintron\tcategorized
gene\tintergenic\tcategorized
"""
    h = parse_hierarchy(_write(tmp_path, real_hierarchy))
    issues = validate_hierarchy(h)
    assert issues == [], f"unexpected issues: {issues}"
