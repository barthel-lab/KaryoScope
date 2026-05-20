"""Unit tests for :mod:`karyoscope.core.io.colors`."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.core.io.colors import (
    ColorsError,
    colors_for_set,
    parse_colors,
    validate_colors,
)
from karyoscope.core.io.hierarchy import Hierarchy, HierarchyRow


class TestParseColors:
    def test_basic(self, tmp_path: Path) -> None:
        p = tmp_path / "colors.txt"
        p.write_text(
            "feature_set\tfeature\tcolor\n"
            "chromosome\tchr1\t#1f77b4\n"
            "chromosome\tchr2\t#ff7f0e\n"
            "region\trA\t#2ca02c\n"
        )
        c = parse_colors(p)
        assert c == {
            "chromosome": {"chr1": "#1f77b4", "chr2": "#ff7f0e"},
            "region": {"rA": "#2ca02c"},
        }

    def test_same_feature_different_sets(self, tmp_path: Path) -> None:
        # Same feature name in two feature sets gets two distinct colours.
        p = tmp_path / "colors.txt"
        p.write_text(
            "feature_set\tfeature\tcolor\n"
            "chromosome\tacrocentric\t#111111\n"
            "region\tacrocentric\t#222222\n"
        )
        c = parse_colors(p)
        assert c["chromosome"]["acrocentric"] == "#111111"
        assert c["region"]["acrocentric"] == "#222222"

    def test_missing_file(self, tmp_path: Path) -> None:
        with pytest.raises(ColorsError, match="not found"):
            parse_colors(tmp_path / "missing.txt")

    def test_empty_file(self, tmp_path: Path) -> None:
        p = tmp_path / "empty.txt"
        p.write_text("")
        with pytest.raises(ColorsError, match="empty"):
            parse_colors(p)

    def test_bad_header(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.txt"
        p.write_text("wrong\tcolumns\there\nchromosome\tchr1\t#000000\n")
        with pytest.raises(ColorsError, match="header"):
            parse_colors(p)

    def test_wrong_column_count(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.txt"
        p.write_text("feature_set\tfeature\tcolor\nchromosome\tchr1\n")
        with pytest.raises(ColorsError, match="3 tab-separated columns"):
            parse_colors(p)

    def test_blank_value_rejected(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.txt"
        p.write_text("feature_set\tfeature\tcolor\nchromosome\tchr1\t\n")
        with pytest.raises(ColorsError, match="non-empty"):
            parse_colors(p)


class TestColorsForSet:
    def test_novel_sentinel_present(self) -> None:
        c = {"region": {"rA": "#2ca02c"}}
        out = colors_for_set(c, "region")
        assert out["novel"] == "#ffffff"
        assert out["rA"] == "#2ca02c"

    def test_user_can_override_novel(self) -> None:
        # If the user genuinely declares a colour for "novel", it wins.
        c = {"region": {"rA": "#2ca02c", "novel": "#abcdef"}}
        out = colors_for_set(c, "region")
        assert out["novel"] == "#abcdef"

    def test_unknown_feature_set_only_has_novel(self) -> None:
        out = colors_for_set({"region": {"rA": "#000"}}, "does_not_exist")
        assert out == {"novel": "#ffffff"}


class TestValidateColors:
    """``validate_colors`` cross-checks hierarchy.tsv against colors.tsv.

    Every node in the hierarchy (children + parents, so including
    the ``categorized`` root which only appears as a parent) must
    have a colour. The ``novel`` sentinel is exempt -- it's always
    rendered white by the renderer.
    """

    def _region_hierarchy(self) -> Hierarchy:
        return Hierarchy(
            rows=[
                HierarchyRow("region", "centromeric", "categorized"),
                HierarchyRow("region", "aSat", "centromeric"),
                HierarchyRow("region", "rA", "aSat"),
            ]
        )

    def test_complete_colors_pass(self) -> None:
        h = self._region_hierarchy()
        colors = {
            "region": {
                "categorized": "#FFF",
                "centromeric": "#EEE",
                "aSat": "#DDD",
                "rA": "#CCC",
            }
        }
        assert validate_colors(h, colors) == []

    def test_missing_leaf_flagged(self) -> None:
        h = self._region_hierarchy()
        colors = {
            "region": {
                "categorized": "#FFF",
                "centromeric": "#EEE",
                "aSat": "#DDD",
                # rA missing
            }
        }
        issues = validate_colors(h, colors)
        assert len(issues) == 1
        assert "rA" in issues[0]

    def test_missing_categorized_flagged(self) -> None:
        # `categorized` only appears as a parent (root). Validation
        # must still require it since smoothing can promote to it.
        h = self._region_hierarchy()
        colors = {
            "region": {
                # categorized missing
                "centromeric": "#EEE",
                "aSat": "#DDD",
                "rA": "#CCC",
            }
        }
        issues = validate_colors(h, colors)
        assert len(issues) == 1
        assert "categorized" in issues[0]

    def test_missing_top_level_grouping_flagged(self) -> None:
        # Mirrors the production CHM13 bug: repeat/repeat or
        # acrocentric/acrocentric where the top-level group's name
        # matches the feature set name and was accidentally omitted
        # from colors.tsv.
        h = Hierarchy(
            rows=[
                HierarchyRow("repeat", "repeat", "categorized"),
                HierarchyRow("repeat", "LINE", "repeat"),
            ]
        )
        colors = {
            "repeat": {
                "categorized": "#FFF",
                "LINE": "#00F",
                # repeat missing
            }
        }
        issues = validate_colors(h, colors)
        assert any("'repeat'" in issue for issue in issues)

    def test_missing_feature_set_entirely_flagged(self) -> None:
        h = self._region_hierarchy()
        colors: dict[str, dict[str, str]] = {}  # no feature sets at all
        issues = validate_colors(h, colors)
        # Every node should be flagged.
        assert len(issues) == 4  # categorized, centromeric, aSat, rA

    def test_novel_in_hierarchy_is_exempt(self) -> None:
        # Defensive: if "novel" somehow appears in hierarchy.tsv
        # (it shouldn't, but the parser doesn't reject it), it's
        # exempt because the renderer special-cases it to white.
        h = Hierarchy(
            rows=[
                HierarchyRow("region", "rA", "categorized"),
                HierarchyRow("region", "novel", "categorized"),
            ]
        )
        colors = {"region": {"categorized": "#FFF", "rA": "#CCC"}}
        # "novel" is missing from colors but should NOT be flagged.
        issues = validate_colors(h, colors)
        assert issues == []

    def test_per_feature_set_isolation(self) -> None:
        # Same feature name in two sets with one missing in only one.
        h = Hierarchy(
            rows=[
                HierarchyRow("region", "shared", "categorized"),
                HierarchyRow("chromosome", "shared", "categorized"),
            ]
        )
        colors = {
            "region": {"categorized": "#FFF", "shared": "#CCC"},
            "chromosome": {"categorized": "#FFF"},  # `shared` missing here
        }
        issues = validate_colors(h, colors)
        assert len(issues) == 1
        assert "'chromosome'" in issues[0]
        assert "'shared'" in issues[0]
