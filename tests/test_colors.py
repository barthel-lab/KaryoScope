"""Unit tests for :mod:`karyoscope.core.io.colors`."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.core.io.colors import (
    ColorsError,
    colors_for_set,
    parse_colors,
    parse_colors_and_groups,
    validate_colors,
    validate_legend_groups,
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


class TestParseColorsAndGroups:
    """The optional 4th ``legend_group`` column.

    Databases shipped before the column existed are 3-column, so the parser
    must accept exactly 3 or exactly 4 and nothing in between -- a 3-column
    file has to keep parsing byte-identically.
    """

    def test_three_column_file_declares_no_groups(self, tmp_path: Path) -> None:
        p = tmp_path / "colors.txt"
        p.write_text(
            "feature_set\tfeature\tcolor\nchromosome\tchr1\t#1f77b4\nchromosome\tchr2\t#ff7f0e\n"
        )
        colors, groups = parse_colors_and_groups(p)
        assert colors == {"chromosome": {"chr1": "#1f77b4", "chr2": "#ff7f0e"}}
        assert groups == {}

    def test_four_column_file_collects_groups(self, tmp_path: Path) -> None:
        p = tmp_path / "colors.txt"
        p.write_text(
            "feature_set\tfeature\tcolor\tlegend_group\n"
            "cytoband\tp11.1\t#000000\tgpos100\n"
            "cytoband\tq21.3\t#000000\tgpos100\n"
            "cytoband\tp13.2\t#ffffff\tgneg\n"
        )
        colors, groups = parse_colors_and_groups(p)
        assert colors["cytoband"] == {
            "p11.1": "#000000",
            "q21.3": "#000000",
            "p13.2": "#ffffff",
        }
        assert groups == {"cytoband": {"p11.1": "gpos100", "q21.3": "gpos100", "p13.2": "gneg"}}

    def test_parse_colors_ignores_the_group_column(self, tmp_path: Path) -> None:
        # The many callers that only want colours must be unaffected by the
        # column's presence.
        p = tmp_path / "colors.txt"
        p.write_text(
            "feature_set\tfeature\tcolor\tlegend_group\ncytoband\tp11.1\t#000000\tgpos100\n"
        )
        assert parse_colors(p) == {"cytoband": {"p11.1": "#000000"}}

    def test_blank_group_is_not_a_group(self, tmp_path: Path) -> None:
        # An empty 4th cell means "ungrouped", not a group named "".
        p = tmp_path / "colors.txt"
        p.write_text(
            "feature_set\tfeature\tcolor\tlegend_group\n"
            "region\trA\t#2ca02c\t\n"
            "region\trB\t#d62728\tsat\n"
        )
        _, groups = parse_colors_and_groups(p)
        assert groups == {"region": {"rB": "sat"}}

    def test_five_column_header_rejected(self, tmp_path: Path) -> None:
        p = tmp_path / "colors.txt"
        p.write_text(
            "feature_set\tfeature\tcolor\tlegend_group\textra\nregion\trA\t#2ca02c\tsat\tx\n"
        )
        with pytest.raises(ColorsError, match="unexpected header"):
            parse_colors_and_groups(p)

    def test_row_must_match_the_declared_column_count(self, tmp_path: Path) -> None:
        # A 4-column header with a 3-column row is a malformed file, not an
        # implicitly ungrouped feature.
        p = tmp_path / "colors.txt"
        p.write_text("feature_set\tfeature\tcolor\tlegend_group\nregion\trA\t#2ca02c\n")
        with pytest.raises(ColorsError, match="expected 4 tab-separated columns"):
            parse_colors_and_groups(p)


class TestValidateLegendGroups:
    """One legend row is one swatch, so a group must not span two colours."""

    def test_consistent_groups_pass(self) -> None:
        colors = {"cytoband": {"p11.1": "#000000", "q21.3": "#000000"}}
        groups = {"cytoband": {"p11.1": "gpos100", "q21.3": "gpos100"}}
        assert validate_legend_groups(colors, groups) == []

    def test_three_column_file_has_nothing_to_validate(self) -> None:
        assert validate_legend_groups({"region": {"rA": "#2ca02c"}}, {}) == []

    def test_group_spanning_two_colours_flagged(self) -> None:
        colors = {"cytoband": {"p11.1": "#000000", "q21.3": "#ffffff"}}
        groups = {"cytoband": {"p11.1": "gpos100", "q21.3": "gpos100"}}
        issues = validate_legend_groups(colors, groups)
        assert len(issues) == 1
        assert "gpos100" in issues[0]
        assert "#000000" in issues[0] and "#ffffff" in issues[0]

    def test_two_groups_may_share_a_colour(self) -> None:
        # Legible (two labels, one swatch), just redundant -- deliberately allowed.
        colors = {"cytoband": {"a": "#000000", "b": "#000000"}}
        groups = {"cytoband": {"a": "gpos100", "b": "gpos75"}}
        assert validate_legend_groups(colors, groups) == []

    def test_missing_colour_is_not_this_validator_s_business(self) -> None:
        # validate_colors reports absent colours; this one must not double-report.
        colors: dict[str, dict[str, str]] = {"cytoband": {}}
        groups = {"cytoband": {"p11.1": "gpos100"}}
        assert validate_legend_groups(colors, groups) == []


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


class TestCustomColorsFile:
    """A user-supplied ``karyotype --colors`` file is parsed + validated exactly
    like a database ``colors.tsv``, so it fails *gracefully* (a clean
    ``ColorsError`` / ``KaryotypeError``, caught by the CLI) rather than crashing
    mid-render when it is malformed or missing features."""

    def test_malformed_custom_file_raises_colorserror(self, tmp_path: Path) -> None:
        p = tmp_path / "custom.tsv"
        # Missing the colour column -> rejected before any rendering.
        p.write_text("feature_set\tfeature\tcolor\ncytoband\t14q12\n")
        with pytest.raises(ColorsError):
            parse_colors(p)

    def test_custom_file_missing_features_is_flagged(self, tmp_path: Path) -> None:
        # A custom palette that omits some hierarchy nodes is caught by
        # validate_colors (which karyotype_run turns into a KaryotypeError).
        h = Hierarchy(
            rows=[
                HierarchyRow("cytoband", "14q12", "categorized"),
                HierarchyRow("cytoband", "21q21.1", "categorized"),
            ]
        )
        p = tmp_path / "custom.tsv"
        p.write_text(
            "feature_set\tfeature\tcolor\n"
            "cytoband\tcategorized\t#FFFFFF\n"
            "cytoband\t14q12\t#a67c5b\n"  # 21q21.1 deliberately omitted
        )
        issues = validate_colors(h, parse_colors(p))
        assert any("21q21.1" in issue for issue in issues)
