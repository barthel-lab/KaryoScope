"""Unit tests for :mod:`karyoscope.core.io.colors`."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.core.io.colors import ColorsError, colors_for_set, parse_colors


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
