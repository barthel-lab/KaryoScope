"""Tests for :mod:`karyoscope.core.io.emit`."""

from __future__ import annotations

from pathlib import Path

import karyoscope.core.build as build_mod
from karyoscope.core.io import emit
from karyoscope.core.io.colors import parse_colors, parse_colors_and_groups
from karyoscope.core.io.hierarchy import parse_hierarchy
from karyoscope.manifest import parse_manifest


def test_auto_palette_distinct_and_deterministic() -> None:
    p1 = emit.auto_palette(5)
    p2 = emit.auto_palette(5)
    assert p1 == p2  # deterministic
    assert len(set(p1)) == 5  # distinct
    assert all(c.startswith("#") and len(c) == 7 for c in p1)
    assert emit.auto_palette(0) == []


def test_assign_colors_precedence() -> None:
    colors = emit.assign_colors(
        nodes={"categorized", "LINE", "SINE", "bg"},
        leaves=["LINE", "SINE", "bg"],
        background="bg",
        provided={"LINE": "#123456"},
    )
    assert colors["LINE"] == "#123456"  # provided wins
    assert colors["bg"] == emit.BACKGROUND_COLOR  # background grey
    assert colors["categorized"] == emit.STRUCTURAL_COLOR  # internal/root
    assert colors["SINE"] not in {"#123456", emit.BACKGROUND_COLOR, emit.STRUCTURAL_COLOR}


def test_writers_roundtrip_through_parsers(tmp_path: Path) -> None:
    emit.write_hierarchy_tsv(
        tmp_path / "hierarchy.tsv",
        [("repeat", "LINE", "categorized"), ("repeat", "bg", "categorized")],
    )
    hier = parse_hierarchy(tmp_path / "hierarchy.tsv")
    assert set(hier.nodes("repeat")) == {"LINE", "bg", "categorized"}

    emit.write_colors_tsv(
        tmp_path / "colors.tsv",
        [("repeat", "LINE", "#111111"), ("repeat", "bg", "#808080")],
    )
    colors = parse_colors(tmp_path / "colors.tsv")
    assert colors["repeat"]["LINE"] == "#111111"


def test_colors_tsv_stays_three_column_without_groups(tmp_path: Path) -> None:
    # Every database shipped before the optional column existed is 3-column,
    # and readers accept exactly 3 or exactly 4 -- so a build that doesn't ask
    # for grouping must keep producing the same 3-column file.
    p = tmp_path / "colors.tsv"
    emit.write_colors_tsv(p, [("repeat", "LINE", "#111111"), ("repeat", "bg", "#808080")])
    lines = p.read_text().splitlines()
    assert lines[0] == "feature_set\tfeature\tcolor"
    assert all(len(line.split("\t")) == 3 for line in lines)


def test_colors_tsv_emits_group_column_when_any_row_has_one(tmp_path: Path) -> None:
    p = tmp_path / "colors.tsv"
    emit.write_colors_tsv(
        p,
        [
            ("cytoband", "p11.1", "#000000", "gpos100"),
            ("cytoband", "q21.3", "#000000", "gpos100"),
            ("cytoband", "categorized", "#B0C4DE", ""),
        ],
    )
    lines = p.read_text().splitlines()
    assert lines[0] == "feature_set\tfeature\tcolor\tlegend_group"
    assert all(len(line.split("\t")) == 4 for line in lines)
    # Round-trips through the reader, groups intact.
    colors, groups = parse_colors_and_groups(p)
    assert colors["cytoband"]["p11.1"] == "#000000"
    assert groups["cytoband"] == {"p11.1": "gpos100", "q21.3": "gpos100"}


def test_colors_tsv_accepts_mixed_row_widths(tmp_path: Path) -> None:
    # 3- and 4-tuples together: the ungrouped row gets a blank 4th cell rather
    # than a ragged line that would fail the reader's column-count check.
    p = tmp_path / "colors.tsv"
    emit.write_colors_tsv(p, [("region", "rA", "#2ca02c"), ("region", "rB", "#d62728", "sat")])
    assert all(len(line.split("\t")) == 4 for line in p.read_text().splitlines())
    _, groups = parse_colors_and_groups(p)
    assert groups == {"region": {"rB": "sat"}}


def test_feature_set_hierarchy_txt_is_headerless_child_parent(tmp_path: Path) -> None:
    p = tmp_path / "features.repeat.hierarchy.txt"
    emit.write_feature_set_hierarchy_txt(p, [("LINE", "categorized"), ("SINE", "categorized")])
    assert p.read_text() == "LINE\tcategorized\nSINE\tcategorized\n"


def test_priorities_file_format(tmp_path: Path) -> None:
    p = tmp_path / "prio.txt"
    emit.write_priorities_file(p, {"categorized": 0, "LINE": 1})
    assert p.read_text() == "categorized 0\nLINE 1\n"


def test_manifest_omits_features_and_parses_as_hks(tmp_path: Path) -> None:
    m = emit.build_manifest_dict(
        db_id="HKS_x",
        version="1.0.0",
        karyoscope_min_version="1.1.0",
        basename="index/features",
        s=31,
        feature_sets=["repeat", "gene"],
        roles={"chromosome_assignment": "gene"},
    )
    assert "features" not in m
    emit.write_manifest(tmp_path / "manifest.yaml", m)
    parsed = parse_manifest(tmp_path / "manifest.yaml")
    assert parsed.features is None
    assert parsed.index.type == "hks"
    assert parsed.feature_sets == ["repeat", "gene"]
    assert parsed.roles == {"chromosome_assignment": "gene"}


# --- build's per-set colours file (the input side of the same column) ------


class TestBuildParseSetColors:
    """``build --colors NAME=file`` must accept the same optional 4th column
    that a database's ``colors.tsv`` carries -- otherwise a user can supply
    legend groups and have them silently dropped, which is what used to happen.
    """

    def test_two_column_form(self, tmp_path: Path) -> None:
        p = tmp_path / "c.txt"
        p.write_text("LINE\t#111111\nSINE\t#222222\n")
        colors, groups = build_mod._parse_set_colors(p)
        assert colors == {"LINE": "#111111", "SINE": "#222222"}
        assert groups == {}

    def test_three_column_form_with_header(self, tmp_path: Path) -> None:
        p = tmp_path / "c.txt"
        p.write_text("feature_set\tfeature\tcolor\nrepeat\tLINE\t#111111\n")
        colors, groups = build_mod._parse_set_colors(p)
        assert colors == {"LINE": "#111111"}
        assert groups == {}

    def test_four_column_form_keeps_groups(self, tmp_path: Path) -> None:
        p = tmp_path / "c.txt"
        p.write_text(
            "feature_set\tfeature\tcolor\tlegend_group\n"
            "cytoband\tp11.1\t#000000\tgpos100\n"
            "cytoband\tq21.3\t#000000\tgpos100\n"
            "cytoband\tp13.2\t#ffffff\tgneg\n"
        )
        colors, groups = build_mod._parse_set_colors(p)
        assert colors["p11.1"] == "#000000"
        assert groups == {"p11.1": "gpos100", "q21.3": "gpos100", "p13.2": "gneg"}

    def test_four_column_header_is_skipped(self, tmp_path: Path) -> None:
        # The header's last cell is 'legend_group', not 'color', so the old
        # header check would have parsed it as a feature named 'feature'.
        p = tmp_path / "c.txt"
        p.write_text(
            "feature_set\tfeature\tcolor\tlegend_group\ncytoband\tp11.1\t#000000\tgpos100\n"
        )
        colors, _ = build_mod._parse_set_colors(p)
        assert "feature" not in colors

    def test_blank_group_cell_is_ungrouped(self, tmp_path: Path) -> None:
        p = tmp_path / "c.txt"
        p.write_text(
            "feature_set\tfeature\tcolor\tlegend_group\nregion\trA\t#2ca02c\t\nregion\trB\t#d62728\tsat\n"
        )
        colors, groups = build_mod._parse_set_colors(p)
        assert set(colors) == {"rA", "rB"}
        assert groups == {"rB": "sat"}
