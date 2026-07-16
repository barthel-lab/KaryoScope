"""Tests for :mod:`karyoscope.core.io.emit`."""

from __future__ import annotations

from pathlib import Path

from karyoscope.core.io import emit
from karyoscope.core.io.colors import parse_colors
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
