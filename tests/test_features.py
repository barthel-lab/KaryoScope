"""Tests for :mod:`karyoscope.core.io.features`."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.core.io.features import (
    NOVEL_NAME,
    Features,
    FeaturesError,
    parse_features,
    render_feature,
)

_SAMPLE = """\
featureID\tchromosome\tregion\tgene
1\tchr1\trA\tintergenic
2\tchr1\trB\tgene_X
3\tchr2\trC\tintergenic
"""


def _write(tmp_path: Path, text: str, name: str = "features.tsv") -> Path:
    p = tmp_path / name
    p.write_text(text)
    return p


# --- happy path -------------------------------------------------------


def test_parse_minimal(tmp_path: Path) -> None:
    f = parse_features(_write(tmp_path, _SAMPLE))
    assert f.feature_sets == ["chromosome", "region", "gene"]
    assert set(f.table.keys()) == {1, 2, 3}


def test_lookup_returns_per_set_value(tmp_path: Path) -> None:
    f = parse_features(_write(tmp_path, _SAMPLE))
    assert f.feature_in_set(1, "chromosome") == "chr1"
    assert f.feature_in_set(1, "region") == "rA"
    assert f.feature_in_set(1, "gene") == "intergenic"
    assert f.feature_in_set(2, "region") == "rB"
    assert f.feature_in_set(3, "chromosome") == "chr2"


def test_lookup_returns_none_for_missing_feature_id(tmp_path: Path) -> None:
    f = parse_features(_write(tmp_path, _SAMPLE))
    assert f.feature_in_set(999, "chromosome") is None


def test_lookup_raises_on_unknown_feature_set(tmp_path: Path) -> None:
    f = parse_features(_write(tmp_path, _SAMPLE))
    with pytest.raises(FeaturesError, match="unknown feature set"):
        f.feature_in_set(1, "no_such_set")


def test_parse_tolerates_trailing_blank_lines(tmp_path: Path) -> None:
    f = parse_features(_write(tmp_path, _SAMPLE + "\n\n"))
    assert len(f.table) == 3


def test_parse_real_dummy_db(unpacked_dummy_db: Path) -> None:
    """The committed dummy db's features.tsv parses cleanly."""
    f = parse_features(unpacked_dummy_db / "features.tsv")
    assert f.feature_sets == ["chromosome", "region"]
    assert f.feature_in_set(1, "chromosome") == "chr1"
    assert f.feature_in_set(1, "region") == "rA"
    assert f.feature_in_set(3, "chromosome") == "chr2"
    assert f.feature_in_set(3, "region") == "rC"


# --- render_feature ---------------------------------------------------


def test_render_feature_id_zero_is_novel(tmp_path: Path) -> None:
    f = parse_features(_write(tmp_path, _SAMPLE))
    assert render_feature(0, "chromosome", f) == NOVEL_NAME
    assert render_feature(0, "gene", f) == NOVEL_NAME


def test_render_feature_known_id(tmp_path: Path) -> None:
    f = parse_features(_write(tmp_path, _SAMPLE))
    assert render_feature(1, "chromosome", f) == "chr1"
    assert render_feature(2, "region", f) == "rB"


def test_render_feature_missing_id_raises(tmp_path: Path) -> None:
    """A non-zero id absent from the table is a hard error.

    'Unknown' can be a legitimate feature name in real KaryoScope
    databases (e.g., in a repeats set), so silently mapping missing
    ids to 'Unknown' would be ambiguous and would hide genuine
    database / index mismatches.
    """
    f = parse_features(_write(tmp_path, _SAMPLE))
    with pytest.raises(FeaturesError, match="feature id 99 is not in features"):
        render_feature(99, "chromosome", f)


# --- error paths ------------------------------------------------------


def test_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FeaturesError, match="not found"):
        parse_features(tmp_path / "nope.tsv")


def test_empty_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FeaturesError, match="empty"):
        parse_features(_write(tmp_path, ""))


def test_wrong_first_column_name_raises(tmp_path: Path) -> None:
    bad = "id\tchromosome\tregion\n1\tchr1\trA\n"
    with pytest.raises(FeaturesError, match="featureID"):
        parse_features(_write(tmp_path, bad))


def test_header_without_feature_sets_raises(tmp_path: Path) -> None:
    bad = "featureID\n1\n"
    with pytest.raises(FeaturesError, match="no feature-set columns"):
        parse_features(_write(tmp_path, bad))


def test_duplicate_feature_set_in_header_raises(tmp_path: Path) -> None:
    bad = "featureID\tchromosome\tchromosome\n1\tchr1\tchr2\n"
    with pytest.raises(FeaturesError, match="duplicate feature-set"):
        parse_features(_write(tmp_path, bad))


def test_wrong_column_count_raises(tmp_path: Path) -> None:
    bad = "featureID\tchromosome\tregion\n1\tchr1\n"  # missing region column
    with pytest.raises(FeaturesError, match="tab-separated"):
        parse_features(_write(tmp_path, bad))


def test_non_integer_feature_id_raises(tmp_path: Path) -> None:
    bad = "featureID\tchromosome\nfoo\tchr1\n"
    with pytest.raises(FeaturesError, match="integer"):
        parse_features(_write(tmp_path, bad))


def test_feature_id_zero_in_file_raises(tmp_path: Path) -> None:
    """featureID 0 is reserved for the implicit 'novel' sentinel."""
    bad = "featureID\tchromosome\n0\tnovel\n"
    with pytest.raises(FeaturesError, match="featureID 0 is reserved"):
        parse_features(_write(tmp_path, bad))


def test_duplicate_feature_id_raises(tmp_path: Path) -> None:
    bad = "featureID\tchromosome\n1\tchr1\n1\tchr2\n"
    with pytest.raises(FeaturesError, match="duplicate featureID"):
        parse_features(_write(tmp_path, bad))


# --- direct construction ----------------------------------------------


def test_construct_features_directly() -> None:
    f = Features(
        feature_sets=["x", "y"],
        table={1: {"x": "a", "y": "b"}, 2: {"x": "c", "y": "d"}},
    )
    assert f.feature_in_set(1, "x") == "a"
    assert f.feature_in_set(2, "y") == "d"
