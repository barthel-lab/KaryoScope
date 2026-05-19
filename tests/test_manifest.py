"""Tests for :mod:`karyoscope.manifest`."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.exceptions import DatabaseLayoutError, ManifestError
from karyoscope.manifest import parse_manifest, validate_database_layout


def _write_minimal_db(dir_: Path) -> Path:
    """Build a minimal, valid KaryoScope database layout. Returns dir_."""
    dir_.mkdir(parents=True, exist_ok=True)
    (dir_ / "manifest.yaml").write_text(
        "id: test_db\n"
        "version: 1.0.0\n"
        "karyoscope_min_version: 0.1.0\n"
        "index:\n"
        "  type: kmc\n"
        "  basename: index/features\n"
        "hierarchy: hierarchy.tsv\n"
        "features: features.tsv\n"
        "colors: colors.txt\n"
        "kmer:\n"
        "  size: 21\n"
        "  type: fixed\n"
        "  max_size: 21\n"
        "feature_sets: [chromosome]\n"
    )
    (dir_ / "hierarchy.tsv").write_text("feature_set\tfeature\tparent\n")
    (dir_ / "features.tsv").write_text("feature_set\tfeature\tfeature_id\n")
    (dir_ / "colors.txt").write_text("feature_set\tfeature\tcolor\n")
    (dir_ / "index").mkdir(exist_ok=True)
    (dir_ / "index" / "features.kmc_pre").write_bytes(b"\x00" * 8)
    (dir_ / "index" / "features.kmc_suf").write_bytes(b"\x00" * 8)
    return dir_


# --- parse_manifest -----------------------------------------------------


def test_parse_manifest_minimal(tmp_path: Path) -> None:
    db = _write_minimal_db(tmp_path / "db")
    m = parse_manifest(db / "manifest.yaml")

    assert m.id == "test_db"
    assert m.version == "1.0.0"
    assert m.index.type == "kmc"
    assert m.index.basename == "index/features"
    assert m.kmer.size == 21
    assert m.feature_sets == ["chromosome"]


def test_parse_manifest_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(ManifestError, match="not found"):
        parse_manifest(tmp_path / "nope.yaml")


@pytest.mark.parametrize(
    "field",
    [
        "id",
        "version",
        "karyoscope_min_version",
        "index",
        "hierarchy",
        "features",
        "colors",
        "kmer",
        "feature_sets",
    ],
)
def test_parse_manifest_missing_required_field(tmp_path: Path, field: str) -> None:
    db = _write_minimal_db(tmp_path / "db")
    text = (db / "manifest.yaml").read_text()
    # Remove the line that starts with `<field>:`.
    new_text = "\n".join(line for line in text.splitlines() if not line.startswith(f"{field}:"))
    (db / "manifest.yaml").write_text(new_text)

    with pytest.raises(ManifestError):
        parse_manifest(db / "manifest.yaml")


def test_parse_manifest_unsupported_index_type(tmp_path: Path) -> None:
    db = _write_minimal_db(tmp_path / "db")
    (db / "manifest.yaml").write_text(
        (db / "manifest.yaml").read_text().replace("type: kmc", "type: nonsense")
    )
    with pytest.raises(ManifestError, match="unsupported index type"):
        parse_manifest(db / "manifest.yaml")


def test_parse_manifest_bad_kmer_type(tmp_path: Path) -> None:
    db = _write_minimal_db(tmp_path / "db")
    (db / "manifest.yaml").write_text(
        (db / "manifest.yaml").read_text().replace("type: fixed", "type: wibble")
    )
    with pytest.raises(ManifestError, match=r"kmer\.type"):
        parse_manifest(db / "manifest.yaml")


def test_parse_manifest_empty_feature_sets(tmp_path: Path) -> None:
    db = _write_minimal_db(tmp_path / "db")
    (db / "manifest.yaml").write_text(
        (db / "manifest.yaml").read_text().replace("feature_sets: [chromosome]", "feature_sets: []")
    )
    with pytest.raises(ManifestError, match="feature_sets"):
        parse_manifest(db / "manifest.yaml")


# --- validate_database_layout ------------------------------------------


def test_validate_database_layout_succeeds_on_minimal(tmp_path: Path) -> None:
    db = _write_minimal_db(tmp_path / "db")
    m = validate_database_layout(db)
    assert m.id == "test_db"
    assert m.directory == db.resolve()


def test_validate_database_layout_missing_hierarchy(tmp_path: Path) -> None:
    db = _write_minimal_db(tmp_path / "db")
    (db / "hierarchy.tsv").unlink()
    with pytest.raises(DatabaseLayoutError, match="hierarchy"):
        validate_database_layout(db)


def test_validate_database_layout_missing_kmc_suf(tmp_path: Path) -> None:
    db = _write_minimal_db(tmp_path / "db")
    (db / "index" / "features.kmc_suf").unlink()
    with pytest.raises(DatabaseLayoutError, match=r"\.kmc_suf"):
        validate_database_layout(db)


def test_validate_database_layout_path_traversal(tmp_path: Path) -> None:
    db = _write_minimal_db(tmp_path / "db")
    # Replace the colors path with one that escapes the db dir.
    manifest_path = db / "manifest.yaml"
    text = manifest_path.read_text().replace(
        "colors: colors.txt",
        "colors: ../escape.txt",
    )
    manifest_path.write_text(text)
    # Create the escape target so the test isn't fooled by file-not-found.
    (db.parent / "escape.txt").write_text("x")

    with pytest.raises(DatabaseLayoutError, match="escapes database directory"):
        validate_database_layout(db)


def test_validate_real_dummy_db(unpacked_dummy_db: Path) -> None:
    """The committed dummy db tarball should validate cleanly."""
    m = validate_database_layout(unpacked_dummy_db)
    assert m.id == "KS_dummy_test_v1"
    assert m.index.type == "kmc"
    assert "chromosome" in m.feature_sets
    assert "region" in m.feature_sets
