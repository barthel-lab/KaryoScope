"""Tests for :mod:`karyoscope.manifest`."""

from __future__ import annotations

import re
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


#: Feature sets used by the multi-set HKS dummy. Deliberately more than one so
#: the layout checks are exercised for *every* declared feature set, not just
#: the first — mirroring the real HKS database, which has one .hksf and one
#: .hierarchy.txt per feature set.
_HKS_FEATURE_SETS: tuple[str, ...] = ("chromosome", "region", "repeat")

#: Per-feature-set file suffixes the HKS layout requires (relative to basename).
_HKS_PER_FS_SUFFIXES: tuple[str, ...] = (".hksf", ".hierarchy.txt")


def _write_minimal_hks_db(dir_: Path, feature_sets: tuple[str, ...] = _HKS_FEATURE_SETS) -> Path:
    """Build a minimal, valid HKS-backed KaryoScope database layout. Returns dir_.

    Writes one ``features.<fs>.hksf`` and one ``features.<fs>.hierarchy.txt``
    per declared feature set, plus the shared ``features.hksb`` base index, so
    the manifest-driven per-feature-set layout checks have something to verify.
    """
    dir_.mkdir(parents=True, exist_ok=True)
    fs_yaml = "".join(f"  - {fs}\n" for fs in feature_sets)
    (dir_ / "manifest.yaml").write_text(
        "id: test_hks_db\n"
        "version: 1.0.0\n"
        "karyoscope_min_version: 0.1.0\n"
        "index:\n"
        "  type: hks\n"
        "  basename: index/features\n"
        "hierarchy: hierarchy.tsv\n"
        "features: features.tsv\n"
        "colors: colors.txt\n"
        "kmer:\n"
        "  size: 31\n"
        "  type: fixed\n"
        "  max_size: 31\n"
        "feature_sets:\n" + fs_yaml
    )
    (dir_ / "hierarchy.tsv").write_text("feature_set\tfeature\tparent\n")
    (dir_ / "features.tsv").write_text("feature_set\tfeature\tfeature_id\n")
    (dir_ / "colors.txt").write_text("feature_set\tfeature\tcolor\n")
    (dir_ / "index").mkdir(exist_ok=True)
    (dir_ / "index" / "features.hksb").write_bytes(b"\x00" * 8)
    for fs in feature_sets:
        (dir_ / "index" / f"features.{fs}.hksf").write_bytes(b"\x00" * 8)
        (dir_ / "index" / f"features.{fs}.hierarchy.txt").write_text("child\tcategorized\n")
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


def test_parse_manifest_hks(tmp_path: Path) -> None:
    db = _write_minimal_hks_db(tmp_path / "db")
    m = parse_manifest(db / "manifest.yaml")

    assert m.id == "test_hks_db"
    assert m.index.type == "hks"
    assert m.index.basename == "index/features"
    assert m.kmer.size == 31
    assert m.feature_sets == list(_HKS_FEATURE_SETS)


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


def test_validate_database_layout_succeeds_on_minimal_hks(tmp_path: Path) -> None:
    db = _write_minimal_hks_db(tmp_path / "db")
    m = validate_database_layout(db)
    assert m.id == "test_hks_db"
    assert m.index.type == "hks"
    assert set(m.feature_sets) == set(_HKS_FEATURE_SETS)


def test_validate_database_layout_missing_hksb(tmp_path: Path) -> None:
    db = _write_minimal_hks_db(tmp_path / "db")
    (db / "index" / "features.hksb").unlink()
    with pytest.raises(DatabaseLayoutError, match=r"\.hksb"):
        validate_database_layout(db)


@pytest.mark.parametrize("feature_set", _HKS_FEATURE_SETS)
@pytest.mark.parametrize("suffix", _HKS_PER_FS_SUFFIXES)
def test_validate_database_layout_missing_hks_per_fs_file(
    tmp_path: Path, feature_set: str, suffix: str
) -> None:
    """Every declared feature set's .hksf and .hierarchy.txt must be required.

    Manifest-driven: the validator reads ``feature_sets`` from the manifest and
    checks the matching per-feature-set files exist. Removing any one of them —
    for any feature set, not just the first — must fail validation.
    """
    db = _write_minimal_hks_db(tmp_path / "db")
    (db / "index" / f"features.{feature_set}{suffix}").unlink()
    with pytest.raises(DatabaseLayoutError, match=re.escape(suffix)):
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


def _write_hks_db_optfeatures(dir_: Path, *, with_features: bool) -> Path:
    """A minimal valid HKS database; features.tsv is optional for hks."""
    dir_.mkdir(parents=True, exist_ok=True)
    features_line = "features: features.tsv\n" if with_features else ""
    (dir_ / "manifest.yaml").write_text(
        "id: hks_db\n"
        "version: 1.0.0\n"
        "karyoscope_min_version: 1.0.0\n"
        "index:\n"
        "  type: hks\n"
        "  basename: index/features\n"
        "hierarchy: hierarchy.tsv\n"
        f"{features_line}"
        "colors: colors.tsv\n"
        "kmer:\n"
        "  size: 31\n"
        "  type: fixed\n"
        "  max_size: 31\n"
        "feature_sets: [repeat]\n"
    )
    (dir_ / "hierarchy.tsv").write_text("feature_set\tchild\tparent\n")
    (dir_ / "colors.tsv").write_text("feature_set\tfeature\tcolor\n")
    if with_features:
        (dir_ / "features.tsv").write_text("featureID\trepeat\n")
    (dir_ / "index").mkdir(exist_ok=True)
    (dir_ / "index" / "features.hksb").write_bytes(b"\x00" * 8)
    (dir_ / "index" / "features.repeat.hksf").write_bytes(b"\x00" * 8)
    (dir_ / "index" / "features.repeat.hierarchy.txt").write_text("")
    return dir_


def test_hks_manifest_without_features_is_valid(tmp_path: Path) -> None:
    db = _write_hks_db_optfeatures(tmp_path / "db", with_features=False)
    manifest = validate_database_layout(db)
    assert manifest.features is None
    assert manifest.index.type == "hks"


def test_hks_manifest_with_features_still_accepted(tmp_path: Path) -> None:
    db = _write_hks_db_optfeatures(tmp_path / "db", with_features=True)
    manifest = validate_database_layout(db)
    assert manifest.features == "features.tsv"


def test_kmc_manifest_still_requires_features(tmp_path: Path) -> None:
    db = tmp_path / "db"
    db.mkdir()
    (db / "manifest.yaml").write_text(
        "id: kmc_db\n"
        "version: 1.0.0\n"
        "karyoscope_min_version: 1.0.0\n"
        "index:\n"
        "  type: kmc\n"
        "  basename: index/features\n"
        "hierarchy: hierarchy.tsv\n"
        "colors: colors.tsv\n"
        "kmer: { size: 21, type: fixed, max_size: 21 }\n"
        "feature_sets: [chromosome]\n"
    )
    with pytest.raises(ManifestError, match="features"):
        parse_manifest(db / "manifest.yaml")
