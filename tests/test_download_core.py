"""Tests for :mod:`karyoscope.download`."""

from __future__ import annotations

import tarfile
from pathlib import Path

import pytest

from karyoscope.download import install_database, is_installed
from karyoscope.exceptions import (
    ChecksumError,
    DatabaseLayoutError,
    IncompatibleVersionError,
)
from karyoscope.installed import load
from karyoscope.registry import DatabaseEntry, TaxonomyEntry


def _entry_from_dummy(url: str, sha256: str, **overrides: object) -> DatabaseEntry:
    """Build a DatabaseEntry pointing at the dummy db tarball."""
    base = dict(
        id="KS_dummy_test_v1",
        title="Dummy",
        version="1.0.0",
        karyoscope_min_version="0.1.0",
        taxonomy=[TaxonomyEntry(genus="Synthetic", species="testus")],
        feature_sets=["chromosome", "region"],
        size_gb=0.0001,
        source="synthetic",
        url=url,
        sha256=sha256,
        kmer_size=21,
        kmer_type="fixed",
        kmer_max_size=21,
    )
    base.update(overrides)
    return DatabaseEntry(**base)  # type: ignore[arg-type]


# --- Happy path ---------------------------------------------------------


def test_install_database_extracts_and_records(
    tmp_path: Path, dummy_db_url: str, dummy_db_sha256: str
) -> None:
    db_root = tmp_path / "db_root"
    entry = _entry_from_dummy(dummy_db_url, dummy_db_sha256)

    target = install_database(entry, db_root, show_progress=False)

    # Directory exists with expected files.
    assert target.is_dir()
    assert (target / "manifest.yaml").is_file()
    assert (target / "index" / "features.kmc_pre").is_file()

    # Recorded in installed.json.
    state = load(db_root)
    assert "KS_dummy_test_v1" in state.databases
    rec = state.databases["KS_dummy_test_v1"]
    assert rec.version == "1.0.0"
    assert rec.source_url == dummy_db_url


def test_is_installed_round_trip(tmp_path: Path, dummy_db_url: str, dummy_db_sha256: str) -> None:
    db_root = tmp_path / "db_root"
    entry = _entry_from_dummy(dummy_db_url, dummy_db_sha256)
    assert is_installed(db_root, entry.id) is False
    install_database(entry, db_root, show_progress=False)
    assert is_installed(db_root, entry.id) is True


def test_install_is_idempotent_without_force(
    tmp_path: Path, dummy_db_url: str, dummy_db_sha256: str
) -> None:
    db_root = tmp_path / "db_root"
    entry = _entry_from_dummy(dummy_db_url, dummy_db_sha256)
    install_database(entry, db_root, show_progress=False)
    first_mtime = (db_root / entry.id / "manifest.yaml").stat().st_mtime

    # A second install without force=True should be a no-op.
    install_database(entry, db_root, show_progress=False)
    second_mtime = (db_root / entry.id / "manifest.yaml").stat().st_mtime
    assert first_mtime == second_mtime


def test_install_force_reinstalls(tmp_path: Path, dummy_db_url: str, dummy_db_sha256: str) -> None:
    db_root = tmp_path / "db_root"
    entry = _entry_from_dummy(dummy_db_url, dummy_db_sha256)
    install_database(entry, db_root, show_progress=False)

    # Touch a file we expect to be removed by --force.
    sentinel = db_root / entry.id / "extra_file.txt"
    sentinel.write_text("should be wiped")

    install_database(entry, db_root, show_progress=False, force=True)
    assert not sentinel.exists()


# --- Failure modes ------------------------------------------------------


def test_install_rejects_bad_checksum(tmp_path: Path, dummy_db_url: str) -> None:
    db_root = tmp_path / "db_root"
    bad = "f" * 64
    entry = _entry_from_dummy(dummy_db_url, bad)
    with pytest.raises(ChecksumError):
        install_database(entry, db_root, show_progress=False)

    # Nothing should have been recorded.
    assert load(db_root).databases == {}


def test_install_skip_checksum(tmp_path: Path, dummy_db_url: str) -> None:
    db_root = tmp_path / "db_root"
    entry = _entry_from_dummy(dummy_db_url, "garbage")
    # Should succeed despite wrong checksum, since we're not verifying.
    install_database(entry, db_root, show_progress=False, verify_checksum=False)
    assert is_installed(db_root, entry.id)


def test_install_rejects_incompatible_version(
    tmp_path: Path, dummy_db_url: str, dummy_db_sha256: str
) -> None:
    db_root = tmp_path / "db_root"
    entry = _entry_from_dummy(dummy_db_url, dummy_db_sha256, karyoscope_min_version="999.0.0")
    with pytest.raises(IncompatibleVersionError, match=r"999\.0\.0"):
        install_database(entry, db_root, show_progress=False)


def test_install_rejects_path_traversal_archive(tmp_path: Path) -> None:
    """An archive with a '..' entry must be rejected before extraction."""
    archive = tmp_path / "evil.tar.gz"
    db_root = tmp_path / "db_root"
    db_root.mkdir()

    # Build an evil tarball with one good top-level entry and one '..' entry.
    payload_dir = tmp_path / "build"
    payload_dir.mkdir()
    (payload_dir / "ok.txt").write_text("safe")
    with tarfile.open(archive, "w:gz") as tar:
        # Top-level dir with the correct name (avoids the "wrong top-level" check).
        info = tarfile.TarInfo("evil_db")
        info.type = tarfile.DIRTYPE
        info.mode = 0o755
        tar.addfile(info)
        # The '..' escape attempt.
        bad_info = tarfile.TarInfo("evil_db/../escape.txt")
        bad_info.size = 5
        tar.addfile(bad_info, fileobj=__import__("io").BytesIO(b"pwned"))

    # Compute the SHA-256 so we can build a valid-looking entry.
    import hashlib

    sha = hashlib.sha256(archive.read_bytes()).hexdigest()

    entry = _entry_from_dummy(archive.absolute().as_uri(), sha, id="evil_db")

    with pytest.raises(DatabaseLayoutError, match=r"'\.\.'"):
        install_database(entry, db_root, show_progress=False)


def test_install_rejects_symlink_archive(tmp_path: Path) -> None:
    """An archive containing a symlink must be rejected."""
    archive = tmp_path / "evil.tar.gz"
    db_root = tmp_path / "db_root"
    db_root.mkdir()

    with tarfile.open(archive, "w:gz") as tar:
        info = tarfile.TarInfo("evil_db")
        info.type = tarfile.DIRTYPE
        info.mode = 0o755
        tar.addfile(info)
        link_info = tarfile.TarInfo("evil_db/sneaky_link")
        link_info.type = tarfile.SYMTYPE
        link_info.linkname = "/etc/passwd"
        tar.addfile(link_info)

    import hashlib

    sha = hashlib.sha256(archive.read_bytes()).hexdigest()

    entry = _entry_from_dummy(archive.absolute().as_uri(), sha, id="evil_db")

    with pytest.raises(DatabaseLayoutError, match="special files"):
        install_database(entry, db_root, show_progress=False)


def test_install_rejects_wrong_toplevel_dir(tmp_path: Path) -> None:
    """The single top-level entry must match the database id."""
    archive = tmp_path / "wrong.tar.gz"
    db_root = tmp_path / "db_root"
    db_root.mkdir()

    with tarfile.open(archive, "w:gz") as tar:
        info = tarfile.TarInfo("wrong_name")
        info.type = tarfile.DIRTYPE
        info.mode = 0o755
        tar.addfile(info)

    import hashlib

    sha = hashlib.sha256(archive.read_bytes()).hexdigest()
    entry = _entry_from_dummy(archive.absolute().as_uri(), sha, id="expected_id")

    with pytest.raises(DatabaseLayoutError, match="expected top-level"):
        install_database(entry, db_root, show_progress=False)
