"""Tests for :mod:`karyoscope.download`."""

from __future__ import annotations

import tarfile
from pathlib import Path

import pytest

from karyoscope import diskspace
from karyoscope.download import (
    _looks_like_sha256,
    _looks_like_url,
    _looks_like_version,
    install_database,
    is_installed,
)
from karyoscope.exceptions import (
    ChecksumError,
    DatabaseLayoutError,
    FetchError,
    IncompatibleVersionError,
    InsufficientDiskSpaceError,
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


# --- Registry hygiene checks (fire before any network or filesystem work) -


class TestLooksLikeUrl:
    """Coarse URL format check used by the registry hygiene guard."""

    @pytest.mark.parametrize(
        "url",
        [
            "https://example.com/db.tar.gz",
            "http://example.com/db.tar.gz",
            "file:///local/db.tar.gz",  # tests + cached-local-download path
        ],
    )
    def test_accepts_known_schemes(self, url: str) -> None:
        assert _looks_like_url(url)

    @pytest.mark.parametrize(
        "url",
        [
            "PLACEHOLDER",
            "",
            "example.com/db.tar.gz",  # missing scheme
            "ftp://example.com/db.tar.gz",  # not supported by fetch
            "s3://my-bucket/db.tar.gz",
            "/local/path",
        ],
    )
    def test_rejects_bad_values(self, url: str) -> None:
        assert not _looks_like_url(url)


class TestLooksLikeSha256:
    def test_accepts_lowercase_hex(self) -> None:
        assert _looks_like_sha256("a" * 64)

    def test_accepts_uppercase_hex(self) -> None:
        assert _looks_like_sha256("A" * 64)

    def test_accepts_mixed_case_hex(self) -> None:
        assert _looks_like_sha256("aB" * 32)

    @pytest.mark.parametrize(
        "s",
        [
            "PLACEHOLDER",
            "",
            "a" * 63,  # too short
            "a" * 65,  # too long
            "g" * 64,  # non-hex char
            "a" * 63 + "z",
        ],
    )
    def test_rejects_bad_values(self, s: str) -> None:
        assert not _looks_like_sha256(s)


class TestLooksLikeVersion:
    @pytest.mark.parametrize(
        "v",
        ["1.0.0", "0.1.0.dev0", "2.0.0a1", "10.20.30"],
    )
    def test_accepts_real_versions(self, v: str) -> None:
        assert _looks_like_version(v)

    @pytest.mark.parametrize(
        "v",
        ["PLACEHOLDER", "", "v1.0.0", "latest", ".1.0"],
    )
    def test_rejects_bad_values(self, v: str) -> None:
        assert not _looks_like_version(v)


def test_install_rejects_placeholder_url(tmp_path: Path) -> None:
    """A registry entry with ``PLACEHOLDER`` for ``url`` is detected
    BEFORE any network or filesystem work. No need to mock the network
    -- the validation fires up front."""
    db_root = tmp_path / "db_root"
    entry = _entry_from_dummy(url="PLACEHOLDER", sha256="a" * 64)
    with pytest.raises(FetchError, match="usable download URL"):
        install_database(entry, db_root, show_progress=False)


def test_install_rejects_placeholder_url_protects_existing_install(
    tmp_path: Path,
) -> None:
    """Critical safety property: the URL validation fires BEFORE the
    rmtree step that prepares the target directory for force-reinstall.
    A malformed registry entry must never destroy an existing install
    just because the URL turned out to be bogus.
    """
    db_root = tmp_path / "db_root"
    # Pre-create a directory the installer would normally remove.
    existing = db_root / "KS_dummy_test_v1"
    existing.mkdir(parents=True)
    sentinel = existing / "important.txt"
    sentinel.write_text("do not delete me")

    entry = _entry_from_dummy(url="PLACEHOLDER", sha256="a" * 64)
    with pytest.raises(FetchError):
        install_database(entry, db_root, show_progress=False, force=True)
    # The pre-existing directory and file are intact.
    assert sentinel.is_file()
    assert sentinel.read_text() == "do not delete me"


def test_install_rejects_placeholder_sha256(tmp_path: Path) -> None:
    db_root = tmp_path / "db_root"
    entry = _entry_from_dummy(url="https://example.com/db.tar.gz", sha256="PLACEHOLDER")
    with pytest.raises(FetchError, match="invalid SHA-256"):
        install_database(entry, db_root, show_progress=False)


def test_install_skips_sha256_check_when_verification_disabled(
    tmp_path: Path,
) -> None:
    """When --no-checksum is passed, the SHA-256 hygiene check is
    skipped (the sha256 isn't used in that mode anyway). The URL
    check still fires; we test that here by giving a placeholder URL
    so we know the URL check ran but didn't get to the sha256 step."""
    db_root = tmp_path / "db_root"
    entry = _entry_from_dummy(url="PLACEHOLDER", sha256="PLACEHOLDER")
    # URL check fires first; sha256 check is moot in --no-checksum mode.
    with pytest.raises(FetchError, match="usable download URL"):
        install_database(entry, db_root, show_progress=False, verify_checksum=False)


def test_install_rejects_placeholder_min_version(tmp_path: Path) -> None:
    """A PLACEHOLDER karyoscope_min_version silently parses to (0,)
    in :func:`_check_version_compatibility`, which always passes. The
    hygiene check catches this so a malformed entry doesn't bypass
    the compat guard entirely."""
    db_root = tmp_path / "db_root"
    entry = _entry_from_dummy(
        url="https://example.com/db.tar.gz",
        sha256="a" * 64,
        karyoscope_min_version="PLACEHOLDER",
    )
    with pytest.raises(FetchError, match="invalid karyoscope_min_version"):
        install_database(entry, db_root, show_progress=False)


def test_install_version_compat_fires_before_hygiene(tmp_path: Path) -> None:
    """When both checks would fail, the version-compat error wins --
    it's a more user-actionable message (you can upgrade KaryoScope)
    than the registry-hygiene one (which mostly means "wait")."""
    db_root = tmp_path / "db_root"
    entry = _entry_from_dummy(
        url="PLACEHOLDER",  # would fail hygiene check
        sha256="PLACEHOLDER",
        karyoscope_min_version="999.0.0",  # would fail compat check
    )
    # Compat error fires first.
    with pytest.raises(IncompatibleVersionError, match=r"999\.0\.0"):
        install_database(entry, db_root, show_progress=False)


# --- Free-space precheck ------------------------------------------------


def test_install_refuses_when_the_database_root_is_too_small(
    tmp_path: Path,
    dummy_db_url: str,
    dummy_db_sha256: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The reported bug: a 13 GB download that needs 36 GB of free space.

    The old code only learned this at the very end, when tarfile hit
    ENOSPC after a 25-minute transfer had already completed.
    """
    db_root = tmp_path / "db"
    entry = _entry_from_dummy(dummy_db_url, dummy_db_sha256, size_gb=22.7, download_size_gb=13.3)
    monkeypatch.setattr(diskspace, "free_bytes", lambda _p: 12 * 1024**3)

    with pytest.raises(InsufficientDiskSpaceError) as excinfo:
        install_database(entry, db_root, show_progress=False)
    message = str(excinfo.value)
    assert "KS_dummy_test_v1" in message
    # The peak, not either individual size, is what the user must have free.
    assert "36" in message or "37" in message
    assert not is_installed(db_root, entry.id)


def test_space_check_is_skippable(
    tmp_path: Path,
    dummy_db_url: str,
    dummy_db_sha256: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """--no-space-check exists for entries whose declared sizes are wrong."""
    db_root = tmp_path / "db"
    entry = _entry_from_dummy(dummy_db_url, dummy_db_sha256, size_gb=999.0)
    monkeypatch.setattr(diskspace, "free_bytes", lambda _p: 1024**3)

    target = install_database(entry, db_root, show_progress=False, check_space=False)
    assert target.is_dir()
    assert is_installed(db_root, entry.id)


def test_space_check_runs_before_an_existing_install_is_removed(
    tmp_path: Path,
    dummy_db_url: str,
    dummy_db_sha256: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A doomed --force reinstall must not destroy the working install first.

    install_database rmtree's the target directory before downloading, so
    ordering the check after it would leave a user with neither the old
    database nor the new one.
    """
    db_root = tmp_path / "db"
    entry = _entry_from_dummy(dummy_db_url, dummy_db_sha256)
    install_database(entry, db_root, show_progress=False)
    target = db_root / entry.id
    assert target.is_dir()

    huge = _entry_from_dummy(dummy_db_url, dummy_db_sha256, size_gb=999.0)
    monkeypatch.setattr(diskspace, "free_bytes", lambda _p: 1024**3)
    with pytest.raises(InsufficientDiskSpaceError):
        install_database(huge, db_root, show_progress=False, force=True)

    assert target.is_dir(), "the existing install was destroyed by a check that came too late"
    assert is_installed(db_root, entry.id)


def test_space_check_credits_the_install_being_replaced(
    tmp_path: Path,
    dummy_db_url: str,
    dummy_db_sha256: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Reinstalling the same database frees its own space first.

    Without the credit, `--force` on a nearly-full disk would be rejected
    even though the replacement is the same size as what it removes.
    """
    db_root = tmp_path / "db"
    entry = _entry_from_dummy(dummy_db_url, dummy_db_sha256)
    install_database(entry, db_root, show_progress=False)
    existing = diskspace.directory_size(db_root / entry.id)
    assert existing > 0

    # Free space is zero; only the credit for the outgoing install can
    # make room for the incoming one.
    big = _entry_from_dummy(
        dummy_db_url,
        dummy_db_sha256,
        size_gb=existing / diskspace.GB / 4,
        download_size_gb=existing / diskspace.GB / 4,
    )
    monkeypatch.setattr(diskspace, "free_bytes", lambda _p: 0)
    install_database(big, db_root, show_progress=False, force=True)
    assert is_installed(db_root, entry.id)


def test_enospc_during_extraction_is_reported_as_a_space_problem(
    tmp_path: Path,
    dummy_db_url: str,
    dummy_db_sha256: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A disk that fills up mid-extract must not surface as a tarfile traceback."""
    import errno

    from karyoscope import download as download_module

    db_root = tmp_path / "db"
    entry = _entry_from_dummy(dummy_db_url, dummy_db_sha256)

    def boom(*_args: object, **_kwargs: object) -> Path:
        raise OSError(errno.ENOSPC, "No space left on device", str(db_root / "x"))

    monkeypatch.setattr(download_module, "_safe_extract_tar", boom)
    with pytest.raises(InsufficientDiskSpaceError) as excinfo:
        install_database(entry, db_root, show_progress=False)
    assert "ran out of disk space" in str(excinfo.value)
