"""High-level download orchestration.

This module composes the lower-level pieces (:mod:`registry`,
:mod:`_fetch`, :mod:`manifest`, :mod:`installed`) into the operations the
``karyoscope download`` CLI command exposes:

* :func:`install_database` — download an entry, verify it, extract it,
  validate the resulting layout, and record it in installed.json.
* :func:`is_installed` — quick check based on installed.json.

Anything more complex (filtering by tag, listing, removing) lives in the
CLI command itself, since it's mostly presentation logic.
"""

from __future__ import annotations

import logging
import tarfile
from pathlib import Path

from karyoscope import diskspace
from karyoscope._fetch import fetch, sha256_file
from karyoscope._version import __version__
from karyoscope.exceptions import (
    DatabaseLayoutError,
    FetchError,
    IncompatibleVersionError,
    KaryoscopeError,
)
from karyoscope.installed import InstalledRecord, load, now_iso, record_install
from karyoscope.manifest import check_install_readiness, validate_database_layout
from karyoscope.registry import DatabaseEntry
from karyoscope.versions import at_least

logger = logging.getLogger(__name__)


def is_installed(db_root: Path, db_id: str) -> bool:
    """Quick check: is this database id present in installed.json?"""
    return db_id in load(db_root).databases


def _check_version_compatibility(entry: DatabaseEntry) -> None:
    """Raise :class:`IncompatibleVersionError` if this KaryoScope is too old."""
    if not at_least(__version__, entry.karyoscope_min_version):
        raise IncompatibleVersionError(
            f"database '{entry.id}' requires KaryoScope >= "
            f"{entry.karyoscope_min_version}, but this is {__version__}. "
            "Upgrade KaryoScope to install this database."
        )


#: URL schemes the downloader knows how to handle. Anything not
#: matching these prefixes (including the literal "PLACEHOLDER" we
#: ship in incomplete registry entries) gets rejected up front with
#: an actionable message rather than an opaque urllib error.
#: Matches the schemes :func:`karyoscope._fetch.fetch` supports;
#: keep in sync if that grows new ones.
_VALID_URL_PREFIXES: tuple[str, ...] = ("http://", "https://", "file://")


def _looks_like_url(s: str) -> bool:
    return any(s.startswith(prefix) for prefix in _VALID_URL_PREFIXES)


def _looks_like_sha256(s: str) -> bool:
    """64-char lowercase or uppercase hex digest."""
    return len(s) == 64 and all(c in "0123456789abcdefABCDEF" for c in s)


def _looks_like_version(s: str) -> bool:
    """Coarse "starts with a digit" version-string check.

    Doesn't try to be strict PEP 440 -- the goal is just to catch
    ``"PLACEHOLDER"``, empty-string, and typo cases that would
    otherwise silently parse to ``(0,)`` (defeating the compat check).
    Real versions always start with a digit.
    """
    return bool(s) and s[0].isdigit()


def _check_registry_entry_publishable(entry: DatabaseEntry, *, verify_checksum: bool) -> None:
    """Validate that a registry entry has usable download metadata.

    Complementary to :func:`_check_version_compatibility`: that
    function validates the *installed* KaryoScope is new enough for
    the entry; this one validates that the entry itself has the
    fields it needs to be downloaded at all.

    Catches three "registry-side" issues with one actionable error
    each, so users see a clear "the registry entry isn't ready yet"
    message instead of an opaque urllib / HTTP failure:

    * ``url`` is missing, empty, ``"PLACEHOLDER"``, or doesn't have
      a known URL scheme.
    * ``karyoscope_min_version`` is missing, ``"PLACEHOLDER"``, or
      doesn't start with a digit (would otherwise silently parse to
      ``(0,)`` and bypass the compat check entirely).
    * ``sha256`` is missing, ``"PLACEHOLDER"``, or not a 64-hex
      digest. Skipped when the caller passed ``--no-checksum``,
      since the sha256 isn't used in that mode.

    All three checks fire *before* the destructive ``rmtree`` step
    in :func:`install_database`, so a malformed entry can never
    destroy an existing install.
    """
    if not _looks_like_url(entry.url):
        raise FetchError(
            f"database {entry.id!r} doesn't have a usable download URL "
            f"(got {entry.url!r}). The registry entry may not yet be "
            f"finalised -- check "
            f"https://github.com/barthel-lab/KaryoScope-registry for "
            f"publication status, or wait until the upload completes."
        )
    if not _looks_like_version(entry.karyoscope_min_version):
        raise FetchError(
            f"database {entry.id!r} has an invalid karyoscope_min_version "
            f"(got {entry.karyoscope_min_version!r}). The registry entry "
            f"may not yet be finalised -- check the registry repo."
        )
    if verify_checksum and not _looks_like_sha256(entry.sha256):
        raise FetchError(
            f"database {entry.id!r} has an invalid SHA-256 "
            f"(got {entry.sha256!r}). Either wait for the registry "
            f"entry to be finalised, or pass --no-checksum to skip "
            f"verification (not recommended for production use)."
        )


def _safe_extract_tar(archive: Path, dest_dir: Path, expected_top_level: str) -> Path:
    """Extract a tarball into ``dest_dir``, with safety checks.

    Ensures that:

    * All entries are under a single top-level directory named
      ``expected_top_level``.
    * No entry attempts to escape that directory via ``..`` or absolute
      paths.
    * No entry is a hardlink, symlink, device, FIFO, or other special type
      (only regular files and directories are allowed; this mirrors what a
      well-formed KaryoScope database archive ever contains).

    Returns the path of the extracted directory.
    """
    if not tarfile.is_tarfile(archive):
        raise DatabaseLayoutError(f"{archive} is not a tar archive")

    with tarfile.open(archive, "r:*") as tar:
        members = tar.getmembers()
        if not members:
            raise DatabaseLayoutError(f"{archive} is empty")

        for m in members:
            # Disallow anything other than regular files and directories.
            if not (m.isreg() or m.isdir()):
                raise DatabaseLayoutError(
                    f"refusing to extract {m.name!r}: special files (links, "
                    "devices, etc.) are not allowed in KaryoScope database archives"
                )

            # Normalize and check for traversal.
            name = m.name.lstrip("/")
            if ".." in Path(name).parts:
                raise DatabaseLayoutError(
                    f"refusing to extract {m.name!r}: contains '..' path component"
                )

            # Require single top-level directory.
            parts = Path(name).parts
            if not parts or parts[0] != expected_top_level:
                raise DatabaseLayoutError(
                    f"archive entry {m.name!r} is not under expected top-level "
                    f"directory {expected_top_level!r}"
                )

        # Python 3.12+ adds a `filter` argument that does similar checks at the
        # tarfile level; we pass `data` for defense-in-depth in addition to our
        # own checks above. (The filter would also accept the archive even
        # without our checks; we keep both because we want to fail with our
        # own clearer error messages.)
        try:
            tar.extractall(dest_dir, filter="data")  # type: ignore[arg-type]
        except TypeError:  # pragma: no cover — Python < 3.12 fallback
            tar.extractall(dest_dir)

    extracted = dest_dir / expected_top_level
    if not extracted.is_dir():
        raise DatabaseLayoutError(f"after extraction, expected directory {extracted} not found")
    return extracted


def staged_archive_path(db_root: Path, db_id: str) -> Path:
    """Where :func:`install_database` stages a database's ``.tar.gz``.

    Public so the CLI can report on a leftover archive without duplicating
    the naming convention.
    """
    return db_root / f".{db_id}.tar.gz"


def _archive_is_reusable(
    archive_path: Path, entry: DatabaseEntry, *, verify_checksum: bool
) -> bool:
    """Whether a staged archive from an earlier run can be extracted as-is.

    An install that gets as far as extraction has already paid for the
    whole download. When extraction then fails — a full disk, an
    interrupted run — throwing the archive away costs the user the entire
    transfer again (25 minutes for the human databases). A SHA-256 match
    against the registry proves the bytes are exactly what a fresh
    download would produce, so there is nothing to gain by re-fetching.

    With ``--no-checksum`` there is no way to tell a good archive from a
    truncated one, so we reuse it (consistent with what that flag means
    elsewhere) but say so — if extraction then fails on a corrupt file,
    the error points at the archive.
    """
    if not archive_path.is_file():
        return False

    size = archive_path.stat().st_size
    if not verify_checksum:
        logger.warning(
            "reusing staged archive %s (%s) without verifying it: --no-checksum was "
            "passed. Delete it and re-run if extraction fails.",
            archive_path,
            diskspace.format_bytes(size),
        )
        return True

    logger.info(
        "found a staged archive from an earlier run (%s); verifying its SHA-256 "
        "before deciding whether to re-download",
        diskspace.format_bytes(size),
    )
    actual = sha256_file(archive_path)
    if actual.lower() == entry.sha256.lower():
        logger.info("staged archive verified; skipping the download")
        return True

    logger.warning(
        "staged archive %s does not match the registry SHA-256 (expected %s, got %s); "
        "discarding it and downloading again",
        archive_path,
        entry.sha256.lower(),
        actual,
    )
    return False


def _check_install_space(entry: DatabaseEntry, db_root: Path, *, skip: bool) -> None:
    """Verify ``db_root`` can hold the archive *and* its extracted form.

    Installing peaks at ``download + installed`` bytes: :func:`fetch` writes
    the whole ``.tar.gz`` into ``db_root`` first, and it is only unlinked
    after :func:`_safe_extract_tar` returns. For the human databases that is
    ~33 GB (KMC) or ~36 GB (HKS) — well above what the archive size alone
    suggests, which is exactly the trap this check exists to close.

    Two adjustments make the number honest rather than merely conservative:

    * a partially-downloaded archive from an interrupted run is credited
      back, since :func:`fetch` resumes into it rather than re-fetching;
    * an existing install being replaced is credited back, since it is
      removed before the download starts.
    """
    archive_path = staged_archive_path(db_root, entry.id)
    partial_path = Path(str(archive_path) + ".part")
    credit = 0
    for staged in (archive_path, partial_path):
        if staged.is_file():
            credit += staged.stat().st_size
    target_dir = db_root / entry.id
    credit += diskspace.directory_size(target_dir)

    hint_lines = [
        "Free up space, or install elsewhere with:",
        "    karyoscope download --db-root /path/on/a/larger/disk " + entry.id,
        "  (set $KARYOSCOPE_DB to that path to make it the default).",
    ]
    if not entry.download_size_is_declared:
        hint_lines.append(
            "  This registry entry does not declare a download_size_gb, so the "
            "archive was assumed to be the same size as the extracted database. "
            "Refresh the registry with --refresh-registry."
        )
    hint_lines.append("  Pass --no-resource-check to attempt the install anyway.")

    diskspace.require_free_space(
        db_root,
        entry.peak_install_bytes,
        what=f"installing {entry.id}",
        credit_bytes=credit,
        estimated=not entry.download_size_is_declared,
        hint="\n".join(hint_lines),
        skip=skip,
    )


def install_database(
    entry: DatabaseEntry,
    db_root: Path,
    *,
    verify_checksum: bool = True,
    show_progress: bool = True,
    force: bool = False,
    check_space: bool = True,
) -> Path:
    """Download, verify, extract, and register a single database.

    Parameters
    ----------
    entry
        The registry entry describing the database to install.
    db_root
        Root directory under which databases live (e.g., ``~/.karyoscope/db/``).
    verify_checksum
        If True (the default), the downloaded archive's SHA-256 is checked
        against the value in the registry entry. Set to False only for
        debugging.
    show_progress
        Whether to display tqdm progress bars during download.
    force
        If True, re-download and re-extract even if the database appears to
        be already installed.
    check_space
        If True (the default), refuse to start when ``db_root``'s filesystem
        cannot hold the archive and its extracted form at once. Set to False
        to attempt the install regardless.

    Returns
    -------
    Path
        The directory where the database is installed.

    Raises
    ------
    IncompatibleVersionError
        If the database requires a newer KaryoScope than this one.
    InsufficientDiskSpaceError
        If ``check_space`` is True and there isn't room for the install, or
        if the filesystem fills up part-way through it anyway.
    ChecksumError
        If verification fails (only raised when ``verify_checksum=True``).
    DatabaseLayoutError
        If the extracted archive is malformed or fails layout validation.
    FetchError
        On network or filesystem errors during the download.
    """
    _check_version_compatibility(entry)
    # Registry-hygiene check fires BEFORE the rmtree step below, so a
    # malformed registry entry can never destroy an existing install.
    _check_registry_entry_publishable(entry, verify_checksum=verify_checksum)
    db_root.mkdir(parents=True, exist_ok=True)
    target_dir = db_root / entry.id

    if not force and target_dir.is_dir() and is_installed(db_root, entry.id):
        logger.debug("%s already installed at %s; skipping", entry.id, target_dir)
        return target_dir

    # Space check goes here: after the "already installed" short-circuit (no
    # point checking for work we won't do) but before the rmtree below, so a
    # doomed install can never destroy a working one.
    _check_install_space(entry, db_root, skip=not check_space)

    # If the directory exists but isn't recorded (or force=True), clean it up
    # so we don't merge old and new file sets.
    if target_dir.exists():
        if not target_dir.is_dir():
            raise KaryoscopeError(
                f"{target_dir} exists but is not a directory; refusing to overwrite"
            )
        # Use _safe_remove via the installed module's uninstall? No — we don't
        # want to depend on installed.json being consistent here. Just rm -rf
        # the target directory.
        import shutil

        logger.debug("removing existing directory %s", target_dir)
        shutil.rmtree(target_dir)

    # Stage the download to a temp file inside db_root.
    archive_path = staged_archive_path(db_root, entry.id)
    try:
        if _archive_is_reusable(archive_path, entry, verify_checksum=verify_checksum):
            logger.info("reusing staged archive %s", archive_path)
        else:
            logger.info("fetching %s from %s", entry.id, entry.url)
            fetch(
                entry.url,
                archive_path,
                expected_sha256=entry.sha256 if verify_checksum else None,
                show_progress=show_progress,
            )
            if verify_checksum:
                logger.debug("SHA-256 verified: %s", entry.sha256)
            else:
                logger.warning("SHA-256 verification skipped for %s", entry.id)

        logger.debug("extracting %s into %s", archive_path.name, db_root)
        _safe_extract_tar(archive_path, db_root, expected_top_level=entry.id)
    except BaseException as exc:
        # Everything from here to the `else` runs only when the install
        # failed. Two clean-up decisions, and they go in opposite
        # directions on purpose:
        #
        #   - KEEP the archive. It is the expensive part (25 minutes for
        #     the human databases) and, when checksummed, provably
        #     complete. Deleting it here is what made a failed extraction
        #     cost a second full download.
        #   - DISCARD whatever was extracted. A half-written database
        #     directory is unusable, is deleted at the top of the next
        #     attempt anyway, and — when the failure was ENOSPC — is
        #     sitting on exactly the space needed to retry.
        reclaimed = diskspace.directory_size(target_dir)
        if reclaimed:
            import shutil

            shutil.rmtree(target_dir, ignore_errors=True)
            logger.warning(
                "removed the partially-extracted %s, reclaiming %s",
                target_dir,
                diskspace.format_bytes(reclaimed),
            )
        if archive_path.is_file():
            logger.warning(
                "keeping the downloaded archive at %s (%s); re-running "
                "`karyoscope download %s` will verify and reuse it instead of "
                "downloading again. Delete it to reclaim that space.",
                archive_path,
                diskspace.format_bytes(archive_path.stat().st_size),
                entry.id,
            )
        # The up-front check uses registry-declared sizes; a filesystem that
        # was already close to full, or shared with another job, can still
        # fill up here. Report that as a disk-space problem rather than a
        # bare "[Errno 28]" traceback out of tarfile.
        if isinstance(exc, OSError):
            raise diskspace.reframe_enospc(
                exc, what=f"installing {entry.id}", path=db_root
            ) from exc
        raise
    else:
        # Only now is the archive genuinely disposable.
        archive_path.unlink(missing_ok=True)

    # Validate the extracted layout. If this fails, leave the broken directory
    # on disk so the user can inspect; just don't record the install.
    logger.debug("validating database layout at %s", target_dir)
    try:
        manifest = validate_database_layout(target_dir)
    except DatabaseLayoutError:
        raise

    # Cross-validate hierarchy <-> colors. The hierarchy parser is
    # structural; this is the semantic check that every node a
    # downstream command might see in BED output has a colour.
    # Failing here surfaces malformed community-built databases at
    # install time rather than karyotype time. We keep the extracted
    # directory on disk for inspection but refuse to register the
    # install in installed.json.
    color_issues = check_install_readiness(target_dir, manifest)
    if color_issues:
        raise KaryoscopeError(
            f"database {entry.id!r} extracted at {target_dir} but FAILED colors "
            f"validation; refusing to register the install:\n  - "
            + "\n  - ".join(color_issues)
            + "\nFix colors.tsv (or use a different database) and re-run "
            "`karyoscope download`."
        )

    # Record the install.
    logger.debug("recording install of %s in installed.json", entry.id)
    record_install(
        db_root,
        entry.id,
        InstalledRecord(
            version=entry.version,
            installed_at=now_iso(),
            source_url=entry.url,
            source_sha256=entry.sha256,
            directory=entry.id,
            registry_doi=entry.doi,
        ),
    )
    logger.info("installed %s to %s", entry.id, target_dir)
    return target_dir
