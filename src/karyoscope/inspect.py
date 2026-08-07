"""Inspecting databases and database archives without installing them.

The logic behind ``karyoscope info``: directory sizing and the streaming
archive staging that lets :func:`karyoscope.manifest.validate_database_layout`
run against a multi-gigabyte tarball without unpacking it. Lives outside
the command layer so it is importable and testable without click, and
imports nothing heavy (``karyoscope info`` must stay fast to start).

Tar members are attacker-controlled, and two callers vet them with two
deliberately different dispositions: :func:`tar_member_issue` names what
is wrong with a member, ``download``'s extractor *raises* on it, and
:func:`stage_archive` here *skips and warns* — inspection should describe
a bad archive, not refuse to look at it.
"""

from __future__ import annotations

import contextlib
import logging
import tarfile
from pathlib import Path, PurePosixPath

logger = logging.getLogger(__name__)

#: Members whose full contents are staged by :func:`stage_archive` (the
#: text files layout validation actually reads); everything else becomes
#: an empty placeholder.
INLINE_SUFFIXES = (".yaml", ".yml", ".tsv", ".txt")
MAX_INLINE_BYTES = 128 * 1024 * 1024


def tar_member_issue(member: tarfile.TarInfo) -> str | None:
    """Why ``member`` must not be materialised, or ``None`` if it may be.

    The shared classification for KaryoScope database archives: only
    regular files and directories are ever legitimate, and no path may
    climb out of the extraction root. The caller decides what to do
    about a flagged member (raise on install, skip-and-warn on inspect).
    """
    if not (member.isreg() or member.isdir()):
        return (
            "special files (links, devices, etc.) are not allowed in KaryoScope database archives"
        )
    if ".." in PurePosixPath(member.name.lstrip("/")).parts:
        return "contains '..' path component"
    return None


def safe_relpath(name: str) -> PurePosixPath | None:
    """Return ``name`` as a relative path, or None if it escapes the root.

    Absolute paths and ``..`` components would let an archive write
    outside the staging directory; ``.`` components are dropped.
    """
    pure = PurePosixPath(name)
    if pure.is_absolute():
        return None
    parts = [p for p in pure.parts if p not in (".", "")]
    if any(p == ".." for p in parts):
        return None
    return PurePosixPath(*parts) if parts else None


def dir_size(path: Path) -> int:
    """Recursively compute the total size of files under ``path``."""
    total = 0
    for f in path.rglob("*"):
        if f.is_file():
            with contextlib.suppress(OSError):
                total += f.stat().st_size
    return total


def stage_archive(path: Path, staging: Path) -> tuple[set[str], int, int]:
    """Stream ``path`` once, reproducing its shape under ``staging``.

    Small text members are written out in full; every other regular file
    becomes an empty placeholder, so ``validate_database_layout`` can
    run its existence checks without unpacking gigabytes of index. Returns
    the set of top-level names, the number of regular files, and their
    total uncompressed size.

    Unsafe or special members are skipped with a warning rather than
    raised on — see the module docstring for why inspection fails soft.
    """
    top_level: set[str] = set()
    n_files = 0
    total_bytes = 0

    # Stream mode ("r|*"): one forward pass, no seeking back into the
    # compressed stream, so a multi-GB archive is read but never written.
    with tarfile.open(path, mode="r|*") as tar:
        for member in tar:
            rel = safe_relpath(member.name)
            if rel is None:
                logger.warning("skipping archive member with unsafe path: %s", member.name)
                continue
            issue = tar_member_issue(member)
            if issue is not None:
                logger.warning("skipping archive member %s: %s", member.name, issue)
                continue
            top_level.add(rel.parts[0])
            if member.isdir():
                (staging / rel).mkdir(parents=True, exist_ok=True)
                continue

            n_files += 1
            total_bytes += member.size
            dest = staging / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            inline = rel.name.lower().endswith(INLINE_SUFFIXES) and member.size <= MAX_INLINE_BYTES
            if inline:
                src = tar.extractfile(member)
                if src is None:  # pragma: no cover — isfile() said otherwise
                    dest.touch()
                    continue
                with src, dest.open("wb") as out:
                    while chunk := src.read(1 << 20):
                        out.write(chunk)
            else:
                dest.touch()

    return top_level, n_files, total_bytes
