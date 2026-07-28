"""Free-disk-space accounting for commands that write large artifacts.

KaryoScope routinely writes tens of gigabytes: installing the human
database costs ~36 GB at peak (a ~13 GB archive plus its ~23 GB extracted
form, both on disk at once), and a six-feature-set ``annotate`` of a
diploid human assembly writes ~29 GB of BED. When the filesystem fills up
mid-run the failure surfaces as a bare ``OSError: [Errno 28] No space left
on device`` from deep inside a write loop — after the user has already
spent 25 minutes downloading, or 20 minutes annotating.

This module provides the two halves of a better story:

* :func:`require_free_space` — an up-front check. Called before the
  expensive step, it fails immediately with the numbers the user needs
  (required, available, shortfall, which filesystem).
* :func:`enospc_error` / :func:`reframe_enospc` — after-the-fact
  translation. Any ``ENOSPC`` that still escapes gets re-raised as an
  :class:`~karyoscope.exceptions.InsufficientDiskSpaceError` naming the
  filesystem that filled up, instead of a traceback.

Units
=====

Sizes are reported in decimal GB (10^9 bytes), matching the registry's
``size_gb`` fields and the numbers Zenodo shows. Note this differs from
``df -h`` on macOS/Linux, which reports binary GiB (2^30 bytes) — 36 GB is
34 GiB. Messages say "GB" and the discrepancy is under 8%, which is well
inside the margin these estimates carry anyway.
"""

from __future__ import annotations

import errno
import logging
import os
import shutil
from pathlib import Path

from karyoscope.exceptions import InsufficientDiskSpaceError

logger = logging.getLogger(__name__)

#: Extra headroom added to every requirement, as a fraction. Filesystems
#: behave badly when driven to exactly zero free bytes, sizes in the
#: registry are rounded, and the OS needs room for its own metadata.
DEFAULT_MARGIN = 0.05

#: Bytes per decimal GB. Everything user-facing is expressed in these.
GB = 1_000_000_000


def format_bytes(n: float) -> str:
    """Render a byte count as a human-readable decimal-unit string."""
    n = float(n)
    for unit, scale in (("TB", GB * 1000), ("GB", GB), ("MB", 1_000_000), ("kB", 1000)):
        if abs(n) >= scale:
            return f"{n / scale:.1f} {unit}"
    return f"{int(n)} B"


def _nearest_existing(path: Path) -> Path:
    """Walk up from ``path`` until we reach a directory that exists.

    ``shutil.disk_usage`` needs an existing path, but callers routinely ask
    about an output directory that hasn't been created yet. Its parent lives
    on the same filesystem in every case that matters here.
    """
    p = path.resolve() if path.is_absolute() else (Path.cwd() / path).resolve()
    while not p.exists():
        parent = p.parent
        if parent == p:  # reached the root without finding anything
            return Path(p.anchor or ".")
        p = parent
    return p


def free_bytes(path: Path) -> int:
    """Return the free bytes on the filesystem holding ``path``.

    ``path`` need not exist yet; the nearest existing ancestor is used.
    Returns ``0`` if the filesystem cannot be queried at all (which makes
    callers fail loudly rather than silently skipping the check).
    """
    target = _nearest_existing(path)
    try:
        return shutil.disk_usage(target).free
    except OSError as exc:  # pragma: no cover — platform-specific
        logger.debug("could not stat filesystem for %s: %s", target, exc)
        return 0


def directory_size(path: Path) -> int:
    """Total apparent size of every regular file under ``path``, in bytes.

    Used to credit space that a command is about to free (e.g. the existing
    database directory that ``download --force`` removes before extracting
    the replacement). Symlinks are not followed and unreadable entries are
    skipped; both would only make the estimate less conservative.
    """
    total = 0
    if not path.is_dir():
        return 0
    for root, _dirs, files in os.walk(path, followlinks=False):
        root_path = Path(root)
        for name in files:
            try:
                total += (root_path / name).lstat().st_size
            except OSError:
                continue
    return total


def require_free_space(
    path: Path,
    required_bytes: float,
    *,
    what: str,
    margin: float = DEFAULT_MARGIN,
    credit_bytes: int = 0,
    estimated: bool = False,
    hint: str | None = None,
    skip: bool = False,
) -> None:
    """Fail unless the filesystem holding ``path`` can hold ``required_bytes``.

    Parameters
    ----------
    path
        A file or directory on the filesystem that will be written to. It
        does not need to exist yet.
    required_bytes
        Peak bytes the operation needs simultaneously on that filesystem.
    what
        Short description of the operation, used in the message — e.g.
        ``"installing HKS_human_CHM13_v2"``.
    margin
        Fractional headroom added on top of ``required_bytes``.
    credit_bytes
        Bytes the operation is about to free on the same filesystem before
        it needs the space (e.g. an old install that gets removed first).
    estimated
        Whether ``required_bytes`` is a projection rather than a known
        figure. Only changes the wording, so the user knows how much to
        trust the number.
    hint
        Extra command-specific advice appended to the error.
    skip
        If True, log the numbers at INFO and return without raising. Wired
        to the ``--no-resource-check`` flags.

    Raises
    ------
    InsufficientDiskSpaceError
        If free space (plus ``credit_bytes``) is below the requirement.
    """
    needed = int(required_bytes * (1.0 + margin))
    available = free_bytes(path) + credit_bytes
    filesystem = _nearest_existing(path)

    if skip:
        logger.info(
            "space check skipped for %s: need ~%s, %s free on %s",
            what,
            format_bytes(needed),
            format_bytes(available),
            filesystem,
        )
        return

    logger.debug(
        "space check for %s: need ~%s (incl. %.0f%% margin), %s available on %s",
        what,
        format_bytes(needed),
        margin * 100,
        format_bytes(available),
        filesystem,
    )
    if available >= needed:
        return

    qualifier = "an estimated " if estimated else ""
    lines = [
        f"not enough free disk space for {what}.",
        f"  Needs:     {qualifier}{format_bytes(needed)} "
        f"(includes a {margin * 100:.0f}% safety margin)",
        f"  Available: {format_bytes(available)} on {filesystem}",
        f"  Short by:  {format_bytes(needed - available)}",
    ]
    if hint:
        lines.append("")
        lines.append(hint)
    raise InsufficientDiskSpaceError("\n".join(lines))


def enospc_error(
    exc: OSError, *, what: str, path: Path | None = None
) -> InsufficientDiskSpaceError:
    """Build a readable error for an ``ENOSPC`` that escaped mid-run.

    The originating ``OSError`` says only "No space left on device"; this
    adds what we were doing, which filesystem is full, and how much is
    left on it (usually ~0, but a quota or a reserved-blocks setting can
    make it non-zero, which is itself worth showing).
    """
    target = Path(getattr(exc, "filename", None) or path or Path.cwd())
    filesystem = _nearest_existing(target)
    return InsufficientDiskSpaceError(
        f"ran out of disk space while {what}.\n"
        f"  Filesystem: {filesystem} ({format_bytes(free_bytes(filesystem))} free)\n"
        f"  Failed on:  {target}\n"
        "Free up space (or point the command at a larger filesystem) and "
        "re-run. Partial output files from this run can be deleted."
    )


def is_enospc(exc: BaseException) -> bool:
    """True if ``exc`` is an OSError meaning "the disk is full"."""
    return isinstance(exc, OSError) and exc.errno in (errno.ENOSPC, errno.EDQUOT)


def reframe_enospc(exc: OSError, *, what: str, path: Path | None = None) -> BaseException:
    """Return the exception to raise for ``exc``, translating ENOSPC.

    Non-ENOSPC errors are returned unchanged so callers can ``raise`` the
    result unconditionally.
    """
    if is_enospc(exc):
        return enospc_error(exc, what=what, path=path)
    return exc
