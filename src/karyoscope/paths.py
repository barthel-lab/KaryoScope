"""Filesystem path resolution for KaryoScope.

The default location for installed databases follows this precedence:

1. An explicit path passed by the caller (typically via a `--db` CLI flag).
2. The `KARYOSCOPE_DB` environment variable.
3. `~/.karyoscope/db/` (created lazily on first use).

Commands that need to find an installed database should go through
``default_db_root()`` rather than hardcoding any of the above.
"""

from __future__ import annotations

import os
from pathlib import Path

#: The environment variable users can set to override the default db root.
DB_ENV_VAR = "KARYOSCOPE_DB"


def _fallback_db_root() -> Path:
    """The fallback database root, computed lazily so $HOME changes apply."""
    return Path.home() / ".karyoscope" / "db"


def default_db_root(explicit: Path | str | None = None) -> Path:
    """Return the directory that contains installed KaryoScope databases.

    Precedence:

    1. ``explicit`` if provided (usually from a ``--db`` CLI flag).
    2. The ``KARYOSCOPE_DB`` environment variable.
    3. ``~/.karyoscope/db/`` as a last resort.

    The returned path is not guaranteed to exist; callers that intend to
    write should call :meth:`pathlib.Path.mkdir` with ``parents=True,
    exist_ok=True`` as needed.
    """
    if explicit is not None:
        return Path(explicit).expanduser().resolve()

    env_value = os.environ.get(DB_ENV_VAR)
    if env_value:
        return Path(env_value).expanduser().resolve()

    return _fallback_db_root()


def ensure_db_root(explicit: Path | str | None = None) -> Path:
    """Like :func:`default_db_root`, but creates the directory if missing."""
    root = default_db_root(explicit)
    root.mkdir(parents=True, exist_ok=True)
    return root


def installed_databases(db_root: Path | None = None) -> list[Path]:
    """Return a sorted list of directories under ``db_root`` that look like
    installed databases (i.e., contain a ``manifest.yaml`` at their top level).

    If ``db_root`` does not exist, returns an empty list.
    """
    root = db_root if db_root is not None else default_db_root()
    if not root.is_dir():
        return []
    return sorted(p for p in root.iterdir() if p.is_dir() and (p / "manifest.yaml").is_file())
