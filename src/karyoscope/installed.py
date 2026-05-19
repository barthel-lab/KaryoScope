"""Tracking of locally installed databases via ``installed.json``.

The file lives at ``<db_root>/installed.json`` and records what's installed,
when it was installed, where it came from, and its expected SHA-256. This
enables ``karyoscope download --status``, ``--remove``, integrity checks,
and detection of "already installed, skip" cases.

The on-disk format is JSON (not YAML) because it is machine-generated and
machine-consumed; users are not expected to hand-edit it.
"""

from __future__ import annotations

import contextlib
import json
import os
import shutil
import tempfile
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from karyoscope.exceptions import KaryoscopeError

INSTALLED_FILENAME = "installed.json"
SCHEMA_VERSION = 1


@dataclass
class InstalledRecord:
    """One row of installed.json — the install state of a single database."""

    version: str
    installed_at: str  # ISO-8601 UTC timestamp
    source_url: str
    source_sha256: str
    directory: str  # relative to db_root
    registry_doi: str | None = None


@dataclass
class InstalledState:
    """The parsed contents of installed.json."""

    schema_version: int = SCHEMA_VERSION
    databases: dict[str, InstalledRecord] = field(default_factory=dict)


def _path(db_root: Path) -> Path:
    return db_root / INSTALLED_FILENAME


def now_iso() -> str:
    """Return an ISO-8601 timestamp in UTC, suitable for ``installed_at``."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load(db_root: Path) -> InstalledState:
    """Load installed.json, returning an empty state if the file does not exist.

    Tolerates missing files and corrupt files (the latter by treating them as
    empty, with the original preserved as ``installed.json.corrupt`` so a
    user can recover anything manually if they need to). This avoids the
    situation where one bad write blocks all future operations.
    """
    path = _path(db_root)
    if not path.is_file():
        return InstalledState()

    try:
        raw = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        backup = path.with_suffix(path.suffix + ".corrupt")
        with contextlib.suppress(OSError):  # pragma: no cover
            path.rename(backup)
        return InstalledState()

    if not isinstance(raw, dict) or raw.get("schema_version") != SCHEMA_VERSION:
        return InstalledState()

    state = InstalledState()
    for db_id, rec_raw in (raw.get("databases") or {}).items():
        if not isinstance(rec_raw, dict):
            continue
        try:
            state.databases[db_id] = InstalledRecord(
                version=rec_raw["version"],
                installed_at=rec_raw["installed_at"],
                source_url=rec_raw["source_url"],
                source_sha256=rec_raw["source_sha256"],
                directory=rec_raw["directory"],
                registry_doi=rec_raw.get("registry_doi"),
            )
        except KeyError:
            # Skip malformed entries silently rather than failing the whole load.
            continue
    return state


def save(db_root: Path, state: InstalledState) -> None:
    """Write installed.json atomically.

    Writes to a temp file in ``db_root`` and renames into place, so a crash
    mid-write does not corrupt the existing file.
    """
    db_root.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": state.schema_version,
        "databases": {db_id: asdict(rec) for db_id, rec in state.databases.items()},
    }
    # Use a tempfile in the same directory so rename is atomic on the same filesystem.
    fd, tmp_name = tempfile.mkstemp(dir=db_root, prefix=".installed_", suffix=".tmp")
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(payload, f, indent=2, sort_keys=True)
            f.write("\n")
        tmp_path.replace(_path(db_root))
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def record_install(db_root: Path, db_id: str, record: InstalledRecord) -> None:
    """Add (or replace) an entry in installed.json."""
    state = load(db_root)
    state.databases[db_id] = record
    save(db_root, state)


def remove_install_record(db_root: Path, db_id: str) -> bool:
    """Remove a database's entry from installed.json. Returns True if removed."""
    state = load(db_root)
    if db_id not in state.databases:
        return False
    del state.databases[db_id]
    save(db_root, state)
    return True


def uninstall(db_root: Path, db_id: str) -> bool:
    """Remove a database's files and its installed.json entry.

    Returns True if the database was found and removed, False if it was not
    installed. Refuses to remove anything outside of ``db_root`` as a basic
    safety measure.

    Raises
    ------
    KaryoscopeError
        If the recorded directory escapes ``db_root`` (which would indicate
        a tampered or corrupted installed.json).
    """
    state = load(db_root)
    record = state.databases.get(db_id)
    if record is None:
        return False

    target = (db_root / record.directory).resolve()
    db_root_resolved = db_root.resolve()
    try:
        target.relative_to(db_root_resolved)
    except ValueError as e:
        raise KaryoscopeError(
            f"refusing to remove '{target}': installed.json entry for "
            f"'{db_id}' points outside the database root"
        ) from e

    if target.is_dir():
        shutil.rmtree(target)

    remove_install_record(db_root, db_id)
    return True
