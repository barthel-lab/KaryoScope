"""``karyoscope register`` — register a manually-placed database.

``karyoscope download`` is the normal way to install a database: it
fetches a tarball, extracts it under the database root, validates it, and
records it in ``installed.json`` so the data commands can find it. But
databases are sometimes produced locally — built by hand or copied from
another machine — and unpacked into the database root directly. Such a
database is valid on disk yet invisible to ``annotate``, ``bin``,
``scaffold``, ``centromeres``, and ``karyotype``, because those commands
resolve databases through ``installed.json`` only (see
:func:`karyoscope.core.annotate.resolve_database`).

``register`` closes that gap: point it at an already-present database
directory and it writes the ``installed.json`` entry, deriving everything
it can from the manifest. The ``version`` and ``id`` come from
``manifest.yaml``; ``installed_at`` is the time of registration (matching
how ``download`` records it); ``source_url`` is recorded as ``local`` and
there is no archive to checksum, so ``source_sha256`` is left empty.
"""

from __future__ import annotations

import logging
from pathlib import Path

import click

from karyoscope import installed as _installed
from karyoscope import paths
from karyoscope.exceptions import (
    DatabaseLayoutError,
    KaryoscopeError,
    ManifestError,
)
from karyoscope.installed import InstalledRecord, now_iso
from karyoscope.manifest import validate_database_layout

logger = logging.getLogger(__name__)


def _resolve_target_dir(db_root: Path, target: str) -> Path:
    """Resolve TARGET (a database id or a filesystem path) to a directory.

    An existing path is used as-is; otherwise TARGET is treated as a
    database id located directly under ``db_root``.
    """
    as_path = Path(target).expanduser()
    if as_path.exists():
        return as_path.resolve()
    return (db_root / target).resolve()


def _register(db_root: Path, target: str, *, force: bool) -> None:
    db_root = db_root.resolve()
    db_dir = _resolve_target_dir(db_root, target)

    if not db_dir.is_dir():
        raise DatabaseLayoutError(
            f"no database directory found for {target!r} (looked at {db_dir}). "
            f"Place the database under {db_root} first, or pass a path to it."
        )

    # installed.json stores the directory as a path relative to db_root, and
    # the data commands resolve it as `db_root / directory`. So the database
    # must live inside the root.
    try:
        rel = db_dir.relative_to(db_root)
    except ValueError as e:
        raise DatabaseLayoutError(
            f"database directory {db_dir} is not inside the database root {db_root}. "
            "Move it under the root (or pass --db-root pointing at its parent) and retry."
        ) from e

    # Validates the manifest and that the referenced files exist; returns the
    # parsed manifest. This is the same check `download` runs at install time.
    manifest = validate_database_layout(db_dir)
    db_id = manifest.id

    if db_dir.name != db_id:
        logger.warning(
            "directory name %r does not match manifest id %r; the database "
            "layout spec expects them to match. Recording under id %r.",
            db_dir.name,
            db_id,
            db_id,
        )

    state = _installed.load(db_root)
    already = db_id in state.databases
    if already and not force:
        existing = state.databases[db_id]
        raise KaryoscopeError(
            f"database {db_id!r} is already registered (directory "
            f"{existing.directory!r}). Re-run with --force to overwrite the entry."
        )

    installed_at = now_iso()
    _installed.record_install(
        db_root,
        db_id,
        InstalledRecord(
            version=manifest.version,
            installed_at=installed_at,
            source_url="local",
            source_sha256="",
            directory=str(rel),
            registry_doi=None,
        ),
    )

    click.echo(f"{'Re-registered' if already else 'Registered'} {db_id} v{manifest.version}")
    click.echo(f"  Directory:    {db_dir}")
    click.echo(f"  Installed at: {installed_at}")
    click.echo("  Source:       local")
    click.echo(f"Run `karyoscope info {db_id}` to inspect it.")


@click.command(
    help="Register a database already present under the database root so the "
    "data commands can use it.",
)
@click.argument("target", metavar="DATABASE_ID_OR_PATH")
@click.option(
    "--db-root",
    "db_root_arg",
    type=click.Path(file_okay=False, path_type=Path),
    help="Override the database root directory (default: $KARYOSCOPE_DB or ~/.karyoscope/db/).",
)
@click.option(
    "--force",
    is_flag=True,
    help="Overwrite an existing installed.json entry for this database.",
)
def cmd(target: str, db_root_arg: Path | None, force: bool) -> None:
    """Record an on-disk database in installed.json so commands can use it.

    \b
    Examples:
        karyoscope register KS_human_CHM13_cytoband
        karyoscope register ./KS_human_CHM13_cytoband --db-root /data/ksdb
    """
    db_root = paths.default_db_root(db_root_arg)
    try:
        _register(db_root, target, force=force)
    except (DatabaseLayoutError, ManifestError, KaryoscopeError) as e:
        raise click.ClickException(str(e)) from e
