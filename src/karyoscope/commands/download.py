"""``karyoscope download`` — acquire and manage pre-built databases.

This is the user-facing CLI command. The heavy lifting lives in
:mod:`karyoscope.registry`, :mod:`karyoscope.download`, and
:mod:`karyoscope.installed`; this module is presentation.

Supported actions (mutually exclusive):

* default (no action flag): install one or more databases
* ``--list``: print available databases
* ``--info ID``: print detailed info about a database
* ``--status``: print locally installed databases
* ``--remove ID``: uninstall a database
"""

from __future__ import annotations

from pathlib import Path

import click

from karyoscope import diskspace, paths
from karyoscope import download as _download
from karyoscope import installed as _installed
from karyoscope import registry as _registry
from karyoscope.commands._options import (
    resolve_db_root_flag,
    resolve_resource_check_flag,
)
from karyoscope.exceptions import (
    DatabaseNotFoundError,
    KaryoscopeError,
)

# --- Formatting helpers -----------------------------------------------------


def _format_size_gb(gb: float) -> str:
    """Render a size in GB with one decimal.

    Deliberately not :func:`karyoscope.diskspace.format_bytes`: registry
    entries declare their sizes in GB, and the listing keeps every row in
    that one unit so entries can be compared down a column (format_bytes
    would render a small database in MB).
    """
    return f"{gb:.1f} GB"


def _format_bytes_gb(n: float) -> str:
    """Render a byte count in GB with one decimal."""
    return _format_size_gb(n / diskspace.GB)


def _format_sizes(entry: _registry.DatabaseEntry) -> str:
    """Render an entry's download and on-disk sizes.

    These differ enough to matter — the HKS human database is a 13 GB
    download that unpacks to 23 GB — so both are shown rather than the one
    unlabelled figure we used to print. Entries that predate
    ``download_size_gb`` only have the on-disk number to show.
    """
    if not entry.download_size_is_declared:
        return f"{_format_bytes_gb(entry.installed_size_bytes)} on disk"
    return (
        f"{_format_bytes_gb(entry.download_size_bytes)} download, "
        f"{_format_bytes_gb(entry.installed_size_bytes)} on disk"
    )


def _format_db_short(entry: _registry.DatabaseEntry) -> str:
    """One-line summary used by ``--list``."""
    default_tag = " (default)" if entry.is_default else ""
    community_tag = " [community]" if entry.source_category == "community" else ""
    taxa = ", ".join(t.common_name or t.species for t in entry.taxonomy)
    return (
        f"{entry.id}{default_tag}{community_tag}  "
        f"v{entry.version}  "
        f"{_format_sizes(entry)}  "
        f"({taxa})"
    )


def _filter_databases(
    entries: list[_registry.DatabaseEntry],
    organism: str | None,
    tag: str | None,
    include_community: bool,
) -> list[_registry.DatabaseEntry]:
    """Apply --organism / --tag / --community filters to a list of entries."""
    out: list[_registry.DatabaseEntry] = []
    for e in entries:
        if not include_community and e.source_category == "community":
            continue
        if organism is not None:
            needle = organism.lower()
            haystack = []
            for t in e.taxonomy:
                if t.common_name:
                    haystack.append(t.common_name.lower())
                haystack.append(t.species.lower())
                haystack.append(t.genus.lower())
            if not any(needle in h for h in haystack):
                continue
        if tag is not None and tag not in e.tags:
            continue
        out.append(e)
    return out


# --- Sub-actions ------------------------------------------------------------


def _action_list(
    db_root: Path,
    registry_url: str,
    organism: str | None,
    tag: str | None,
    include_community: bool,
    refresh: bool,
) -> None:
    registry = _registry.load_registry(db_root, registry_url, refresh=refresh)
    entries = _filter_databases(registry.databases, organism, tag, include_community)
    if not entries:
        click.echo("No databases match the given filters.")
        return
    click.echo(f"Found {len(entries)} database(s):")
    for entry in entries:
        click.echo(f"  {_format_db_short(entry)}")
    click.echo("\nUse `karyoscope download --info <ID>` to see details for a specific database.")


def _action_info(
    db_root: Path,
    registry_url: str,
    db_id: str,
    refresh: bool,
) -> None:
    registry = _registry.load_registry(db_root, registry_url, refresh=refresh)
    entry = registry.find(db_id)
    if entry is None:
        raise DatabaseNotFoundError(f"database '{db_id}' is not in the registry")

    click.echo(f"{entry.id}")
    click.echo(f"  Title: {entry.title}")
    if entry.description:
        click.echo(f"  Description: {entry.description}")
    click.echo(f"  Version: {entry.version}")
    if entry.release_date:
        click.echo(f"  Release date: {entry.release_date}")
    click.echo(f"  KaryoScope min version: {entry.karyoscope_min_version}")
    click.echo(f"  Size on disk after install: {_format_bytes_gb(entry.installed_size_bytes)}")
    if entry.download_size_is_declared:
        click.echo(f"  Download size: {_format_bytes_gb(entry.download_size_bytes)} (.tar.gz)")
        click.echo(
            f"  Free space needed to install: "
            f"{_format_bytes_gb(entry.peak_install_bytes)} "
            "(archive and extracted database coexist until extraction finishes)"
        )
    else:
        click.echo("  Download size: not declared in the registry")
    taxa_str = "; ".join(
        f"{t.genus} {t.species}" + (f" ({t.common_name})" if t.common_name else "")
        for t in entry.taxonomy
    )
    click.echo(f"  Taxonomy: {taxa_str}")
    click.echo(
        f"  k-mer: size={entry.kmer_size}, type={entry.kmer_type}, max={entry.kmer_max_size}"
    )
    click.echo(f"  Feature sets: {', '.join(entry.feature_sets)}")
    if entry.roles:
        roles_str = ", ".join(f"{k}={v}" for k, v in entry.roles.items())
        click.echo(f"  Roles: {roles_str}")
    if entry.tags:
        click.echo(f"  Tags: {', '.join(entry.tags)}")
    click.echo(f"  Source: {entry.source}")
    click.echo(f"  URL: {entry.url}")
    click.echo(f"  SHA-256: {entry.sha256}")
    if entry.doi:
        click.echo(f"  DOI: {entry.doi}")
    if entry.citation:
        click.echo(f"  Citation: {entry.citation}")
    if entry.maintained_by:
        click.echo(f"  Maintained by: {entry.maintained_by}")
    click.echo(f"  Category: {entry.source_category}")
    if entry.is_default:
        click.echo(
            "  This is the default database (installed when `karyoscope download` "
            "is run with no arguments)."
        )


def _action_status(db_root: Path) -> None:
    state = _installed.load(db_root)
    click.echo(f"Database root: {db_root}")
    if not state.databases:
        click.echo("No databases installed.")
        return
    click.echo(f"\nInstalled databases ({len(state.databases)}):")
    for db_id, rec in sorted(state.databases.items()):
        click.echo(f"  {db_id}")
        click.echo(f"    Version: {rec.version}")
        click.echo(f"    Installed: {rec.installed_at}")
        click.echo(f"    Source: {rec.source_url}")


def _action_remove(db_root: Path, db_id: str, yes: bool) -> None:
    state = _installed.load(db_root)
    record = state.databases.get(db_id)
    if record is None:
        raise DatabaseNotFoundError(f"database '{db_id}' is not installed")

    target = db_root / record.directory
    if not yes:
        click.echo(f"About to remove database '{db_id}' from {target}.")
        click.confirm("Continue?", abort=True)

    if _installed.uninstall(db_root, db_id):
        click.echo(f"Removed {db_id}.")
    else:  # pragma: no cover — guarded by the earlier check
        raise DatabaseNotFoundError(f"database '{db_id}' was not installed")


def _action_install(
    db_root: Path,
    registry_url: str,
    db_ids: tuple[str, ...],
    *,
    refresh: bool,
    force: bool,
    verify_checksum: bool,
    show_progress: bool,
    check_space: bool,
) -> None:
    registry = _registry.load_registry(db_root, registry_url, refresh=refresh)

    # Resolve "no args" to the default database.
    if not db_ids:
        default = registry.default_database()
        if default is None:
            raise DatabaseNotFoundError(
                "no database id was given and the registry does not mark any "
                "database as default. Run `karyoscope download --list` to see "
                "available databases, then `karyoscope download <ID>` to install."
            )
        db_ids = (default.id,)

    # Resolve all ids before doing any downloads, so we fail fast on typos.
    entries: list[_registry.DatabaseEntry] = []
    for db_id in db_ids:
        entry = registry.find(db_id)
        if entry is None:
            raise DatabaseNotFoundError(f"database '{db_id}' is not in the registry")
        entries.append(entry)

    for entry in entries:
        if not force and _download.is_installed(db_root, entry.id):
            click.echo(f"{entry.id}: already installed, skipping (use --force to reinstall)")
            continue

        click.echo(f"Installing {entry.id} v{entry.version} ({_format_sizes(entry)})...")
        # Spelled out before the progress bar starts, because the peak is the
        # number that matters and it is larger than either figure above: the
        # archive is not deleted until extraction succeeds.
        click.echo(
            f"  Needs {_format_bytes_gb(entry.peak_install_bytes)} free during install; "
            f"{diskspace.format_bytes(diskspace.free_bytes(db_root))} available "
            f"on {db_root}"
        )
        # An archive left by a failed run is reused if it verifies. Say so
        # here rather than only in the log: without it, the pause while we
        # hash 13 GB looks like a hang, and the absence of a progress bar
        # looks like nothing is happening. The verdict itself is logged by
        # install_database, which is where the hashing happens.
        staged = _download.staged_archive_path(db_root, entry.id)
        if staged.is_file():
            click.echo(
                f"  Found an archive from an earlier run "
                f"({diskspace.format_bytes(staged.stat().st_size)}); "
                "verifying it — if it matches, the download is skipped."
            )
        target = _download.install_database(
            entry,
            db_root,
            verify_checksum=verify_checksum,
            show_progress=show_progress,
            force=force,
            check_space=check_space,
        )
        click.echo(f"  Installed to {target}")


# --- Click command definition -----------------------------------------------


@click.command(
    help="Download and manage pre-built KaryoScope databases.",
    no_args_is_help=False,
)
@click.argument("database_ids", nargs=-1, metavar="[DATABASE_ID ...]")
@click.option(
    "--list",
    "list_mode",
    is_flag=True,
    help="List databases available in the registry, then exit.",
)
@click.option(
    "--info",
    "info_id",
    metavar="DATABASE_ID",
    help="Show detailed information about a single database, then exit.",
)
@click.option(
    "--status",
    is_flag=True,
    help="Show locally installed databases, then exit.",
)
@click.option(
    "--remove",
    "remove_id",
    metavar="DATABASE_ID",
    help="Uninstall a locally installed database, then exit.",
)
@click.option(
    "--organism",
    metavar="NAME",
    help="Filter --list by common name / genus / species (case-insensitive substring).",
)
@click.option(
    "--tag",
    metavar="TAG",
    help="Filter --list by tag (e.g. 'test', 'reference').",
)
@click.option(
    "--community",
    is_flag=True,
    help="Include community-contributed databases in listings.",
)
@click.option(
    "--db-root",
    "db_root_arg",
    type=click.Path(file_okay=False, path_type=Path),
    help="Override the database root directory (default: $KARYOSCOPE_DB or ~/.karyoscope/db/).",
)
@click.option(
    "--db",
    "db_alias",
    type=click.Path(file_okay=False, path_type=Path),
    hidden=True,
    help="Deprecated alias for --db-root.",
)
@click.option(
    "--registry-url",
    metavar="URL",
    default=None,
    help="Override the registry URL (advanced; for testing or private registries).",
)
@click.option(
    "--refresh-registry",
    is_flag=True,
    help="Force a fresh fetch of the registry, ignoring any cached copy.",
)
@click.option(
    "--force",
    is_flag=True,
    help="Re-download and re-install even if the database is already present.",
)
@click.option(
    "--no-checksum",
    is_flag=True,
    help="Skip SHA-256 verification (not recommended; useful for debugging).",
)
@click.option(
    "--no-resource-check",
    is_flag=True,
    help="Install even if the database root looks too small to hold the archive "
    "and its extracted contents. Only useful when the registry's declared "
    "sizes are wrong for your copy of the database. (Named for consistency with "
    "the other commands; `download` only ever checks disk.)",
)
@click.option(
    "--no-space-check",
    is_flag=True,
    hidden=True,
    help="Deprecated alias for --no-resource-check.",
)
@click.option(
    "-y",
    "--yes",
    is_flag=True,
    help="Assume 'yes' to interactive prompts (e.g. --remove).",
)
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    help="Suppress progress bars.",
)
def cmd(
    database_ids: tuple[str, ...],
    list_mode: bool,
    info_id: str | None,
    status: bool,
    remove_id: str | None,
    organism: str | None,
    tag: str | None,
    community: bool,
    db_root_arg: Path | None,
    db_alias: Path | None,
    registry_url: str | None,
    refresh_registry: bool,
    force: bool,
    no_checksum: bool,
    no_resource_check: bool,
    no_space_check: bool,
    yes: bool,
    quiet: bool,
) -> None:
    """Download a database, list what's available, inspect installations.

    \b
    Common uses:
        # Download the default database
        karyoscope download

        # Download a specific database by ID
        karyoscope download KS_human_CHM13_v2

        # List what's available
        karyoscope download --list
        karyoscope download --list --organism human

        # See what's installed locally
        karyoscope download --status

        # Remove an installed database
        karyoscope download --remove KS_mouse_v1
    """
    skip_checks = resolve_resource_check_flag(no_resource_check, no_space_check, command="download")
    db_root = paths.ensure_db_root(resolve_db_root_flag(db_root_arg, db_alias, command="download"))
    effective_registry_url = registry_url or _registry.DEFAULT_REGISTRY_URL

    # The action flags are mutually exclusive at the conceptual level.
    chosen = [flag for flag in (list_mode, bool(info_id), status, bool(remove_id)) if flag]
    if len(chosen) > 1:
        raise click.UsageError(
            "Only one of --list, --info, --status, --remove may be given at a time."
        )

    try:
        if list_mode:
            _action_list(
                db_root,
                effective_registry_url,
                organism,
                tag,
                include_community=community,
                refresh=refresh_registry,
            )
            return

        if info_id is not None:
            _action_info(db_root, effective_registry_url, info_id, refresh=refresh_registry)
            return

        if status:
            _action_status(db_root)
            return

        if remove_id is not None:
            _action_remove(db_root, remove_id, yes=yes)
            return

        # Default action: install one or more databases.
        _action_install(
            db_root,
            effective_registry_url,
            database_ids,
            refresh=refresh_registry,
            force=force,
            verify_checksum=not no_checksum,
            show_progress=not quiet,
            check_space=not skip_checks,
        )
    except KaryoscopeError as e:
        # Convert known KaryoScope errors to clean user-facing messages.
        raise click.ClickException(str(e)) from e
