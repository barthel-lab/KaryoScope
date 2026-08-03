"""``karyoscope info`` — inspect databases and installation state.

Three modes:

* ``karyoscope info`` (no arg): list installed databases with one-line
  summaries. Equivalent to ``karyoscope download --status`` but with
  more detail (size on disk, manifest summary).
* ``karyoscope info <DATABASE_ID>``: detailed report on a single
  installed database — the parsed manifest, feature-set counts from
  ``hierarchy.tsv``, size on disk, and registry provenance.
* ``karyoscope info <path>``: probe a filesystem path. If it's a
  directory containing a ``manifest.yaml``, treat it like a database
  (useful for inspecting databases not installed via ``download``).

The output is plain text. Machine-readable output (``--json``) is
deferred to a later stage when there's a concrete consumer.
"""

from __future__ import annotations

import contextlib
import logging
from pathlib import Path

import click

from karyoscope import installed as _installed
from karyoscope import paths
from karyoscope.commands._options import resolve_db_root_flag
from karyoscope.core.io.hierarchy import (
    HierarchyError,
    parse_hierarchy,
    validate_hierarchy,
)
from karyoscope.exceptions import (
    DatabaseLayoutError,
    DatabaseNotFoundError,
    KaryoscopeError,
    ManifestError,
)
from karyoscope.manifest import Manifest, validate_database_layout

logger = logging.getLogger(__name__)


# --- formatting helpers ----------------------------------------------------


def _format_size(num_bytes: float) -> str:
    """Format a byte count using SI-style units (KB, MB, GB)."""
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if num_bytes < 1024.0 or unit == "TB":
            return f"{num_bytes:.1f} {unit}"
        num_bytes /= 1024.0
    return f"{num_bytes:.1f} TB"  # pragma: no cover (unreachable)


def _dir_size(path: Path) -> int:
    """Recursively compute the total size of files under ``path``."""
    total = 0
    for f in path.rglob("*"):
        if f.is_file():
            with contextlib.suppress(OSError):
                total += f.stat().st_size
    return total


# --- listing all installed databases ---------------------------------------


def _list_installed(db_root: Path) -> None:
    """Default mode: list all installed databases with brief details."""
    state = _installed.load(db_root)
    click.echo(f"Database root: {db_root}")

    if not db_root.exists():
        click.echo("(root does not exist yet; nothing installed)")
        return

    if not state.databases:
        # Fall back to filesystem listing in case installed.json is missing
        # — a user might have unpacked a tarball manually.
        loose = paths.installed_databases(db_root)
        if not loose:
            click.echo("No installed databases found.")
            click.echo("Run `karyoscope download --list` to see what's available.")
            return
        click.echo(
            f"\nFound {len(loose)} database(s) on disk but none recorded in "
            "installed.json (perhaps unpacked manually):"
        )
        for db in loose:
            click.echo(f"  {db.name}  ({_format_size(_dir_size(db))})")
        return

    click.echo(f"\nInstalled databases ({len(state.databases)}):")
    for db_id, rec in sorted(state.databases.items()):
        db_dir = db_root / rec.directory
        size = _format_size(_dir_size(db_dir)) if db_dir.is_dir() else "missing"
        click.echo(f"  {db_id}")
        click.echo(f"    Version:   {rec.version}")
        click.echo(f"    Installed: {rec.installed_at}")
        click.echo(f"    Size:      {size}")
        click.echo(f"    Source:    {rec.source_url}")


# --- detailed view of one database -----------------------------------------


def _print_manifest_summary(manifest: Manifest, db_dir: Path) -> None:
    """Print the manifest's contents in a tidy plain-text format."""
    click.echo(f"  Version:                {manifest.version}")
    click.echo(f"  KaryoScope min version: {manifest.karyoscope_min_version}")
    click.echo(
        f"  k-mer:                  size={manifest.kmer.size}, "
        f"type={manifest.kmer.type}, max={manifest.kmer.max_size}"
    )
    click.echo(f"  Index type:             {manifest.index.type}")
    if manifest.index.basename:
        click.echo(f"  Index basename:         {manifest.index.basename}")
    if manifest.roles:
        roles = ", ".join(f"{k}={v}" for k, v in sorted(manifest.roles.items()))
        click.echo(f"  Roles:                  {roles}")
    else:
        # Roles are optional: only scaffold/centromeres/karyotype consult them.
        # A database without roles (e.g. a cytoband database) is still valid for
        # annotate/bin, so we note the absence rather than treating it as an error.
        click.echo(
            "  Roles:                  (none declared — fine for annotate/bin; "
            "scaffold/centromeres/karyotype fall back to default feature-set names)"
        )
    if manifest.smoothing:
        for k, v in sorted(manifest.smoothing.items()):
            click.echo(f"  Smoothing ({k}): {v}")
    click.echo(f"  Size on disk:           {_format_size(_dir_size(db_dir))}")


def _print_hierarchy_summary(manifest: Manifest, db_dir: Path) -> None:
    """Print feature-set counts derived from hierarchy.tsv.

    Also runs :func:`validate_hierarchy` and prints any issues as
    warnings. ``info`` is read-only inspection — we don't refuse to
    show what we can just because the hierarchy is malformed. The
    ``annotate`` command treats the same issues as hard errors when
    smoothing is requested.
    """
    hierarchy_path = db_dir / manifest.hierarchy
    try:
        hierarchy = parse_hierarchy(hierarchy_path)
    except HierarchyError as e:
        click.echo(f"  Hierarchy:              <failed to parse: {e}>")
        return

    # Cross-validate against features.tsv when available — surfaces the
    # most useful "this database is internally inconsistent" cases.
    features_columns: dict[str, set[str]] | None = None
    try:
        from karyoscope.core.io.features import parse_features

        features = parse_features(db_dir / manifest.features)
        features_columns = {
            fs: {row[fs] for row in features.table.values()} for fs in features.feature_sets
        }
    except Exception:
        # If features.tsv is missing or malformed the manifest-validation
        # step already complained; just skip the cross-check here.
        features_columns = None

    issues = validate_hierarchy(hierarchy, feature_columns=features_columns)

    # Cross-validate colors.tsv against the hierarchy. info is
    # read-only inspection; treat both checks as warnings rather
    # than erroring. ``download`` and ``karyotype`` enforce the
    # same checks as hard errors.
    color_issues: list[str] = []
    try:
        from karyoscope.core.io.colors import (
            parse_colors_and_groups,
            validate_colors,
            validate_legend_groups,
        )

        colors, legend_groups = parse_colors_and_groups(db_dir / manifest.colors)
        color_issues = validate_colors(hierarchy, colors)
        color_issues += validate_legend_groups(colors, legend_groups)
    except Exception as e:
        color_issues = [f"colors.tsv could not be parsed: {e}"]

    counts = hierarchy.count_by_feature_set()
    click.echo("  Feature sets:")
    for fs in manifest.feature_sets:
        n = counts.get(fs, 0)
        click.echo(f"    {fs}: {n} edges")

    if issues:
        click.echo("  Hierarchy warnings:")
        for issue in issues:
            click.echo(f"    ! {issue}")
    if color_issues:
        click.echo("  Colors warnings:")
        for issue in color_issues:
            click.echo(f"    ! {issue}")


def _show_database(db_root: Path, db_id: str) -> None:
    """Detailed report on a single installed database."""
    state = _installed.load(db_root)
    record = state.databases.get(db_id)
    if record is None:
        raise DatabaseNotFoundError(
            f"database {db_id!r} is not installed. "
            "Run `karyoscope info` (no argument) to list installed databases."
        )

    db_dir = db_root / record.directory
    if not db_dir.is_dir():
        raise DatabaseNotFoundError(
            f"database {db_id!r} is recorded in installed.json but its directory "
            f"{db_dir} is missing."
        )

    click.echo(f"{db_id}")
    click.echo(f"  Path:                   {db_dir}")
    click.echo(f"  Installed at:           {record.installed_at}")
    click.echo(f"  Source URL:             {record.source_url}")
    if record.registry_doi:
        click.echo(f"  Registry DOI:           {record.registry_doi}")

    # Validate the layout and print manifest details. Any layout problems
    # become user-visible errors rather than crashes.
    manifest = validate_database_layout(db_dir)
    _print_manifest_summary(manifest, db_dir)
    _print_hierarchy_summary(manifest, db_dir)


# --- inspecting an ad-hoc path --------------------------------------------


def _show_path(path: Path) -> None:
    """Probe an arbitrary filesystem path."""
    click.echo(f"Path: {path}")
    if not path.exists():
        click.echo("  Does not exist.")
        return

    if path.is_dir():
        manifest_path = path / "manifest.yaml"
        if manifest_path.is_file():
            click.echo("  Type: KaryoScope database directory")
            try:
                manifest = validate_database_layout(path)
            except (ManifestError, DatabaseLayoutError) as e:
                click.echo(f"  Layout valid: NO ({e})")
                click.echo(f"  Size on disk: {_format_size(_dir_size(path))}")
                return
            click.echo(f"  Database id:  {manifest.id}")
            _print_manifest_summary(manifest, path)
            _print_hierarchy_summary(manifest, path)
        else:
            click.echo("  Type: directory")
            click.echo(f"  Size: {_format_size(_dir_size(path))}")
        return

    if path.is_file():
        size = _format_size(path.stat().st_size)
        click.echo(f"  Type: file ({size})")
        return

    click.echo("  Type: special (not a regular file or directory)")


# --- click command ---------------------------------------------------------


@click.command(
    help="Inspect installed databases or a given path.",
    no_args_is_help=False,
)
@click.argument("target", required=False)
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
def cmd(target: str | None, db_root_arg: Path | None, db_alias: Path | None) -> None:
    """Show information about installed databases or a given path.

    \b
    Examples:
        karyoscope info                     # list installed databases
        karyoscope info KS_human_CHM13_v2   # details on one database
        karyoscope info ./some/path/        # probe a filesystem path
    """
    db_root = paths.ensure_db_root(resolve_db_root_flag(db_root_arg, db_alias, command="info"))

    try:
        if target is None:
            _list_installed(db_root)
            return

        # If the argument looks like an existing path, treat it as one.
        # Otherwise assume it's a database id.
        as_path = Path(target).expanduser()
        if as_path.exists():
            _show_path(as_path)
            return
        if "/" in target or target.startswith("."):
            # Path-shaped but doesn't exist — let the user know.
            _show_path(as_path)
            return

        _show_database(db_root, target)
    except (
        DatabaseNotFoundError,
        DatabaseLayoutError,
        ManifestError,
        KaryoscopeError,
    ) as e:
        raise click.ClickException(str(e)) from e
