"""``karyoscope info`` — inspect databases, files, and installation state.

In v0.1.0.dev this command supports a limited subset of the planned behavior:

* ``karyoscope info`` (no argument) — show installed databases and the default
  database root.
* ``karyoscope info databases`` — same as above (explicit form).
* ``karyoscope info <path>`` — basic info about a path on disk (existence,
  type, and whether it looks like an installed KaryoScope database).

Richer introspection (parsing BED headers, inspecting KMC indexes, etc.) will
land alongside the corresponding subcommands.
"""

from __future__ import annotations

import contextlib
from pathlib import Path

import click

from karyoscope.paths import default_db_root, installed_databases


def _format_size(num_bytes: int) -> str:
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


def _show_installed_databases() -> None:
    root = default_db_root()
    dbs = installed_databases(root)

    click.echo(f"Default database root: {root}")
    if not root.exists():
        click.echo("(root does not exist yet; nothing installed)")
        return

    if not dbs:
        click.echo("No installed databases found.")
        click.echo("Run `karyoscope download --list` to see what's available.")
        return

    click.echo(f"\nInstalled databases ({len(dbs)}):")
    for db in dbs:
        size = _format_size(_dir_size(db))
        click.echo(f"  {db.name}  ({size})")


def _show_path(path: Path) -> None:
    click.echo(f"Path: {path}")
    if not path.exists():
        click.echo("  Does not exist.")
        return

    if path.is_dir():
        manifest = path / "manifest.yaml"
        if manifest.is_file():
            click.echo("  Type: KaryoScope database directory")
            click.echo(f"  Manifest: {manifest}")
            click.echo(f"  Total size: {_format_size(_dir_size(path))}")
        else:
            click.echo("  Type: directory")
            click.echo(f"  Total size: {_format_size(_dir_size(path))}")
        return

    if path.is_file():
        size = _format_size(path.stat().st_size)
        click.echo(f"  Type: file ({size})")
        return

    click.echo("  Type: special (not a regular file or directory)")


@click.command(
    help="Inspect installed databases, file paths, or KaryoScope installation state.",
    no_args_is_help=False,
)
@click.argument("target", required=False)
def cmd(target: str | None) -> None:
    """Show information about installed databases or a given path.

    Examples:

    \b
        karyoscope info
        karyoscope info databases
        karyoscope info ~/.karyoscope/db/KS_human_CHM13_v2/
        karyoscope info my_assembly.chromosome.smoothed.bed.gz
    """
    if target is None or target == "databases":
        _show_installed_databases()
        return

    _show_path(Path(target).expanduser())
