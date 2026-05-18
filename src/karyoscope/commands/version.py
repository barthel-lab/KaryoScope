"""``karyoscope version`` — print version and environment info.

This is the detailed form of ``karyoscope --version``. It reports the
KaryoScope version, the Python interpreter, the install location, the
default database root, and the presence and versions of the external
command-line tools KaryoScope depends on. The output is suitable for
pasting into bug reports.
"""

from __future__ import annotations

import platform
import shutil
import subprocess
import sys
from importlib.metadata import PackageNotFoundError, version

import click

from karyoscope._version import __version__
from karyoscope.paths import default_db_root, installed_databases

#: External tools KaryoScope shells out to. Each entry is
#: (display_name, executable_name, --version flag).
_EXTERNAL_TOOLS: tuple[tuple[str, str, str], ...] = (
    ("KMC", "kmc", "--version"),
    ("bgzip", "bgzip", "--version"),
    ("tabix", "tabix", "--version"),
    ("seqtk", "seqtk", ""),
)

#: Python dependencies whose versions we report.
_PYTHON_DEPS: tuple[str, ...] = (
    "click",
    "drawsvg",
    "cairosvg",
    "requests",
    "pyyaml",
    "tqdm",
    "jsonschema",
)


def _python_dep_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not installed"


def _external_tool_version(executable: str, version_flag: str) -> tuple[str | None, str | None]:
    """Return (path, version_string) for an external tool, or (None, None) if absent.

    The version string is the first non-empty line of output from running the
    tool with ``version_flag``. If ``version_flag`` is empty, the tool is run
    with no arguments (covers tools like seqtk that print usage with a version
    in their header).
    """
    path = shutil.which(executable)
    if path is None:
        return None, None

    cmd = [executable, version_flag] if version_flag else [executable]
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (subprocess.TimeoutExpired, OSError):
        return path, "(could not determine version)"

    # Many of these tools print version info to stderr rather than stdout,
    # and some emit a multi-line usage message; take the first non-empty line.
    for stream in (proc.stdout, proc.stderr):
        for line in stream.splitlines():
            line = line.strip()
            if line:
                return path, line
    return path, "(no version output)"


@click.command(
    help="Print KaryoScope version and environment info (useful for bug reports).",
    no_args_is_help=False,
)
def cmd() -> None:
    """Show detailed version and environment information."""
    click.echo(f"KaryoScope {__version__}")
    click.echo(f"  Python: {sys.version.split()[0]} ({sys.executable})")
    click.echo(f"  Platform: {platform.platform()}")

    db_root = default_db_root()
    n_installed = len(installed_databases(db_root))
    click.echo(f"  Default database root: {db_root}")
    click.echo(f"  Installed databases: {n_installed}")

    click.echo("\nPython dependencies:")
    for name in _PYTHON_DEPS:
        click.echo(f"  {name}: {_python_dep_version(name)}")

    click.echo("\nExternal tools:")
    for display_name, executable, version_flag in _EXTERNAL_TOOLS:
        path, ver = _external_tool_version(executable, version_flag)
        if path is None:
            click.echo(f"  {display_name}: not found on PATH")
        else:
            click.echo(f"  {display_name}: {ver} (at {path})")
