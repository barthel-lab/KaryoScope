"""Console-script entry point, deliberately free of heavy imports.

``karyoscope.cli`` eagerly imports every command module so it can register
the subcommands, and those pull in the package's hard dependencies
(``drawsvg``, ``cairosvg``, ``requests``, ...). When one of them is
missing — a partial install, or a shell using a different interpreter than
the one KaryoScope was installed into — the failure happens *during*
``from karyoscope.cli import main``. That import lives in the console
script pip generates:

    from karyoscope.cli import main
    sys.exit(main())

so it is outside any code KaryoScope controls. The result is a raw
``ModuleNotFoundError`` traceback from **every** command, including
``karyoscope --help`` and ``karyoscope version`` — precisely the two a
user would reach for to diagnose it. The dependency preflight in
:mod:`karyoscope.preflight` cannot help either: it runs inside a command,
long after import has already failed.

Routing the entry point through this module puts a ``try`` around that
import, so a broken install reports itself in a few readable lines. The
message names the interpreter in use, because "wrong Python" is the most
common cause and the least obvious from the traceback alone.
"""

from __future__ import annotations

import sys


def _report_broken_install(exc: ImportError) -> None:
    """Explain an import failure that happened before any command could run."""
    missing = getattr(exc, "name", None) or "a required package"
    print(
        f"Error: KaryoScope's installation is incomplete — could not import "
        f"{missing!r}.\n"
        f"\n"
        f"  KaryoScope itself is importable, so the package is present but at "
        f"least one\n"
        f"  of its dependencies is not — or this shell is using a different "
        f"Python than\n"
        f"  the one KaryoScope was installed into.\n"
        f"\n"
        f"  Interpreter in use: {sys.executable}\n"
        f"\n"
        f"  Fix with one of:\n"
        f"    pip install -e .          # from a KaryoScope checkout\n"
        f"    pip install karyoscope\n"
        f"    conda activate karyoscope # if you meant to use the env\n"
        f"\n"
        f"  Original error: {exc}",
        file=sys.stderr,
    )


def main() -> int:
    """Import and run the CLI, reporting a broken install readably."""
    try:
        from karyoscope.cli import main as run
    except ImportError as exc:
        _report_broken_install(exc)
        return 1
    return run()
