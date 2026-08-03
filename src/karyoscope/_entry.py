"""Console-script entry point, deliberately free of heavy imports.

When one of KaryoScope's dependencies is missing — a partial install, or a
shell using a different interpreter than the one KaryoScope was installed
into — the failure would otherwise surface as a raw ``ModuleNotFoundError``
traceback from the console script pip generates:

    from karyoscope.cli import main
    sys.exit(main())

That import is outside any code KaryoScope controls, so the traceback comes
back from **every** command, including ``karyoscope --help`` and
``karyoscope version`` — precisely the two a user would reach for to diagnose
it. The dependency preflight in :mod:`karyoscope.preflight` cannot help
either: it runs inside a command, long after import has already failed.

Routing the entry point through this module makes a broken install report
itself in a few readable lines. The message names the interpreter in use,
because "wrong Python" is the most common cause and the least obvious from
the traceback alone.

**Why the check is explicit.** It used to be inferred: ``karyoscope.cli``
imported every command module eagerly, those pulled in the hard dependencies,
and wrapping that one import therefore caught a missing one. Subcommands are
now imported on demand — registering them cost ~190 ms of startup for an HTTP
library most commands never touch — so that side-effect is gone, and a missing
dependency would otherwise not surface until mid-dispatch, as exactly the bare
traceback this module exists to prevent. The dependencies are now checked
directly, which is also stricter: the old form caught only what the eager
imports happened to touch, while the list below is exact.
"""

from __future__ import annotations

import importlib.util
import sys

#: Dependencies to verify when the installed metadata cannot be read (e.g. a
#: source checkout that was never pip-installed). Keep in step with
#: ``pyproject.toml``; the metadata path below is preferred precisely because
#: it cannot drift.
_FALLBACK_REQUIREMENTS: tuple[str, ...] = (
    "click",
    "drawsvg",
    "cairosvg",
    "requests",
    "yaml",
    "tqdm",
    "jsonschema",
)

#: Distribution name -> module name, for the cases where they differ.
_MODULE_NAMES: dict[str, str] = {"pyyaml": "yaml"}


def _required_modules() -> list[str]:
    """Module names for our declared runtime dependencies.

    Read from installed metadata so this cannot drift from ``pyproject.toml``.
    Requirements guarded by an ``extra`` marker are skipped: a missing optional
    extra is not a broken install, and reporting it as one would send the user
    chasing a package they never asked for.
    """
    try:
        from importlib.metadata import requires

        declared = requires("karyoscope")
    except Exception:
        declared = None
    if not declared:
        return list(_FALLBACK_REQUIREMENTS)

    modules: list[str] = []
    for raw in declared:
        if "extra ==" in raw:
            continue
        # "requests>=2.31" / "pyyaml (>=6.0)" / "tqdm" -> the bare name.
        name = raw.split(";")[0].strip()
        for sep in ("[", "(", ">", "<", "=", "!", "~", " "):
            name = name.split(sep)[0]
        name = name.strip().lower()
        if name:
            modules.append(_MODULE_NAMES.get(name, name.replace("-", "_")))
    return modules or list(_FALLBACK_REQUIREMENTS)


def _missing_dependencies() -> list[str]:
    """Declared dependencies that are not installed.

    Uses :func:`importlib.util.find_spec`, which locates a module without
    executing it — the point of the exercise, since importing them is the cost
    being avoided. Measured at ~8 ms for the full set, against ~190 ms to
    import ``requests`` alone.
    """
    missing: list[str] = []
    for module in _required_modules():
        try:
            if importlib.util.find_spec(module) is None:
                missing.append(module)
        except (ImportError, ValueError):
            # A parent package that is itself broken -- treat as missing.
            missing.append(module)
    return missing


def _report_broken_install(missing: str, detail: str) -> None:
    """Explain a dependency failure that happened before any command could run."""
    print(
        f"Error: KaryoScope's installation is incomplete — could not import "
        f"{missing}.\n"
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
        f"  Original error: {detail}",
        file=sys.stderr,
    )


def main() -> int:
    """Run the CLI, reporting a broken install readably."""
    missing = _missing_dependencies()
    if missing:
        names = ", ".join(repr(m) for m in missing)
        _report_broken_install(names, f"missing module(s): {names}")
        return 1

    # find_spec proves a module is present, not that it imports cleanly: a
    # truncated install, or a compiled dependency built against the wrong ABI
    # (cairosvg's C bindings being the realistic case), is found and then
    # fails. Both the CLI import and the run are wrapped, because subcommands
    # are now imported during dispatch rather than up front -- so an import
    # failure can surface after `run()` has already been entered.
    try:
        from karyoscope.cli import main as run

        return run()
    except ImportError as exc:
        name = getattr(exc, "name", None)
        _report_broken_install(repr(name) if name else "a required package", str(exc))
        return 1
