"""Up-front checks that a command's external dependencies are available.

KaryoScope shells out to a handful of binaries, and imports a couple of
Python packages that wrap native libraries. Historically each was resolved
at its point of use — which means a `karyotype` run could spend twenty
minutes annotating and then fail because ``seqtk`` was never installed, or
an `annotate` could finish the k-mer query and die at the bgzip step.

This module resolves everything a command needs *before* it starts, and
reports all missing dependencies in one message rather than one per rerun.

The catalogue below is the single place that knows how to obtain each
dependency, so the install hints stay consistent with ``environment.yml``
and the README. :func:`require` is what commands call; :func:`check`
returns the missing set without raising, for reporting contexts like
``karyoscope version``.
"""

from __future__ import annotations

import importlib.util
import logging
import os
import re
import shutil
import subprocess
from collections.abc import Callable
from dataclasses import dataclass

from karyoscope.exceptions import IncompatibleVersionError, MissingDependencyError
from karyoscope.versions import at_least

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Dependency:
    """One external dependency KaryoScope may need.

    ``kind`` is ``"binary"`` (looked up on ``$PATH``, unless the name has
    an entry in :data:`_RESOLVERS`) or ``"python"`` (looked up as an
    importable module).

    ``min_version`` is the oldest release KaryoScope can drive, for the
    binaries whose interface we depend on. It is only enforced when the
    installed version can actually be read — see :func:`outdated`.
    """

    name: str
    kind: str
    purpose: str
    install_hint: str
    min_version: str | None = None


#: Every dependency KaryoScope can require, keyed by the name commands use.
DEPENDENCIES: dict[str, Dependency] = {
    "hks": Dependency(
        name="hks",
        kind="binary",
        purpose="querying and building HKS-backend databases",
        install_hint=(
            "Build it from source and put it on PATH:\n"
            "    git clone --recurse-submodules https://github.com/jnalanko/HKS.git\n"
            '    cargo install --path HKS --root "$CONDA_PREFIX"\n'
            "  or set $KARYOSCOPE_HKS to an existing hks binary."
        ),
        # 0.3.0 added --miss-label and smooth's --no-header, which let hks
        # write KaryoScope's BED format directly. Against 0.2.0 every
        # annotate would die at the first lookup with clap's "unexpected
        # argument '--miss-label'" — true, but it does not say "upgrade hks".
        min_version="0.3.0",
    ),
    "get_featureIDs": Dependency(
        name="get_featureIDs",
        kind="binary",
        purpose="querying KMC-backend databases",
        install_hint=(
            "Compile the bundled helper from the KaryoScope source tree:\n"
            "    make -C native/get_featureIDs\n"
            "  or set $KARYOSCOPE_GET_FEATUREIDS to an existing binary."
        ),
    ),
    "kmc": Dependency(
        name="kmc",
        kind="binary",
        purpose="building KMC-backend databases",
        install_hint="conda install -c bioconda kmc",
    ),
    "bgzip": Dependency(
        name="bgzip",
        kind="binary",
        purpose="compressing output BEDs",
        install_hint=(
            "conda install -c bioconda htslib\n"
            "  or re-run with --no-bgzip to write plain .bed files."
        ),
    ),
    "tabix": Dependency(
        name="tabix",
        kind="binary",
        purpose="indexing compressed BEDs",
        install_hint="conda install -c bioconda htslib",
    ),
    "seqtk": Dependency(
        name="seqtk",
        kind="binary",
        purpose="telomere detection",
        install_hint="conda install -c bioconda seqtk",
    ),
    "samtools": Dependency(
        name="samtools",
        kind="binary",
        purpose="reading BAM input",
        install_hint="conda install -c bioconda samtools",
    ),
    "cairosvg": Dependency(
        name="cairosvg",
        kind="python",
        purpose="rendering PDF/PNG karyotypes",
        install_hint=(
            "pip install cairosvg, and install the native library it wraps:\n"
            "    conda install -c conda-forge cairo\n"
            "  or re-run with --format svg."
        ),
    ),
}


def _hks_binary() -> str:
    from karyoscope.core.io.hks import get_hks_binary

    return get_hks_binary()


def _get_featureids_binary() -> str:
    from karyoscope.core.io.kmc import get_featureids_binary

    return get_featureids_binary()


#: Dependencies whose lookup is more than "is it on $PATH". Both k-mer
#: backends accept an environment override and ``get_featureIDs`` also
#: resolves out of the source tree for editable installs, so we must ask
#: the backend rather than reimplement its search order — a second
#: implementation here would reject working installs.
#:
#: Imported lazily inside each function: :mod:`karyoscope.core.annotate`
#: imports this module, so a top-level import of its siblings would be a
#: cycle.
_RESOLVERS: dict[str, Callable[[], str]] = {
    "hks": _hks_binary,
    "get_featureIDs": _get_featureids_binary,
}


def resolve_binary(name: str) -> str | None:
    """Return the path KaryoScope would use for ``name``, or None if absent.

    The single source of truth for "can KaryoScope run this tool", used by
    the checks below and by ``karyoscope version``. Honours each backend's
    own lookup order rather than assuming ``$PATH``.
    """
    resolver = _RESOLVERS.get(name)
    if resolver is None:
        return shutil.which(name)
    try:
        return resolver()
    except Exception as exc:  # ToolNotFoundError, or an override pointing nowhere
        logger.debug("preflight: %s did not resolve: %s", name, exc)
        return None


#: ``hks`` has no ``--version`` flag; every subcommand logs a banner instead.
_HKS_VERSION_BANNER = re.compile(r"Running hks version\s+(\d[\w.]*)")


def _hks_version(path: str) -> str | None:
    """Read the version out of the banner ``hks`` logs on startup.

    Run with no arguments it prints usage and exits non-zero, but the banner
    is emitted first, which is all we need — so the exit status is ignored.

    The banner goes through ``log::info!``, so ``RUST_LOG`` in the caller's
    environment could suppress it and make a perfectly good hks look
    unidentifiable. We force it back on for this one probe rather than let a
    user's debugging setting change whether KaryoScope runs.
    """
    try:
        proc = subprocess.run(
            [path],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
            env={**os.environ, "RUST_LOG": "info"},
        )
    except (subprocess.TimeoutExpired, OSError) as exc:
        logger.debug("preflight: could not run %s to read its version: %s", path, exc)
        return None

    match = _HKS_VERSION_BANNER.search(proc.stderr) or _HKS_VERSION_BANNER.search(proc.stdout)
    return match.group(1) if match else None


#: How to read the installed version of a binary, keyed as :data:`_RESOLVERS`.
#: A dependency with a ``min_version`` but no probe here is never checked.
_VERSION_PROBES: dict[str, Callable[[str], str | None]] = {
    "hks": _hks_version,
}


def installed_version(name: str) -> str | None:
    """Return the version string of an installed binary, or None if unknown.

    None means "could not determine" — the tool is absent, has no probe, or
    did not identify itself. It never means "too old".
    """
    dep = DEPENDENCIES.get(name)
    probe = _VERSION_PROBES.get(name)
    if dep is None or probe is None:
        return None
    path = resolve_binary(name)
    if path is None:
        return None
    return probe(path)


def _binary_available(dep: Dependency) -> bool:
    return resolve_binary(dep.name) is not None


def _python_available(dep: Dependency) -> bool:
    try:
        return importlib.util.find_spec(dep.name) is not None
    except (ImportError, ValueError):  # pragma: no cover — malformed install
        return False


def _available(dep: Dependency) -> bool:
    return _binary_available(dep) if dep.kind == "binary" else _python_available(dep)


def check(names: list[str]) -> list[Dependency]:
    """Return the dependencies in ``names`` that are not available.

    Unknown names are ignored rather than raising: the caller is asking
    "what's missing", and an unrecognised name is a KaryoScope bug, not a
    user-facing dependency problem.
    """
    missing: list[Dependency] = []
    for name in names:
        dep = DEPENDENCIES.get(name)
        if dep is None:
            logger.debug("preflight: unknown dependency %r, skipping", name)
            continue
        if _available(dep):
            logger.debug("preflight: %s found", name)
        else:
            missing.append(dep)
    return missing


def outdated(names: list[str]) -> list[tuple[Dependency, str]]:
    """Return ``(dependency, installed_version)`` for each one that is too old.

    Only reports a dependency when we have positive evidence of a problem:
    it declares a ``min_version``, we could read what is installed, and that
    is older. An unreadable version is not reported — refusing to run over a
    tool that merely declined to identify itself would be worse than letting
    it try and fail on its own terms.
    """
    stale: list[tuple[Dependency, str]] = []
    for name in names:
        dep = DEPENDENCIES.get(name)
        if dep is None or dep.min_version is None:
            continue
        have = installed_version(name)
        if have is None:
            logger.debug("preflight: could not read %s version; not enforcing minimum", name)
            continue
        if at_least(have, dep.min_version):
            logger.debug("preflight: %s %s satisfies >= %s", name, have, dep.min_version)
        else:
            stale.append((dep, have))
    return stale


def require(names: list[str], *, context: str) -> None:
    """Raise unless every dependency in ``names`` is available and new enough.

    Absence is reported before staleness: a tool that is not installed has no
    version to complain about, and one message about one problem is easier to
    act on than two about the same tool.

    Parameters
    ----------
    names
        Keys into :data:`DEPENDENCIES`.
    context
        What the dependencies are needed for, used in the message — e.g.
        ``"annotate against KS_human_CHM13_v2"``.

    Raises
    ------
    MissingDependencyError
        Listing every missing dependency, with an install hint for each.
    IncompatibleVersionError
        If an installed dependency is older than its ``min_version``.
    """
    missing = check(names)
    if not missing:
        _require_versions(names, context=context)
        return

    plural = "dependencies" if len(missing) > 1 else "dependency"
    lines = [f"missing required {plural} for {context}:"]
    for dep in missing:
        where = "on $PATH" if dep.kind == "binary" else "as a Python module"
        lines.append(f"\n  {dep.name} — not found {where} (needed for {dep.purpose})")
        for hint_line in dep.install_hint.splitlines():
            lines.append(f"    {hint_line}")
    lines.append(
        "\nThe conda environment shipped with KaryoScope (environment.yml) "
        "provides all of these except hks and get_featureIDs, which are "
        "built separately — see the Installation section of the README.\n"
        "Run `karyoscope version` to see what KaryoScope can currently find."
    )
    raise MissingDependencyError("\n".join(lines))


def _require_versions(names: list[str], *, context: str) -> None:
    """Raise :class:`IncompatibleVersionError` if any dependency is too old."""
    stale = outdated(names)
    if not stale:
        return

    lines = [f"outdated dependency for {context}:"]
    for dep, have in stale:
        lines.append(
            f"\n  {dep.name} {have} is installed, but {dep.min_version} or newer "
            f"is required (needed for {dep.purpose})"
        )
        for hint_line in dep.install_hint.splitlines():
            lines.append(f"    {hint_line}")
    lines.append("\nRun `karyoscope version` to see the versions KaryoScope can see.")
    raise IncompatibleVersionError("\n".join(lines))
