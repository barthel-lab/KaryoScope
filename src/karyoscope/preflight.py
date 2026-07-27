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
import shutil
from collections.abc import Callable
from dataclasses import dataclass

from karyoscope.exceptions import MissingDependencyError

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Dependency:
    """One external dependency KaryoScope may need.

    ``kind`` is ``"binary"`` (looked up on ``$PATH``, unless the name has
    an entry in :data:`_RESOLVERS`) or ``"python"`` (looked up as an
    importable module).
    """

    name: str
    kind: str
    purpose: str
    install_hint: str


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


def require(names: list[str], *, context: str) -> None:
    """Raise unless every dependency in ``names`` is available.

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
    """
    missing = check(names)
    if not missing:
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
