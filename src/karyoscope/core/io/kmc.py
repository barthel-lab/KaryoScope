"""Locate and invoke the ``get_featureIDs`` C++ binary.

The C++ helper is what actually queries KMC indexes; this module is the
thin Python adapter that calls it via :mod:`karyoscope.core.external` so
the rest of the codebase doesn't have to know anything about subprocess
plumbing.

Binary lookup order
===================

1. ``$KARYOSCOPE_GET_FEATUREIDS`` — explicit override. Useful in test
   environments, custom installs, or when developing two builds in
   parallel.
2. ``shutil.which("get_featureIDs")`` — anything on ``$PATH``. This is
   the path Bioconda users will hit (the package places the binary in
   the env's ``bin/``).
3. ``<repo>/native/get_featureIDs/build/get_featureIDs`` — the editable
   dev install path. The package source lives at
   ``<repo>/src/karyoscope/...``; we walk up from this file's location
   to find the repo root.

If none of those resolve, :func:`get_featureids_binary` raises
:class:`ToolNotFoundError` with a message telling the user exactly what
to do.
"""

from __future__ import annotations

import logging
import os
import shutil
from pathlib import Path

from karyoscope.core.external import ToolNotFoundError, run_tool

logger = logging.getLogger(__name__)

#: The bare command name we search for on ``$PATH``.
BINARY_NAME = "get_featureIDs"

#: Environment variable users can set to override the binary location.
ENV_OVERRIDE = "KARYOSCOPE_GET_FEATUREIDS"


def _editable_install_candidate() -> Path | None:
    """Try to locate a developer-built binary in the source tree.

    Walks up from this file's location looking for
    ``native/get_featureIDs/build/get_featureIDs``. Returns ``None`` if
    no such path exists — e.g., when running from a wheel install.
    """
    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / "native" / "get_featureIDs" / "build" / "get_featureIDs"
        if candidate.is_file():
            return candidate
        # Stop walking once we reach a filesystem root; otherwise we could
        # walk forever on weird paths.
        if parent == parent.parent:
            return None
    return None


def get_featureids_binary() -> str:
    """Return the path to the ``get_featureIDs`` binary.

    Raises
    ------
    ToolNotFoundError
        If the binary cannot be found via any of the documented lookup
        mechanisms. The error message includes instructions for both
        the source-build and Bioconda paths.
    """
    # 1. Explicit env override.
    env_value = os.environ.get(ENV_OVERRIDE)
    if env_value:
        env_path = Path(env_value)
        if env_path.is_file():
            logger.debug("using %s from %s=%s", BINARY_NAME, ENV_OVERRIDE, env_path)
            return str(env_path)
        # An explicitly-set-but-bad path is a user error worth surfacing
        # rather than silently falling through.
        raise ToolNotFoundError(
            f"{ENV_OVERRIDE} is set to {env_value!r}, but no file exists there. "
            f"Either unset {ENV_OVERRIDE} or point it at the {BINARY_NAME} binary."
        )

    # 2. On $PATH.
    path_binary = shutil.which(BINARY_NAME)
    if path_binary:
        logger.debug("using %s from $PATH at %s", BINARY_NAME, path_binary)
        return path_binary

    # 3. Built in the source tree (editable install).
    src_binary = _editable_install_candidate()
    if src_binary is not None:
        logger.debug("using %s from source tree at %s", BINARY_NAME, src_binary)
        return str(src_binary)

    # Out of options.
    raise ToolNotFoundError(
        f"{BINARY_NAME!r} was not found. To fix this, choose one of:\n"
        "  - Build it from source:\n"
        "      make -C native/get_featureIDs\n"
        "  - Install KaryoScope via Bioconda (planned for v1.0):\n"
        "      conda install -c bioconda karyoscope\n"
        f"  - Set the {ENV_OVERRIDE} environment variable to an explicit "
        "binary path."
    )


def run_get_featureids(
    *,
    db_path: Path,
    input_path: Path,
    output_dir: Path,
    threads: int = 0,
    prefix: str | None = None,
    input_format: str | None = None,
    capture: bool = False,
) -> Path:
    """Invoke ``get_featureIDs`` with the given options.

    Parameters
    ----------
    db_path
        Path to the KMC database. May point at the basename (without
        ``.kmc_pre``/``.kmc_suf`` suffix), at one of the index files
        directly, or at the directory containing the pair.
    input_path
        FASTA/FASTQ input file. Pass ``Path("-")`` for stdin (and set
        ``input_format``).
    output_dir
        Directory where the combined BED will be written. Created if
        absent.
    threads
        Number of worker threads. ``0`` (the default) means "let
        get_featureIDs decide based on hardware concurrency."
    prefix
        Optional override for the output filename prefix. If ``None``,
        ``get_featureIDs`` derives the prefix from the input FASTA
        basename plus the KMC database basename.
    input_format
        ``"fasta"`` or ``"fastq"``. Required only when reading from
        stdin; otherwise the format is inferred from the input filename.
    capture
        If ``True``, capture the subprocess stdout/stderr instead of
        passing through. Useful for tests.

    Returns
    -------
    Path
        Path to the produced ``*.combined.presmoothed.featureIDs.bed``
        file.

    Raises
    ------
    ToolNotFoundError
        If the binary cannot be located.
    ExternalToolError
        If the subprocess exits with a non-zero status.
    """
    binary = get_featureids_binary()

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd: list[str] = [
        binary,
        "--db",
        str(db_path),
        "--input",
        str(input_path),
        "--output",
        str(output_dir),
        "--threads",
        str(threads),
    ]
    if prefix is not None:
        cmd += ["--prefix", prefix]
    if input_format is not None:
        cmd += ["--input-format", input_format]

    run_tool(cmd, capture=capture)

    # Reconstruct the expected output filename so callers don't need to
    # know the internal naming convention.
    if prefix is not None:
        base = prefix
    else:
        # Mirrors get_featureIDs' get_fasta_prefix + kmc basename logic.
        # Strip known FASTA/FASTQ extensions in the same order the C++ does.
        input_basename = input_path.name
        for ext in (
            ".fasta.gz",
            ".fa.gz",
            ".fna.gz",
            ".fastq.gz",
            ".fq.gz",
            ".fasta",
            ".fa",
            ".fna",
            ".fastq",
            ".fq",
        ):
            if input_basename.lower().endswith(ext):
                input_basename = input_basename[: -len(ext)]
                break
        kmc_basename = Path(str(db_path)).name
        # If db_path pointed at a file like "features.kmc_pre", strip the suffix.
        for suffix in (".kmc_pre", ".kmc_suf"):
            if kmc_basename.endswith(suffix):
                kmc_basename = kmc_basename[: -len(suffix)]
                break
        base = f"{input_basename}.{kmc_basename}"

    return output_dir / f"{base}.combined.presmoothed.featureIDs.bed"
