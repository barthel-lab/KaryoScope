"""Locate and invoke the ``hks`` binary for k-mer lookup queries.

HKS (Hierarchical K-mer Sets) is an alternative k-mer index backend that
stores named labels directly in the index, eliminating the integer feature-ID
translation layer required by KMC. Where the KMC backend queries one combined
index and translates integer feature ids through ``features.tsv``, the HKS
backend queries one ``.hksf`` file per feature set and reads label names
straight out of the ``hks lookup`` output.

Binary lookup order
===================

1. ``$KARYOSCOPE_HKS`` — explicit override.
2. ``shutil.which("hks")`` — anything on ``$PATH``.

If neither resolves, :func:`get_hks_binary` raises :class:`ToolNotFoundError`.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from karyoscope.core.external import (
    ExternalToolError,
    ToolNotFoundError,
    require_tool,
    run_tool,
)
from karyoscope.core.io.features import NOVEL_NAME

logger = logging.getLogger(__name__)

BINARY_NAME = "hks"
ENV_OVERRIDE = "KARYOSCOPE_HKS"

#: Label HKS emits for k-mers not found in the index.
_HKS_MISS_LABEL = "none"

#: Default max-gap (bases) passed to ``hks smooth``, matching karyoscope's Python default.
_DEFAULT_SMOOTH_MAX_GAP = 1000

_INPUT_EXTENSIONS: tuple[str, ...] = (
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
    ".bam",
)


def get_hks_binary() -> str:
    """Return the path to the ``hks`` binary.

    Raises
    ------
    ToolNotFoundError
        If the binary cannot be found.
    """
    env_value = os.environ.get(ENV_OVERRIDE)
    if env_value:
        env_path = Path(env_value)
        if env_path.is_file():
            logger.debug("using %s from %s=%s", BINARY_NAME, ENV_OVERRIDE, env_path)
            return str(env_path)
        raise ToolNotFoundError(
            f"{ENV_OVERRIDE} is set to {env_value!r}, but no file exists there. "
            f"Either unset {ENV_OVERRIDE} or point it at the {BINARY_NAME} binary."
        )

    path_binary = shutil.which(BINARY_NAME)
    if path_binary:
        logger.debug("using %s from $PATH at %s", BINARY_NAME, path_binary)
        return path_binary

    raise ToolNotFoundError(
        f"{BINARY_NAME!r} was not found on PATH. "
        f"Build HKS from source (cargo build --release in the HKS repo) "
        f"and place the binary on PATH, or set {ENV_OVERRIDE} to its location."
    )


def _infer_prefix(input_path: Path, db_basename: str) -> str:
    """Compute the output filename prefix for HKS lookup output."""
    input_basename = input_path.name
    for ext in _INPUT_EXTENSIONS:
        if input_basename.lower().endswith(ext):
            input_basename = input_basename[: -len(ext)]
            break
    return f"{input_basename}.{db_basename}"


def convert_hks_tsv_to_bed(tsv_path: Path, bed_path: Path) -> None:
    """Convert HKS TSV output (with header) to a headerless BED file.

    Strips the header line and replaces the HKS miss label ``none`` with
    KaryoScope's ``novel`` sentinel.
    """
    lines = tsv_path.read_text().splitlines(keepends=True)
    with bed_path.open("w") as out:
        for line in lines[1:]:  # skip header
            out.write(line.replace(f"\t{_HKS_MISS_LABEL}\n", f"\t{NOVEL_NAME}\n"))


def run_hks_lookup(
    *,
    base_path: Path,
    feature_set_file: Path,
    k: int,
    input_path: Path,
    output_path: Path,
    threads: int = 0,
    capture: bool = False,
) -> Path:
    """Invoke ``hks lookup`` for one feature set.

    Parameters
    ----------
    base_path
        Path to the HKS base index file (``*.hksb``).
    feature_set_file
        Path to the HKS feature-set file (``*.hksf``) for this feature set.
    k
        K-mer length to query.
    input_path
        FASTA or FASTQ input (plain or gzipped). For BAM inputs, pass the
        ``.bam`` path — it will be converted via ``samtools fasta`` internally.
    output_path
        Where to write the raw HKS lookup TSV (header + ``none`` for misses).
    threads
        Number of worker threads (0 = let HKS decide, effectively ``4`` by default).
    capture
        If ``True``, capture subprocess stdout/stderr instead of passing through.

    Returns
    -------
    Path
        ``output_path``, after writing.

    Raises
    ------
    ToolNotFoundError
        If ``hks`` (or ``samtools`` for BAM inputs) is not found.
    ExternalToolError
        If the subprocess exits with a non-zero status.
    """
    binary = get_hks_binary()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if input_path.suffix.lower() == ".bam":
        return _run_hks_lookup_from_bam(
            binary=binary,
            base_path=base_path,
            feature_set_file=feature_set_file,
            k=k,
            bam_path=input_path,
            output_path=output_path,
            threads=threads,
            capture=capture,
        )

    n_threads = threads if threads > 0 else 4
    cmd: list[str] = [
        binary,
        "lookup",
        "-i", str(base_path),
        "--feature-set-file", str(feature_set_file),
        "-k", str(k),
        "-q", str(input_path),
        "--report-query-names",
        "--report-misses",
        "-t", str(n_threads),
        "-o", str(output_path),
    ]
    logger.debug("running: %s", " ".join(cmd))

    run_tool(cmd, capture=capture)
    return output_path


def _run_hks_lookup_from_bam(
    *,
    binary: str,
    base_path: Path,
    feature_set_file: Path,
    k: int,
    bam_path: Path,
    output_path: Path,
    threads: int,
    capture: bool,
) -> Path:
    """Convert a BAM to FASTA in a temp file and run HKS lookup on it.

    HKS requires a seekable file path, so we materialise the samtools
    output to a temporary FASTA rather than streaming via a named pipe.
    """
    samtools = require_tool(
        "samtools",
        install_hint=(
            "Install samtools to use BAM inputs:\n"
            "  conda install -c bioconda samtools\n"
            "Or convert the BAM to FASTA first:\n"
            "  samtools fasta input.bam | gzip > input.fasta.gz"
        ),
    )

    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp:
        tmp_fasta = Path(tmp.name)

    try:
        logger.debug("converting BAM to FASTA: %s -> %s", bam_path, tmp_fasta)
        samtools_result = subprocess.run(
            [samtools, "fasta", str(bam_path)],
            stdout=tmp_fasta.open("wb"),
            stderr=subprocess.PIPE if capture else None,
            check=False,
        )
        if samtools_result.returncode != 0:
            stderr = samtools_result.stderr.decode() if samtools_result.stderr else ""
            raise ExternalToolError(
                cmd=[samtools, "fasta", str(bam_path)],
                returncode=samtools_result.returncode,
                stderr=stderr,
            )
        return run_hks_lookup(
            base_path=base_path,
            feature_set_file=feature_set_file,
            k=k,
            input_path=tmp_fasta,
            output_path=output_path,
            threads=threads,
            capture=capture,
        )
    finally:
        if tmp_fasta.exists():
            tmp_fasta.unlink()


def run_hks_smooth(
    *,
    hierarchy_file: Path,
    input_path: Path,
    output_path: Path,
    max_gap: int = _DEFAULT_SMOOTH_MAX_GAP,
    threads: int = 0,
    capture: bool = False,
) -> Path:
    """Invoke ``hks smooth`` on a lookup TSV and write the result as a BED file.

    Parameters
    ----------
    hierarchy_file
        Path to the HKS hierarchy file (``*.hierarchy.txt``) for this feature set.
    input_path
        Raw TSV produced by ``hks lookup`` (with header, ``none`` for misses).
    output_path
        Where to write the smoothed BED (headerless, ``novel`` for misses).
    max_gap
        Maximum gap in bases between adjacent intervals that are still
        considered connected for smoothing purposes.
    threads
        Number of worker threads. ``hks smooth`` parallelizes across query
        sequences. ``0`` (the default) lets HKS decide (effectively ``4``).
    capture
        If ``True``, capture subprocess stdout/stderr instead of passing through.

    Returns
    -------
    Path
        ``output_path``, after writing.

    Raises
    ------
    ToolNotFoundError
        If ``hks`` is not found.
    ExternalToolError
        If the subprocess exits with a non-zero status.
    """
    binary = get_hks_binary()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n_threads = threads if threads > 0 else 4

    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as tmp:
        smooth_tsv = Path(tmp.name)

    try:
        cmd: list[str] = [
            binary,
            "smooth",
            "--feature-hierarchy", str(hierarchy_file),
            "-i", str(input_path),
            "-o", str(smooth_tsv),
            "--max-gap", str(max_gap),
            "-t", str(n_threads),
        ]
        logger.debug("running: %s", " ".join(cmd))
        run_tool(cmd, capture=capture)
        convert_hks_tsv_to_bed(smooth_tsv, output_path)
    finally:
        smooth_tsv.unlink(missing_ok=True)

    return output_path
