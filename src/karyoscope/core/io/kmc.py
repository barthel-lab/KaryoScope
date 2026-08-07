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

import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from karyoscope._version import __version__
from karyoscope.core.external import (
    ExternalToolError,
    ToolNotFoundError,
    require_tool,
    run_tool,
)
from karyoscope.exceptions import KaryoscopeError

logger = logging.getLogger(__name__)

#: The bare command name we search for on ``$PATH``.
BINARY_NAME = "get_featureIDs"

#: Filename suffix the C++ helper always appends after the user/derived
#: prefix. Centralised so the combined-BED path is computed in exactly
#: one place.
COMBINED_BED_SUFFIX = ".combined.presmoothed.featureIDs.bed"

#: Suffix of the sidecar completion marker written next to the combined
#: BED (``<combined>.done``). The marker is written only after
#: ``get_featureIDs`` exits 0, so its presence -- together with a
#: matching size/mtime -- is what lets a rerun safely *reuse* an
#: existing combined BED instead of regenerating it. A combined BED left
#: behind by a killed run has no marker (or a mismatched one) and is
#: never silently trusted. See :func:`combined_bed_is_complete`.
_MARKER_SUFFIX = ".done"

#: Bump if the marker payload schema changes incompatibly; an older or
#: unrecognised schema makes :func:`combined_bed_is_complete` fail
#: closed (regenerate) rather than trust a stale layout.
_MARKER_SCHEMA = 1

#: Environment variable users can set to override the binary location.
ENV_OVERRIDE = "KARYOSCOPE_GET_FEATUREIDS"

#: Advice appended when ``get_featureIDs`` is killed by an OOM-like signal.
#: The detection itself lives in :mod:`karyoscope.core.external`, which
#: applies to every tool; only the memory figures are backend-specific, and
#: KMC's are an order of magnitude above HKS's.
KMC_OOM_HINT = (
    "The KMC index for a human-scale database needs ~20-30 GB to load, plus "
    "per-thread working memory.\n"
    "Recommended fixes:\n"
    "  * On a SLURM cluster: request more memory (e.g. --mem=100G) and "
    "limit threads (--cpus-per-task=16).\n"
    "  * On a login node: move to a compute node.\n"
    "  * On any node, pass `karyoscope annotate -t 16 ...` explicitly. "
    "Without --threads, get_featureIDs auto-detects the machine's full "
    "core count, which on shared nodes can be much higher than the "
    "memory allocation supports.\n"
)


def _augment_with_oom_hint(
    *,
    cmd: list[str],
    returncode: int,
    stderr: str = "",
    stdout: str = "",
) -> ExternalToolError:
    """Build an :class:`ExternalToolError` carrying KMC's memory advice.

    Used on the paths that drive ``get_featureIDs`` through raw
    :mod:`subprocess` rather than :func:`run_tool`, so they get the same
    hint. Whether the advice is actually shown is decided in
    :class:`ExternalToolError` from the return code -- ordinary failures
    (malformed input, missing files) never see it, because pointing at
    memory when the cause is something else is worse than saying nothing.
    """
    return ExternalToolError(
        cmd=cmd,
        returncode=returncode,
        stderr=stderr,
        stdout=stdout,
        oom_hint=KMC_OOM_HINT,
    )


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


#: Extensions stripped from the input filename when deriving an
#: output prefix. Order matters: longer suffixes first so they win
#: over the bare-extension fallback. Includes BAM so that BAM inputs
#: (piped through ``samtools fastq``) get a clean prefix.
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
    ".cram",
)

#: Alignment formats the C++ binary cannot read, streamed in through
#: ``samtools fasta`` instead. CRAM additionally needs ``-T <reference>``.
_ALIGNMENT_EXTENSIONS: tuple[str, ...] = (".bam", ".cram")


def _infer_prefix(input_path: Path, db_path: Path) -> str:
    """Compute the prefix ``get_featureIDs`` would derive by default.

    Mirrors the C++ binary's ``get_fasta_prefix`` logic so callers can
    compute the same output filename without having to ask the binary
    to do it. Used both to reconstruct the output path after a normal
    file-input run and to pass ``--prefix`` explicitly when reading
    from stdin (the BAM pipe path).
    """
    input_basename = input_path.name
    for ext in _INPUT_EXTENSIONS:
        if input_basename.lower().endswith(ext):
            input_basename = input_basename[: -len(ext)]
            break
    kmc_basename = Path(str(db_path)).name
    # If db_path pointed at a file like "features.kmc_pre", strip the suffix.
    for suffix in (".kmc_pre", ".kmc_suf"):
        if kmc_basename.endswith(suffix):
            kmc_basename = kmc_basename[: -len(suffix)]
            break
    return f"{input_basename}.{kmc_basename}"


def combined_bed_path(output_dir: Path, prefix: str) -> Path:
    """Return the combined-BED path ``get_featureIDs`` writes for ``prefix``.

    The single source of truth for the output filename, used both to
    locate the result after a run and to test for a reusable one before
    a run (resume).
    """
    return output_dir / f"{prefix}{COMBINED_BED_SUFFIX}"


def combined_marker_path(combined_bed: Path) -> Path:
    """Path of the completion marker sidecar for ``combined_bed``."""
    return combined_bed.with_name(combined_bed.name + _MARKER_SUFFIX)


def write_combined_marker(
    combined_bed: Path,
    *,
    prefix: str,
    db_path: Path,
    input_path: Path,
) -> Path:
    """Write the ``<combined>.done`` marker recording a successful run.

    Records the combined BED's size and mtime so a later run can detect
    a file that was truncated or otherwise modified after the marker was
    written, plus provenance (input, db, version) for debugging. Call
    only after ``get_featureIDs`` has exited 0.
    """
    st = combined_bed.stat()
    payload = {
        "schema": _MARKER_SCHEMA,
        "karyoscope_version": __version__,
        "combined_bed": combined_bed.name,
        "size_bytes": st.st_size,
        "mtime_ns": st.st_mtime_ns,
        "prefix": prefix,
        "db_path": str(db_path),
        "input": str(input_path),
    }
    marker = combined_marker_path(combined_bed)
    marker.write_text(json.dumps(payload, indent=2) + "\n")
    logger.debug("wrote combined-BED completion marker %s", marker)
    return marker


def clear_combined_marker(combined_bed: Path) -> None:
    """Remove the completion marker for ``combined_bed`` if present.

    Called before a (re)run so that a crash mid-write can never leave a
    marker pointing at a half-written BED, and alongside deletion of the
    combined intermediate so no dangling marker is left behind.
    """
    combined_marker_path(combined_bed).unlink(missing_ok=True)


def combined_bed_is_complete(combined_bed: Path) -> bool:
    """True if ``combined_bed`` exists and its completion marker matches.

    Returns ``True`` only when the BED is present, its sidecar marker is
    present and parses, the schema matches, and the recorded size and
    mtime equal the file's current size and mtime. Any mismatch, missing
    file, or unreadable marker returns ``False`` -- failing closed so a
    partial or post-hoc-modified file is regenerated rather than trusted.
    """
    if not combined_bed.is_file():
        return False
    marker = combined_marker_path(combined_bed)
    try:
        data = json.loads(marker.read_text())
    except (OSError, ValueError):
        return False
    if data.get("schema") != _MARKER_SCHEMA:
        return False
    try:
        st = combined_bed.stat()
    except OSError:
        return False
    return data.get("size_bytes") == st.st_size and data.get("mtime_ns") == st.st_mtime_ns


def run_get_featureids(
    *,
    db_path: Path,
    input_path: Path,
    output_dir: Path,
    threads: int = 0,
    prefix: str | None = None,
    input_format: str | None = None,
    reference: Path | None = None,
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
        FASTA, FASTQ, or BAM input file. Plain or gzipped is fine for
        FASTA/FASTQ; ``.bam`` is detected by extension and routed
        through ``samtools fasta`` as a streaming pipe (no temp file
        and no quality strings, since k-mer counting doesn't need
        them). Pass ``Path("-")`` for stdin (and set ``input_format``).
    output_dir
        Directory where the combined BED will be written. Created if
        absent.
    threads
        Number of worker threads. ``0`` (the default) means "let
        get_featureIDs decide based on hardware concurrency."
    prefix
        Optional override for the output filename prefix. If ``None``,
        ``_infer_prefix`` derives the prefix from the input basename
        plus the KMC database basename.
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
        If the binary cannot be located (or, for BAM inputs, if
        ``samtools`` is not on ``$PATH``).
    ExternalToolError
        If the subprocess exits with a non-zero status.
    """
    binary = get_featureids_binary()
    output_dir.mkdir(parents=True, exist_ok=True)

    if prefix is None:
        prefix = _infer_prefix(input_path, db_path)

    out_path = combined_bed_path(output_dir, prefix)
    # Clear any stale completion marker up front: if this run is killed
    # mid-write, the leftover (partial) BED must NOT keep an old marker
    # that a resume would then wrongly trust. The fresh marker is
    # written only after the binary exits 0, below.
    clear_combined_marker(out_path)

    if input_path.suffix.lower() in _ALIGNMENT_EXTENSIONS:
        out_path = _run_get_featureids_piped_from_alignment(
            binary=binary,
            db_path=db_path,
            alignment_path=input_path,
            output_dir=output_dir,
            threads=threads,
            prefix=prefix,
            reference=reference,
            capture=capture,
        )
        write_combined_marker(out_path, prefix=prefix, db_path=db_path, input_path=input_path)
        return out_path

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
        "--prefix",
        prefix,
    ]
    if input_format is not None:
        cmd += ["--input-format", input_format]

    # The hint rides along on every invocation; ExternalToolError shows it
    # only when the return code says the process was killed, so ordinary
    # failures are unaffected.
    run_tool(cmd, capture=capture, oom_hint=KMC_OOM_HINT)
    write_combined_marker(out_path, prefix=prefix, db_path=db_path, input_path=input_path)
    return out_path


def _run_get_featureids_piped_from_alignment(
    *,
    binary: str,
    db_path: Path,
    alignment_path: Path,
    output_dir: Path,
    threads: int,
    prefix: str,
    reference: Path | None,
    capture: bool,
) -> Path:
    """Stream ``samtools fasta <aln>`` into ``get_featureIDs --input -``.

    BAM/CRAM inputs aren't read by the C++ binary directly; we convert
    on the fly via ``samtools fasta`` (not ``fastq`` -- KaryoScope
    only needs the sequence, not the quality string, and FASTA is
    smaller and slightly faster to write). The pipe is streaming,
    so no temp file is created and memory stays bounded. (The HKS
    backend cannot do this -- ``hks lookup`` needs a seekable path --
    which is why it materialises a temp FASTA instead.)

    ``reference`` is passed to ``samtools fasta --reference``. It is required
    for CRAM, whose bases are stored as a diff against the reference used for
    alignment, and ignored for BAM.

    The goal is the *original unaligned read set*: every read's full sequence
    exactly once, as the pre-alignment FASTQ would have held it. ``-F 0x900``
    (stated explicitly rather than left to samtools' default) selects primary
    records only, which delivers exactly that: each read has exactly one
    primary record carrying full-length SEQ, while supplementary records are
    hard-clipped slices of reads the primary already supplied in full.
    Unmapped reads (``0x4``) and duplicates (``0x400``) are NOT excluded by
    this mask -- they are genuine distinct reads and are kept.

    ``-N`` forces the ``/1``/``/2`` mate suffix onto every read name. Aligners
    routinely strip it from QNAME, leaving both mates of a pair sharing a
    byte-identical name; without the suffix they are indistinguishable in the
    output.

    Error handling: if ``get_featureIDs`` exits non-zero we report
    that (it's the more informative failure for the user). If it
    succeeds but ``samtools`` exited non-zero we report the samtools
    error (rare; usually a malformed BAM or a mismatched reference).
    """
    if alignment_path.suffix.lower() == ".cram" and reference is None:
        raise KaryoscopeError(
            f"{alignment_path.name} is a CRAM, which stores bases as a diff against "
            f"the reference it was aligned to, so it cannot be decoded without that "
            f"reference. Pass --reference <genome.fasta> (the same one used for "
            f"alignment)."
        )

    samtools = require_tool(
        "samtools",
        install_hint="Install samtools to use BAM/CRAM inputs:\n"
        "  conda install -c bioconda samtools\n"
        "Or convert the alignment to FASTA first with:\n"
        "  samtools fasta input.bam | gzip > input.fasta.gz",
    )

    samtools_cmd = [samtools, "fasta", "-F", "0x900", "-N"]
    if reference is not None:
        # NOT -T: in `samtools fasta` that is the copy-tags-to-header taglist,
        # and passing a path there silently produces a tag-decorated FASTA with
        # no reference supplied at all.
        samtools_cmd += ["--reference", str(reference)]
    samtools_cmd.append(str(alignment_path))
    getfid_cmd = [
        binary,
        "--db",
        str(db_path),
        "--input",
        "-",
        "--input-format",
        "fasta",
        "--output",
        str(output_dir),
        "--threads",
        str(threads),
        "--prefix",
        prefix,
    ]
    logger.debug("piping: %s | %s", " ".join(samtools_cmd), " ".join(getfid_cmd))

    # samtools stderr goes to a temp file, not a pipe: nothing drains a
    # pipe until get_featureIDs finishes, so a warning-heavy BAM (one
    # line per malformed record) could fill the ~64 KB pipe buffer,
    # block samtools mid-write, stall its stdout, and deadlock the
    # whole pipeline. A file has no such limit.
    with tempfile.TemporaryFile(mode="w+", encoding="utf-8", errors="replace") as errf:
        samtools_proc = subprocess.Popen(
            samtools_cmd,
            stdout=subprocess.PIPE,
            stderr=errf,
            text=True,
        )
        try:
            getfid_proc = subprocess.run(
                getfid_cmd,
                stdin=samtools_proc.stdout,
                capture_output=True,
                text=True,
                check=False,
            )
        finally:
            # Closing the read end before wait() lets samtools see SIGPIPE
            # cleanly if get_featureIDs exited early.
            if samtools_proc.stdout is not None:
                samtools_proc.stdout.close()
        samtools_returncode = samtools_proc.wait()
        errf.seek(0)
        samtools_stderr = errf.read()

    # Order matters: downstream failure (get_featureIDs) is the more
    # actionable error to surface first. Both paths route through
    # _augment_with_oom_hint, which appends the OOM hint only for
    # SIGKILL-like exit codes (plain failures pass through unchanged).
    if getfid_proc.returncode != 0:
        raise _augment_with_oom_hint(
            cmd=getfid_cmd,
            returncode=getfid_proc.returncode,
            stderr=getfid_proc.stderr or "",
            stdout=getfid_proc.stdout or "",
        )
    if samtools_returncode != 0:
        raise _augment_with_oom_hint(
            cmd=samtools_cmd,
            returncode=samtools_returncode,
            stderr=samtools_stderr,
        )

    if not capture:
        # Forward subprocess output to our own stdout/stderr so the
        # user sees it just like with the file-input path.
        if getfid_proc.stdout:
            sys.stdout.write(getfid_proc.stdout)
        if getfid_proc.stderr:
            sys.stderr.write(getfid_proc.stderr)
        if samtools_stderr:
            sys.stderr.write(samtools_stderr)

    return combined_bed_path(output_dir, prefix)
