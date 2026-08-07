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


#: Bytes read per block by :func:`convert_hks_tsv_to_bed`. Large enough
#: that per-block Python overhead vanishes against the memchr-speed scan
#: inside ``bytes.replace``, small enough that the buffer is irrelevant
#: beside the k-mer index this runs alongside. Exposed as a parameter only
#: so tests can shrink it and hammer the block-boundary logic.
_CONVERT_BLOCK_BYTES = 8 << 20  # 8 MiB


def convert_hks_tsv_to_bed(
    tsv_path: Path, bed_path: Path, *, block_bytes: int = _CONVERT_BLOCK_BYTES
) -> None:
    """Convert HKS TSV output (with header) to a headerless BED file.

    Strips the header line and replaces the HKS miss label ``none`` with
    KaryoScope's ``novel`` sentinel.

    Runs twice per feature set — once on the raw lookup TSV to make the
    presmoothed BED, once inside :func:`run_hks_smooth` on the smoothed
    TSV — and on human-scale input those two passes were ~10% of
    ``annotate``'s wall time (129 s of a 21-minute HG002 run), single
    threaded, while the rest of the machine idled. The previous version
    looped in Python over every line of a multi-GB file and paid a decode
    plus an encode per line; this one reads binary blocks and lets
    ``bytes.replace`` do the scan in C.

    Still streams: memory is one block plus at most one partial line, so a
    reads TSV of tens of GB is fine.

    Correctness at block boundaries
    ------------------------------
    Only whole lines are ever handed to ``replace``: each block is cut at
    its **last newline** and the remainder carried into the next block.
    That makes a straddling match impossible, because ``miss`` ends in a
    newline — an occurrence extending past the cut would have to end at a
    newline beyond the last one, which by definition does not exist.
    ``miss`` also contains no interior newline, so it can never match
    across two lines. The carried remainder is a partial final line and
    cannot contain a complete match.

    Byte-for-byte identical output to the line-by-line version on any
    input, including a final line with no trailing newline (whose label,
    lacking the terminating newline, is left alone by both).
    """
    miss = f"\t{_HKS_MISS_LABEL}\n".encode()
    novel = f"\t{NOVEL_NAME}\n".encode()

    with tsv_path.open("rb") as inp, bed_path.open("wb") as out:
        # Drop the header, which may in principle span more than one read.
        pending = b""
        while True:
            block = inp.read(block_bytes)
            if not block:
                return  # empty or header-only: no records to write
            pending += block
            newline = pending.find(b"\n")
            if newline != -1:
                pending = pending[newline + 1 :]
                break

        while True:
            block = inp.read(block_bytes)
            if not block:
                break
            pending += block
            cut = pending.rfind(b"\n")
            if cut == -1:
                continue  # a line longer than one block; keep accumulating
            out.write(pending[: cut + 1].replace(miss, novel))
            pending = pending[cut + 1 :]

        if pending:
            out.write(pending.replace(miss, novel))


def run_hks_lookup(
    *,
    base_path: Path,
    feature_set_file: Path,
    k: int,
    input_path: Path,
    output_path: Path,
    threads: int = 0,
    report_query_names: bool = True,
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
    report_query_names
        If ``True`` (default), pass ``--report-query-names`` so column 1 of the
        output holds each sequence's name -- correct for assemblies, whose
        contig names map to karyotype chromosomes. Pass ``False`` for reads:
        HKS then emits integer query ranks and skips its pre-pass that loads
        every sequence name into memory (~10 GB at hundreds of millions of
        reads). The caller decides based on input type, because a BAM input is
        materialised to a temp ``.fasta`` and cannot be reclassified here.
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
            report_query_names=report_query_names,
            capture=capture,
        )

    n_threads = threads if threads > 0 else 4
    cmd: list[str] = [
        binary,
        "lookup",
        "-i",
        str(base_path),
        "--feature-set-file",
        str(feature_set_file),
        "-k",
        str(k),
        "-q",
        str(input_path),
    ]
    if report_query_names:
        cmd.append("--report-query-names")
    cmd += [
        "--report-misses",
        "-t",
        str(n_threads),
        "-o",
        str(output_path),
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
    report_query_names: bool = True,
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

    # Next to the output, not the system tempdir: the FASTA is input-sized
    # (whole genome or read set), and /tmp on cluster nodes is often small
    # or RAM-backed.
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, dir=output_path.parent) as tmp:
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
            report_query_names=report_query_names,
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

    # Next to the output, not the system tempdir: the smoothed TSV is the
    # size of the lookup TSV (GBs for an assembly, tens of GBs for reads),
    # and /tmp on cluster nodes is often small or RAM-backed.
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False, dir=output_path.parent) as tmp:
        smooth_tsv = Path(tmp.name)

    try:
        cmd: list[str] = [
            binary,
            "smooth",
            "--feature-hierarchy",
            str(hierarchy_file),
            "-i",
            str(input_path),
            "-o",
            str(smooth_tsv),
            "--max-gap",
            str(max_gap),
            "-t",
            str(n_threads),
        ]
        logger.debug("running: %s", " ".join(cmd))
        run_tool(cmd, capture=capture)
        convert_hks_tsv_to_bed(smooth_tsv, output_path)
    finally:
        smooth_tsv.unlink(missing_ok=True)

    return output_path


# --- index construction ----------------------------------------------
#
# The two phases below build a database rather than query one. ``build-base``
# creates the shared s-mer index (``.hksb``) once; ``add-feature-set`` layers a
# named labeling (``.hksf``) on top of it, once per feature set. See the HKS
# README "Build an index" section and :mod:`karyoscope.core.build`.


def run_hks_build_base(
    *,
    output_path: Path,
    s: int,
    input_path: Path | None = None,
    input_file_list: Path | None = None,
    threads: int = 0,
    mem_gigas: int = 8,
    external_memory: Path | None = None,
    forward_only: bool = False,
    capture: bool = False,
) -> Path:
    """Invoke ``hks build-base`` to build the shared base index (``*.hksb``).

    Exactly one of ``input_path`` or ``input_file_list`` must be given: a single
    FASTA/FASTQ, or a file listing one input path per line (both strands are
    added unless ``forward_only``).

    Parameters
    ----------
    output_path
        Where to write the base index. Recommended extension: ``.hksb``.
    s
        Maximum query length (``-s``), up to 256. Every feature set built on
        this base can later be queried at any ``k <= s``.
    threads
        Worker threads. ``0`` (default) lets HKS decide (effectively ``4``).
    mem_gigas
        RAM budget for SBWT construction, in gigabytes (``--mem-gigas``).
    external_memory
        If given, run in external-memory mode using this directory as scratch
        (lower RAM peak, slower; identical output).
    forward_only
        Do not add reverse-complemented k-mers.

    Returns
    -------
    Path
        ``output_path``, after writing.

    Raises
    ------
    ValueError
        If neither or both of ``input_path`` / ``input_file_list`` are given.
    ToolNotFoundError
        If ``hks`` is not found.
    ExternalToolError
        If the subprocess exits with a non-zero status.
    """
    if (input_path is None) == (input_file_list is None):
        raise ValueError("run_hks_build_base requires exactly one of input_path or input_file_list")

    binary = get_hks_binary()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    n_threads = threads if threads > 0 else 4

    cmd: list[str] = [binary, "build-base", "-s", str(s), "-o", str(output_path)]
    if input_path is not None:
        cmd += ["--input", str(input_path)]
    else:
        cmd += ["--input-file-list", str(input_file_list)]
    cmd += ["-t", str(n_threads), "--mem-gigas", str(mem_gigas)]
    if external_memory is not None:
        cmd += ["--external-memory", str(external_memory)]
    if forward_only:
        cmd.append("--forward-only")

    logger.debug("running: %s", " ".join(cmd))
    run_tool(cmd, capture=capture)
    return output_path


def run_hks_add_feature_set(
    *,
    base_path: Path,
    output_path: Path,
    feature_set_name: str,
    feature_names: Path,
    feature_hierarchy: Path,
    feature_file_list: Path | None = None,
    feature_per_seq_file: Path | None = None,
    feature_priorities: Path | None = None,
    variable_k_support: bool = False,
    forward_only: bool = False,
    threads: int = 0,
    capture: bool = False,
) -> Path:
    """Invoke ``hks add-feature-set`` to build one feature-set file (``*.hksf``).

    Exactly one of ``feature_file_list`` (one FASTA per feature) or
    ``feature_per_seq_file`` (a single FASTA, one sequence per feature) must be
    given. All k-mers in those inputs must already be present in ``base_path``.

    ``variable_k_support`` and ``feature_priorities`` are mutually exclusive in
    HKS; passing both raises before invoking the tool.

    Parameters
    ----------
    base_path
        Existing base index (``*.hksb``).
    output_path
        Where to write the feature-set file (``*.hksf``).
    feature_set_name
        Name recorded for this labeling.
    feature_names
        File with one feature name per line, in input order.
    feature_hierarchy
        Edge-list hierarchy file (``child parent`` per line) for this set.
    feature_priorities
        Optional ``<name> <priority>`` file enabling priority-aware LCA. Lower
        value = higher priority; nodes absent from the file default to ``0``.
    variable_k_support
        Enable queries at any ``k <= s``. Cannot be combined with
        ``feature_priorities``.
    forward_only
        Do not add reverse-complemented k-mers.

    Returns
    -------
    Path
        ``output_path``, after writing.

    Raises
    ------
    ValueError
        On invalid input combinations.
    ToolNotFoundError
        If ``hks`` is not found.
    ExternalToolError
        If the subprocess exits with a non-zero status.
    """
    if (feature_file_list is None) == (feature_per_seq_file is None):
        raise ValueError(
            "run_hks_add_feature_set requires exactly one of feature_file_list "
            "or feature_per_seq_file"
        )
    if variable_k_support and feature_priorities is not None:
        raise ValueError("variable_k_support and feature_priorities are mutually exclusive in HKS")

    binary = get_hks_binary()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    n_threads = threads if threads > 0 else 4

    cmd: list[str] = [
        binary,
        "add-feature-set",
        "-i",
        str(base_path),
        "-o",
        str(output_path),
        "--feature-set-name",
        feature_set_name,
        "--feature-names",
        str(feature_names),
        "--feature-hierarchy",
        str(feature_hierarchy),
    ]
    if feature_file_list is not None:
        cmd += ["--feature-file-list", str(feature_file_list)]
    else:
        cmd += ["--feature-per-seq-file", str(feature_per_seq_file)]
    if feature_priorities is not None:
        cmd += ["--feature-priorities", str(feature_priorities)]
    if variable_k_support:
        cmd.append("--variable-k-support")
    if forward_only:
        cmd.append("--forward-only")
    cmd += ["-t", str(n_threads)]

    logger.debug("running: %s", " ".join(cmd))
    run_tool(cmd, capture=capture)
    return output_path


def validate_sibling_priorities(
    parent_of: dict[str, str],
    priority: dict[str, int],
) -> list[str]:
    """Check the HKS constraint that each sibling group is priority-consistent.

    HKS's priority-aware LCA (``PriorityLca::new`` in ``priority_lca.rs``)
    rejects any group of siblings whose priorities are neither all equal nor all
    distinct, because a mix makes ``plca`` non-associative. We validate the same
    rule up front so a misauthored priority file fails with a clear message here
    rather than an opaque crash inside ``hks add-feature-set``.

    Parameters
    ----------
    parent_of
        ``{child: parent}`` for one feature set (the root has no entry).
    priority
        ``{node: priority}`` for every node; callers should default missing
        nodes to ``0`` before calling, matching HKS.

    Returns
    -------
    list[str]
        Human-readable issue strings; empty when every sibling group is valid.
    """
    children_of: dict[str, list[str]] = {}
    for child, parent in parent_of.items():
        children_of.setdefault(parent, []).append(child)

    issues: list[str] = []
    for parent, kids in children_of.items():
        if len(kids) < 2:
            continue
        prios = [priority.get(k, 0) for k in kids]
        if len(set(prios)) == 1:
            continue  # all equal — ok
        if len(set(prios)) == len(prios):
            continue  # all distinct — ok
        # Mixed: report the priorities that are shared by more than one sibling.
        by_prio: dict[int, list[str]] = {}
        for k, p in zip(kids, prios, strict=True):
            by_prio.setdefault(p, []).append(k)
        dups = {p: ks for p, ks in by_prio.items() if len(ks) > 1}
        issues.append(
            f"children of {parent!r} have mixed priorities (must be all-equal or "
            f"all-distinct): duplicates {dups!r}. Assign distinct priorities or "
            f"group the tied siblings under a named parent node."
        )
    return issues
