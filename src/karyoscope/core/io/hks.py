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
from karyoscope.exceptions import KaryoscopeError

logger = logging.getLogger(__name__)

BINARY_NAME = "hks"
ENV_OVERRIDE = "KARYOSCOPE_HKS"

#: Label HKS emits for k-mers not found in the index, unless told otherwise.
#: KaryoScope overrides it with ``--miss-label`` so that HKS writes its own
#: :data:`~karyoscope.core.io.features.NOVEL_NAME` sentinel and no pass over
#: the output is needed to rewrite it.
_HKS_MISS_LABEL = "none"

#: Advice appended when ``hks`` is killed by an OOM-like signal. Detection
#: lives in :mod:`karyoscope.core.external`; only the figures are backend
#: specific, and HKS's are far below KMC's -- telling an HKS user to request
#: 100 GB would be as unhelpful as telling them nothing.
#:
#: The numbers are measured, not estimated: peak RSS is the index itself --
#: the shared base plus one feature set's labeling at a time -- so it is
#: ~10 GB for a human database whatever is being annotated, and barely moves
#: with --threads.
HKS_OOM_HINT = (
    "hks holds the index in memory: the shared base index plus one feature "
    "set's labeling at a time. For a human-scale database that is ~10 GB, "
    "and it does not shrink with fewer --threads or a smaller input.\n"
    "Recommended fixes:\n"
    "  * On a SLURM cluster: request at least 16 GB (e.g. --mem=16G).\n"
    "  * On a login node: move to a compute node.\n"
    "  * Check the database's index size on disk -- `du -sh <db>/index` is a "
    "good lower bound for what one lookup needs resident.\n"
)


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
    ".cram",
)

#: Alignment formats that ``hks`` cannot read and that are therefore
#: materialised to a temp FASTA via ``samtools fasta`` first. CRAM
#: additionally needs the reference its bases were encoded against.
_ALIGNMENT_EXTENSIONS: tuple[str, ...] = (".bam", ".cram")


def _relay_hks_log(result, what: str) -> None:
    """Forward a captured ``hks`` run's own log lines into ours, at DEBUG.

    ``hks`` is invoked with ``capture=True`` so a normal run does not
    interleave two tools' logs on stderr. That also means everything it
    reports about itself — how long the base index and the labeling each
    took to load, how long the query and the smoothing ran — is discarded on
    success, and only ever surfaces in the error message of a run that
    failed. Those are exactly the numbers needed to tell which phase a
    change moved, so ``-vv`` gets them.
    """
    if result is None:
        return
    for stream in (result.stderr, result.stdout):
        for line in (stream or "").splitlines():
            line = line.strip()
            if line:
                logger.debug("hks[%s]: %s", what, line)


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


def estimate_hks_memory_bytes(
    *, db_dir: Path, basename: str, feature_sets: list[str]
) -> int | None:
    """Bytes ``hks`` must hold resident to query ``feature_sets``.

    Unlike the output-size estimate, this is not extrapolated from anything:
    ``hks`` loads the shared base index plus **one** feature set's labeling
    at a time, so the peak is the base plus the largest labeling among those
    requested, and both are files whose size can simply be read.

    That structure is why the figure barely moves with ``--threads`` or with
    input size, and why a single haplotype costs the same as a diploid.
    Measured on the human database: base 6.04 GiB + largest labeling
    3.02 GiB = 9.06 GiB, against a 10.04 GB observed peak — the remainder is
    query buffering, covered by the caller's margin.

    Returns ``None`` if the files cannot be sized, so the caller can fail
    open rather than block on a database it could not inspect.
    """
    base = db_dir / f"{basename}.hksb"
    try:
        total = base.stat().st_size
    except OSError:
        return None

    largest = 0
    for fs in feature_sets:
        try:
            largest = max(largest, (db_dir / f"{basename}.{fs}.hksf").stat().st_size)
        except OSError:
            return None
    if largest == 0:
        return None
    return total + largest


def _infer_prefix(input_path: Path, db_basename: str) -> str:
    """Compute the output filename prefix for HKS lookup output."""
    input_basename = input_path.name
    for ext in _INPUT_EXTENSIONS:
        if input_basename.lower().endswith(ext):
            input_basename = input_basename[: -len(ext)]
            break
    return f"{input_basename}.{db_basename}"


def run_hks_lookup(
    *,
    base_path: Path,
    feature_set_file: Path,
    k: int,
    input_path: Path,
    output_path: Path,
    threads: int = 0,
    report_query_names: bool = True,
    reference: Path | None = None,
    capture: bool = False,
) -> Path:
    """Invoke ``hks lookup`` for one feature set against one input.

    A convenience wrapper over :func:`run_hks_lookup_batch` with a single
    ``(input, output)`` pair — deliberately not a second implementation. The
    two used to be separate, and drifted: the batch one kept asking ``hks``
    for its default output format after the single one had moved to writing
    KaryoScope's BED directly, which git merged without a conflict because
    the edits were in different functions.

    See :func:`run_hks_lookup_batch` for the parameters and for what the
    output actually is (the presmoothed BED, not a TSV awaiting conversion).
    """
    run_hks_lookup_batch(
        base_path=base_path,
        feature_set_file=feature_set_file,
        k=k,
        io_pairs=[(input_path, output_path)],
        threads=threads,
        report_query_names=report_query_names,
        reference=reference,
        capture=capture,
    )
    return output_path


def run_hks_lookup_batch(
    *,
    base_path: Path,
    feature_set_file: Path,
    k: int,
    io_pairs: list[tuple[Path, Path]],
    threads: int = 0,
    report_query_names: bool = True,
    reference: Path | None = None,
    capture: bool = False,
) -> None:
    """Invoke ``hks lookup`` ONCE for a feature set, querying many inputs.

    ``io_pairs`` is a list of ``(input_path, output_tsv_path)``. The (large) base
    index and feature-set file are loaded a single time and every input is queried
    against them in turn, each written to its paired output TSV — so annotating N
    inputs against one feature set costs one index load instead of N.

    ``report_query_names`` is a single ``hks`` flag applied to the whole batch, so
    every input here must want the same setting (the caller groups inputs by type;
    see :func:`karyoscope.core.io.hks._is_reads_input`-based grouping in the
    annotate batch backend). BAM and CRAM inputs are materialised to temp FASTA
    (via ``samtools fasta``) up front, mirroring the single-input path; the temp
    files are removed afterwards.

    ``reference`` is the FASTA the alignment was encoded against, passed to
    ``samtools fasta --reference``. CRAM stores bases as a diff against it, so
    decoding one without it fails (or silently yields wrong sequence when
    samtools resolves a *different* reference from the header's M5 tags via
    ``$REF_PATH``). It is therefore required for CRAM and ignored for BAM, which
    is self-contained.

    The goal is the *original unaligned read set*: every read's full sequence
    exactly once, as the pre-alignment FASTQ would have held it. ``-F 0x900``
    (stated explicitly rather than left to samtools' default) selects primary
    records only, which delivers exactly that:

    * Each read has exactly one primary record, and it carries the full-length
      SEQ — primaries are never hard-clipped.
    * Supplementary records are hard-clipped *slices* of reads whose full
      sequence the primary already supplied, so including them would feed those
      bases through twice.
    * Unmapped reads (``0x4``) and duplicates (``0x400``) are NOT excluded by
      this mask. They are genuine distinct reads and are kept — dropping them
      would under-count.

    ``samtools fasta`` reverse-complements minus-strand records, so what is
    written is the read as sequenced, not as aligned.

    ``-N`` forces the ``/1`` and ``/2`` mate suffix onto every read name.
    Aligners routinely strip it from QNAME (both mates of a pair then share a
    byte-identical name in the file), and samtools' default only restores it
    where the READ1/READ2 flag bits happen to be set. Without the suffix the
    two mates are indistinguishable in the output, which silently destroys any
    downstream fragment-level analysis. Measured on a DRAGEN WGS CRAM: ``-n``
    collapsed 1,996,531 records onto 1,002,139 distinct names; ``-N`` kept all
    1,996,531 distinct.

    Note these temp FASTAs are full-size — a 65x human WGS CRAM expands to
    ~260 GB — and land in ``$TMPDIR``. Point that at node-local scratch, not at
    a shared filesystem.
    """
    if not io_pairs:
        return

    binary = get_hks_binary()
    n_threads = threads if threads > 0 else 4

    tmp_fastas: list[Path] = []
    samtools: str | None = None
    try:
        # Resolve each input to a seekable query path (converting BAM up front).
        query_paths: list[Path] = []
        for input_path, output_path in io_pairs:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            suffix = input_path.suffix.lower()
            if suffix in _ALIGNMENT_EXTENSIONS:
                fmt = suffix.lstrip(".").upper()
                if suffix == ".cram" and reference is None:
                    raise KaryoscopeError(
                        f"{input_path.name} is a CRAM, which stores bases as a diff "
                        f"against the reference it was aligned to, so it cannot be "
                        f"decoded without that reference. Pass --reference "
                        f"<genome.fasta> (the same one used for alignment)."
                    )
                if samtools is None:
                    samtools = require_tool(
                        "samtools",
                        install_hint=(
                            "Install samtools to use BAM/CRAM inputs:\n"
                            "  conda install -c bioconda samtools\n"
                            "Or convert the alignment to FASTA first:\n"
                            "  samtools fasta input.bam | gzip > input.fasta.gz"
                        ),
                    )
                with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp:
                    tmp_fasta = Path(tmp.name)
                tmp_fastas.append(tmp_fasta)
                samtools_cmd = [samtools, "fasta", "-F", "0x900", "-N"]
                if reference is not None:
                    # NOT -T: in `samtools fasta` that is the copy-tags-to-header
                    # taglist, and passing a path there silently produces a
                    # tag-decorated FASTA with no reference supplied at all.
                    samtools_cmd += ["--reference", str(reference)]
                # Give the decode its own threads: on a 30x human CRAM this
                # stage is minutes of wall clock, and it is pure setup before
                # the index load can even begin.
                if threads > 0:
                    samtools_cmd += ["-@", str(threads)]
                samtools_cmd.append(str(input_path))
                logger.debug("converting %s to FASTA: %s -> %s", fmt, input_path, tmp_fasta)
                result = subprocess.run(
                    samtools_cmd,
                    stdout=tmp_fasta.open("wb"),
                    stderr=subprocess.PIPE if capture else None,
                    check=False,
                )
                if result.returncode != 0:
                    stderr = result.stderr.decode() if result.stderr else ""
                    raise ExternalToolError(
                        cmd=samtools_cmd,
                        returncode=result.returncode,
                        stderr=stderr,
                    )
                logger.info(
                    "decoded %s to %.1f GB of FASTA at %s",
                    input_path.name,
                    tmp_fasta.stat().st_size / 1_000_000_000,
                    tmp_fasta,
                )
                query_paths.append(tmp_fasta)
            else:
                query_paths.append(input_path)

        cmd: list[str] = [
            binary,
            "lookup",
            "-i",
            str(base_path),
            "--feature-set-file",
            str(feature_set_file),
            "-k",
            str(k),
        ]
        if report_query_names:
            cmd.append("--report-query-names")
        # Same output shape as the single-input path, and for the same reason:
        # these ARE the presmoothed BEDs, and run_hks_smooth parses the miss
        # token back out of them. See :func:`run_hks_lookup`.
        cmd += [
            "--report-misses",
            "--miss-label",
            NOVEL_NAME,
            "--no-header",
            "-t",
            str(n_threads),
        ]
        for query_path, (_input_path, output_path) in zip(query_paths, io_pairs, strict=True):
            cmd += ["-q", str(query_path), "-o", str(output_path)]

        logger.debug("running (batch, %d queries): %s", len(io_pairs), " ".join(cmd))
        _relay_hks_log(run_tool(cmd, capture=capture, oom_hint=HKS_OOM_HINT), "lookup-batch")
    finally:
        for p in tmp_fastas:
            if p.exists():
                p.unlink()


def run_hks_smooth(
    *,
    hierarchy_file: Path,
    input_path: Path,
    output_path: Path,
    max_gap: int = _DEFAULT_SMOOTH_MAX_GAP,
    threads: int = 0,
    capture: bool = False,
) -> Path:
    """Invoke ``hks smooth`` on a presmoothed BED, writing the smoothed one.

    ``hks smooth`` is told the output shape we want -- headerless, ``novel``
    for misses -- so its output *is* the smoothed BED. Previously it wrote a
    TSV to a temp file that a second pass then rewrote into place, which cost
    a full read plus a full write of a multi-gigabyte file per feature set,
    single threaded, for what amounts to one substitution and a dropped line.

    The miss label is not just cosmetic to smooth: it parses that same token
    out of its input to recognise a miss run. It therefore has to match what
    :func:`run_hks_lookup` wrote, which is why both pass ``NOVEL_NAME`` and
    neither hardcodes HKS's ``none`` default. The same goes for the label
    vocabulary: the input has no header to declare it, so both sides rely on
    ``--report-label-ids`` being absent, i.e. names.

    Parameters
    ----------
    hierarchy_file
        Path to the HKS hierarchy file (``*.hierarchy.txt``) for this feature set.
    input_path
        The presmoothed BED written by :func:`run_hks_lookup` (headerless,
        ``novel`` for misses).
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

    cmd: list[str] = [
        binary,
        "smooth",
        "--feature-hierarchy",
        str(hierarchy_file),
        "-i",
        str(input_path),
        "-o",
        str(output_path),
        "--max-gap",
        str(max_gap),
        "-t",
        str(n_threads),
        "--miss-label",
        NOVEL_NAME,
        "--no-header",
    ]
    logger.debug("running: %s", " ".join(cmd))
    _relay_hks_log(run_tool(cmd, capture=capture, oom_hint=HKS_OOM_HINT), "smooth")

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
    run_tool(cmd, capture=capture, oom_hint=HKS_OOM_HINT)
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
    run_tool(cmd, capture=capture, oom_hint=HKS_OOM_HINT)
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
