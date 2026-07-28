"""Orchestration for the ``annotate`` command.

The high-level flow:

1. Locate the database and validate its layout.
2. Invoke ``get_featureIDs`` (the vendored C++ helper) against the input
   FASTA. The result is a single combined BED whose 4th column holds
   integer feature ids.
3. Read the database's ``features.tsv`` (and ``hierarchy.tsv`` when
   smoothing is enabled) to learn the mapping from feature ids to
   per-feature-set feature names and the relationships between them.
4. For each requested feature set, stream the combined BED through a
   worker pool that translates feature ids to names, produces the
   *presmoothed* BED (adjacent same-name intervals merged), and
   optionally also runs the hierarchy-aware smoothing pass to produce
   the *smoothed* BED.
5. Optionally bgzip the outputs (default on) and delete files the user
   doesn't want (the combined intermediate, the presmoothed BEDs when
   ``--no-keep-presmoothed`` was passed).

Smoothing is implemented in :mod:`karyoscope.core.smooth` and runs in a
``multiprocessing.Pool``. The number of worker processes is taken from
the same ``--threads`` argument that controls the C++ k-mer query
threads.
"""

from __future__ import annotations

import contextlib
import logging
import multiprocessing as mp
import os
import shutil
import sys
import threading
import time
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

from karyoscope import cpus as _cpus
from karyoscope import diskspace, preflight
from karyoscope import installed as _installed
from karyoscope.core.io.bgzip import bgzip_file
from karyoscope.core.io.features import NOVEL_NAME, Features, parse_features, render_feature
from karyoscope.core.io.hierarchy import (
    Hierarchy,
    parse_hierarchy,
    validate_hierarchy,
)
from karyoscope.core.io.hks import (
    convert_bam_to_fasta,
    run_hks_lookup,
    run_hks_lookup_batch,
    run_hks_smooth,
)
from karyoscope.core.io.kmc import (
    clear_combined_marker,
    combined_bed_is_complete,
    combined_bed_path,
    run_get_featureids,
)
from karyoscope.core.smooth import (
    HierarchyIndex,
    chunked_seq_reader,
    make_features_for_worker,
    process_seq_chunk,
    safe_filename,
    worker_initializer,
)
from karyoscope.exceptions import (
    DatabaseNotFoundError,
    KaryoscopeError,
)
from karyoscope.manifest import validate_database_layout
from karyoscope.progress import SILENT, Progress

logger = logging.getLogger(__name__)


#: Input extensions we recognise when deriving the output basename.
#: Order matters — longer extensions first so they win over shorter
#: ones. Includes FASTQ (read by ``get_featureIDs`` directly) and
#: BAM (piped through ``samtools fastq``); see
#: :func:`karyoscope.core.io.kmc.run_get_featureids` for the BAM
#: streaming path.
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


# --- result type ------------------------------------------------------


@dataclass
class AnnotateResult:
    """What :func:`annotate` produces.

    Attributes
    ----------
    presmoothed_paths
        ``{feature_set: path}`` for the per-feature-set presmoothed BEDs.
        Empty when the user passed ``--no-keep-presmoothed``. Path
        extensions are ``.bed.gz`` when bgzipped, ``.bed`` otherwise.
    smoothed_paths
        ``{feature_set: path}`` for the per-feature-set smoothed BEDs.
        Empty when the user passed ``--no-smooth``.
    combined_intermediate
        Path to the C++ combined BED, or ``None`` if it was deleted.
    """

    presmoothed_paths: dict[str, Path] = field(default_factory=dict)
    smoothed_paths: dict[str, Path] = field(default_factory=dict)
    combined_intermediate: Path | None = None

    @property
    def all_output_paths(self) -> list[Path]:
        """All output BED paths in a flat list (for CLI display)."""
        return list(self.presmoothed_paths.values()) + list(self.smoothed_paths.values())


# --- helpers ----------------------------------------------------------


def _resolve_query_k(manifest, k: int | None, db_id: str) -> int:
    """Return the k-mer length to query with, validating an explicit override.

    Without ``k``, uses ``manifest.kmer.size``. An explicit ``k`` is honoured
    only on a variable-k HKS index (``kmer.type == "variable"``), where any
    ``1 <= k <= kmer.max_size`` is valid. On a fixed-k index the only queryable
    length is ``kmer.size``, so any other ``k`` is a hard error pointing at the
    variable-k build option.
    """
    if k is None:
        return manifest.kmer.size
    if k < 1:
        raise KaryoscopeError(f"--k must be >= 1, got {k}")
    if k == manifest.kmer.size:
        return k
    if manifest.kmer.type != "variable":
        raise KaryoscopeError(
            f"database {db_id!r} is a fixed-k index (kmer.type={manifest.kmer.type!r}); "
            f"it can only be queried at k={manifest.kmer.size}. Build a variable-k index "
            f"(`karyoscope build --variable-k`) to query other k values."
        )
    if k > manifest.kmer.max_size:
        raise KaryoscopeError(
            f"--k {k} exceeds this index's maximum queryable k "
            f"(kmer.max_size={manifest.kmer.max_size})."
        )
    return k


def _derive_input_basename(input_path: Path) -> str:
    """Strip recognised FASTA/FASTQ/BAM extensions from a path's filename.

    ``my_assembly.fa.gz`` -> ``my_assembly``,
    ``reads.fastq.gz`` -> ``reads``,
    ``aln.bam`` -> ``aln``.
    Falls back to the raw stem if no known extension matches.
    """
    name = input_path.name
    name_lower = name.lower()
    for ext in _INPUT_EXTENSIONS:
        if name_lower.endswith(ext):
            return name[: -len(ext)]
    return input_path.stem


#: Read-level input extensions. Used to pick the right smoothing
#: dispatch path: read-level inputs (FASTQ, BAM) with preserve-order
#: get the streaming-ordered path (no per-sequence temp files), while
#: assemblies get the per-sequence temp-files path. See
#: :func:`_smooth_all_feature_sets` for the dispatch.
_READS_INPUT_EXTENSIONS: tuple[str, ...] = (
    ".fastq.gz",
    ".fq.gz",
    ".fastq",
    ".fq",
    ".bam",
)


def _is_reads_input(input_path: Path) -> bool:
    """True if ``input_path`` looks like read-level data (FASTQ or BAM).

    Pure extension check. Long-read FASTA files (millions of reads with
    a ``.fasta`` extension) are technically read-level data but will
    return ``False`` here -- users with that case should pass
    ``--no-preserve-order`` to opt out of the per-sequence temp files
    explicitly. We don't try to detect by content because pre-scanning
    a multi-GB input would be slow.
    """
    name_lower = input_path.name.lower()
    return any(name_lower.endswith(ext) for ext in _READS_INPUT_EXTENSIONS)


# --- output-size estimation ------------------------------------------
#
# Annotate is the command most likely to fill a disk: a six-feature-set
# run on a diploid human assembly writes ~29 GB of uncompressed BED, and
# bgzip does not reduce the peak because it runs after every BED has been
# written. The constants below let us say so before the run starts,
# instead of after twenty minutes of work.
#
# Calibration (2026-07, HKS_human_CHM13_v2, HG002 v1.1, 6.00 Gbp diploid,
# 6 feature sets, --no-bgzip):
#
#     presmoothed  21.70 GB  ->  0.603 bytes per base per feature set
#     smoothed      7.26 GB  ->  0.202 bytes per base per feature set
#
# These describe BED content, so they hold for both backends — the
# smoothed/presmoothed output of a KMC database is the same text as an
# HKS one.

#: Bytes of presmoothed BED per input base, per feature set.
PRESMOOTHED_BYTES_PER_BASE = 0.60

#: Bytes of smoothed BED per input base, per feature set.
SMOOTHED_BYTES_PER_BASE = 0.20

#: Transient intermediate allowance, as a multiple of one feature set's
#: presmoothed estimate. Covers the KMC combined BED, and the HKS lookup
#: output in the one case where it is still a temp file. Both are "one row
#: per k-mer run" files dominated by the sequence name and coordinates that
#: a presmoothed BED also carries, so one feature set's worth is the right
#: scale; the 1.5 absorbs the spread between feature sets (the largest set
#: measured 1.48x the six-set average).
#:
#: The HKS backend usually needs none of it. ``hks lookup`` writes the
#: presmoothed BED directly, so when that output is being kept there is no
#: second copy on disk at any point — the allowance only applies under
#: ``--no-keep-presmoothed``, where the same file becomes a temp one that
#: ``hks smooth`` reads and we then delete.
TRANSIENT_INTERMEDIATE_FACTOR = 1.5

#: Uncompressed-to-compressed ratio assumed for a gzipped nucleotide
#: FASTA/FASTQ when no ``.fai`` is available. Measured 3.45x on HG002 v1.1
#: (1.77 GB gz -> 6.10 GB plain); rounded down so the estimate errs low
#: rather than blocking runs that would have fit.
GZIP_EXPANSION_FACTOR = 3.4

#: Fraction of an uncompressed file's bytes that are sequence, by format.
#: FASTA loses ~1.7% to headers and line breaks; FASTQ carries a quality
#: string per base plus two header lines, so barely half its bytes are
#: sequence. BAM stores bases 4-bit-packed alongside compressed qualities,
#: which nets out near 1 base per byte of file.
_SEQUENCE_FRACTION_FASTA = 0.98
_SEQUENCE_FRACTION_FASTQ = 0.49
_BASES_PER_BAM_BYTE = 1.0


def _bases_from_fai(input_path: Path) -> int | None:
    """Total sequence length from a samtools ``.fai`` index, if one exists.

    Exact, and free — assemblies distributed as bgzipped FASTA usually
    ship the index alongside. Returns None when there is no usable index,
    leaving the caller to fall back on file-size heuristics.
    """
    fai = Path(str(input_path) + ".fai")
    if not fai.is_file():
        return None
    total = 0
    try:
        with fai.open() as handle:
            for line in handle:
                fields = line.split("\t")
                if len(fields) < 2:
                    return None
                total += int(fields[1])
    except (OSError, ValueError) as exc:
        logger.debug("could not read %s: %s", fai, exc)
        return None
    return total or None


def estimate_input_bases(input_path: Path) -> int:
    """Estimate the number of sequence bases in ``input_path``.

    Prefers an exact count from a sibling ``.fai``; otherwise scales the
    file size by the format's sequence fraction, expanding first if the
    file is gzipped. The result feeds a disk-space estimate, so being
    within a factor of ~1.2 is entirely adequate.
    """
    exact = _bases_from_fai(input_path)
    if exact is not None:
        logger.debug("input size from %s.fai: %d bases", input_path.name, exact)
        return exact

    try:
        file_size = input_path.stat().st_size
    except OSError:
        return 0

    name = input_path.name.lower()
    if name.endswith(".bam"):
        return int(file_size * _BASES_PER_BAM_BYTE)

    uncompressed = file_size * GZIP_EXPANSION_FACTOR if name.endswith(".gz") else file_size
    is_fastq = any(name.endswith(ext) for ext in (".fastq", ".fq", ".fastq.gz", ".fq.gz"))
    fraction = _SEQUENCE_FRACTION_FASTQ if is_fastq else _SEQUENCE_FRACTION_FASTA
    return int(uncompressed * fraction)


def estimate_output_bytes(
    *,
    input_bases: int,
    n_feature_sets: int,
    keep_presmoothed: bool,
    smooth: bool,
    index_type: str = "kmc",
) -> int:
    """Estimate peak bytes written to the output directory by one annotate run.

    "Peak", not "final": ``--bgzip`` compresses each BED only after all of
    them have been written, so it shrinks the result but not the high-water
    mark. ``--no-keep-presmoothed`` and ``--no-smooth`` do reduce it, and
    are reflected here.

    ``index_type`` decides whether a transient intermediate is counted at
    all. An HKS run keeping its presmoothed output writes no second copy of
    anything, so charging it for one would refuse runs that fit — about
    5 GB's worth on a diploid human assembly. The default is the
    conservative backend, so a caller that does not know errs high.
    """
    per_set = 0.0
    if keep_presmoothed:
        per_set += PRESMOOTHED_BYTES_PER_BASE
    if smooth:
        per_set += SMOOTHED_BYTES_PER_BASE
    outputs = input_bases * per_set * n_feature_sets

    needs_intermediate = index_type != "hks" or not keep_presmoothed
    transient = (
        input_bases * PRESMOOTHED_BYTES_PER_BASE * TRANSIENT_INTERMEDIATE_FACTOR
        if needs_intermediate
        else 0.0
    )
    return int(outputs + transient)


def _annotate_dependencies(*, index_type: str, input_path: Path, bgzip: bool) -> list[str]:
    """External tools this particular annotate run will need.

    Resolved from the database's backend, the input format, and the output
    flags rather than assumed, so a user is never told to install ``kmc``
    for an HKS run or ``samtools`` for a FASTA one.
    """
    needed = ["hks"] if index_type == "hks" else ["get_featureIDs"]
    if input_path.suffix.lower() == ".bam":
        needed.append("samtools")
    if bgzip:
        needed.append("bgzip")
    return needed


def resolve_database(
    db_root: Path,
    db_id: str | None,
) -> tuple[str, Path]:
    """Locate and validate the database to use.

    Returns ``(database_id, database_directory)``. Raises a clean error
    if the user requested an unknown id or if the layout is invalid.

    Shared across :mod:`karyoscope.commands.annotate`,
    :mod:`karyoscope.commands.bin_cmd`, and (eventually)
    :mod:`karyoscope.commands.scaffold`.
    """
    state = _installed.load(db_root)

    if db_id is not None:
        record = state.databases.get(db_id)
        if record is None:
            installed_ids = sorted(state.databases.keys())
            raise DatabaseNotFoundError(
                f"database {db_id!r} is not installed at {db_root}. "
                f"Installed databases: {installed_ids or '(none)'}. "
                "Run `karyoscope download --list` to see what's available."
            )
        db_dir = db_root / record.directory
    else:
        # No explicit id — use the only installed db if there's exactly
        # one, otherwise complain.
        if not state.databases:
            raise DatabaseNotFoundError(
                f"no databases installed at {db_root}. "
                "Install one with `karyoscope download`, or point --db-root at a "
                "directory that already has one (then --db <ID> selects it)."
            )
        if len(state.databases) > 1:
            ids = sorted(state.databases.keys())
            raise DatabaseNotFoundError(
                f"multiple databases installed; pass --db to choose one. Installed: {ids}"
            )
        db_id, record = next(iter(state.databases.items()))
        db_dir = db_root / record.directory

    if not db_dir.is_dir():
        raise DatabaseNotFoundError(
            f"database {db_id!r} is recorded but its directory {db_dir} is missing on disk."
        )

    # Surface layout errors early rather than letting get_featureIDs fail
    # with a less informative message later.
    validate_database_layout(db_dir)

    return db_id, db_dir


def _iter_bed_records(handle: TextIO) -> Iterator[tuple[str, int, int, int]]:
    """Yield ``(seq_name, start, end, feature_id)`` from a 4-column BED stream.

    Skips blank lines and lines starting with ``#``. Raises ``ValueError``
    on malformed numeric fields — the caller should reframe as a clean
    user error.
    """
    for line_no, raw in enumerate(handle, start=1):
        line = raw.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            raise ValueError(
                f"line {line_no}: expected at least 4 tab-separated columns, "
                f"got {len(parts)}: {raw!r}"
            )
        try:
            start = int(parts[1])
            end = int(parts[2])
            feature_id = int(parts[3])
        except ValueError as e:
            raise ValueError(
                f"line {line_no}: could not parse integers in BED record: {raw!r}"
            ) from e
        yield parts[0], start, end, feature_id


def _split_combined_bed(
    combined_bed: Path,
    feature_sets: list[str],
    features: Features,
    output_paths: dict[str, Path],
) -> None:
    """Read ``combined_bed`` and write one per-feature-set BED per requested set.

    Adjacent records on the same sequence whose translated feature name
    matches are merged on the fly (the C++ already does this for
    same-feature-id runs, but translation across feature sets can collapse
    distinct ids into the same name — e.g., two adjacent rows mapping to
    chr1 via different region ids).
    """
    # Per-feature-set running record. Stored as 4-tuples
    # (seq_name, start, end, current_feature_name) or None.
    pending: dict[str, tuple[str, int, int, str] | None] = {fs: None for fs in feature_sets}
    handles: dict[str, TextIO] = {fs: output_paths[fs].open("w") for fs in feature_sets}

    # Flat per-set {id: name} tables, seeded with the id-0 "novel"
    # sentinel: this loop runs once per record per feature set over a
    # combined BED with hundreds of millions of records, so the per-hit
    # cost must be one dict lookup, not render_feature's validation and
    # nested lookups. render_feature stays as the miss path, so an
    # unknown id raises the same FeaturesError as before.
    lookups = [(fs, {0: NOVEL_NAME, **features.names_for_set(fs)}) for fs in feature_sets]

    def flush(fs: str) -> None:
        rec = pending[fs]
        if rec is None:
            return
        seq, start, end, name = rec
        handles[fs].write(f"{seq}\t{start}\t{end}\t{name}\n")
        pending[fs] = None

    try:
        with combined_bed.open() as f:
            for seq_name, start, end, fid in _iter_bed_records(f):
                for fs, id_to_name in lookups:
                    name = id_to_name.get(fid)
                    if name is None:
                        name = render_feature(fid, fs, features)
                    rec = pending[fs]
                    if (
                        rec is not None
                        and rec[0] == seq_name
                        and rec[3] == name
                        and rec[2] == start
                    ):
                        # Extend the running interval.
                        pending[fs] = (rec[0], rec[1], end, rec[3])
                    else:
                        flush(fs)
                        pending[fs] = (seq_name, start, end, name)

        for fs in feature_sets:
            flush(fs)
    finally:
        for h in handles.values():
            h.close()


# --- dead-worker watchdog -------------------------------------------
#
# Background: ``multiprocessing.Pool`` does NOT reassign tasks whose
# worker was killed by an external signal (most commonly SIGKILL from
# the kernel OOM-killer on a memory-constrained node). The main thread
# blocks forever in ``pool.imap_unordered`` waiting for results that
# will never arrive — observed in practice as multi-hour silent hangs
# on whole-genome HG002 runs against the 50 GB SLURM allocation.
#
# This watchdog runs in a daemon thread, polls the pool's worker list
# every couple of seconds, and on detecting a dead or pool-replaced
# worker, writes an actionable stderr message and calls os._exit(137).
# os._exit bypasses Python cleanup (atexit, TemporaryDirectory) —
# acceptable because the user is already in a failure state, and the
# alternative (signalling the main thread while it's blocked deep in
# C code) is unreliable.


def _detect_worker_death(
    pool: mp.pool.Pool,
    initial_pids: set[int],
) -> tuple[list[int], list[int]] | None:
    """Inspect a Pool for worker death; return diagnostic info or ``None``.

    Two failure modes are detected:

    * A current worker has a non-``None`` exitcode -- it died and the
      pool's ``_handle_workers`` thread hasn't replaced it yet.
    * A current worker's PID is not in ``initial_pids`` -- the pool
      already replaced a dead worker, so the new PID is evidence the
      original died.

    Returns ``(died_pids, exitcodes)`` on detection (either list may
    be empty if all dead workers had already been cleaned up by the
    time we looked), or ``None`` when all workers are accounted for.

    Pure function. Reads ``pool._pool`` (documented-private since
    Python 3.4) and each worker's ``pid`` / ``exitcode`` attributes
    -- no thread state, safe to call from anywhere.
    """
    current = pool._pool
    current_pids = {w.pid for w in current if w.pid is not None}
    dead = [w for w in current if w.exitcode is not None]
    new_pids = current_pids - initial_pids
    if not dead and not new_pids:
        return None
    died = sorted({w.pid for w in dead if w.pid is not None})
    if not died:
        died = sorted(initial_pids - current_pids)
    exitcodes = sorted({w.exitcode for w in dead if w.exitcode is not None})
    return died, exitcodes


def _spawn_pool_watchdog(
    pool: mp.pool.Pool,
    *,
    check_interval_s: float = 2.0,
) -> threading.Event:
    """Start a daemon thread that ``os._exit``s if a pool worker dies.

    Returns a ``threading.Event``; the caller must ``set()`` it before
    the pool is legitimately closed (in a ``finally`` block immediately
    around the ``imap_unordered`` loop), so the watchdog stops checking
    before pool teardown intentionally exits the workers and triggers
    false positives.
    """
    stop_event = threading.Event()
    initial_pids = {w.pid for w in pool._pool if w.pid is not None}

    def _watch() -> None:
        while not stop_event.wait(check_interval_s):
            detection = _detect_worker_death(pool, initial_pids)
            if detection is None or stop_event.is_set():
                continue
            died, exitcodes = detection
            msg = (
                "\nFATAL: smoothing worker process(es) died unexpectedly "
                f"(pid(s)={died}, exitcode(s)={exitcodes}).\n"
                "  Most likely cause: the kernel OOM-killer reaped a "
                "worker under memory pressure.\n"
                "  Each smoothing worker holds the full per-sequence "
                "interval list in memory (~0.5 GB per million intervals; "
                "the longest contigs of a human assembly are 5-6 GB each).\n"
                "  Resolve by reducing --threads (e.g. -t 8 for typical "
                "human inputs on 50 GB nodes) or increasing the job's "
                "RAM allocation.\n"
            )
            sys.stderr.write(msg)
            sys.stderr.flush()
            os._exit(137)

    t = threading.Thread(target=_watch, daemon=True, name="ks-pool-watchdog")
    t.start()
    return stop_event


@contextlib.contextmanager
def _quiet_worker_pipe_errors() -> Iterator[None]:
    """Suppress benign ``BrokenPipeError`` / ``EOFError`` from pool threads.

    When a smoothing worker is killed (most often by the OOM-killer),
    ``multiprocessing.Pool``'s internal ``_handle_tasks`` /
    ``_handle_results`` daemon threads raise ``BrokenPipeError`` (or
    ``EOFError``) trying to talk to the dead worker's pipe, in the
    window before the watchdog (:func:`_spawn_pool_watchdog`) fires.
    Python's default ``threading.excepthook`` dumps a full traceback per
    such thread -- a wall of noise on top of the watchdog's single
    actionable FATAL message. We install a scoped hook that drops *only*
    those two exception types and delegates everything else to the
    previous hook, restoring it on exit. Scoped to the smoothing pass,
    where the pool threads are the only plausible source of these.
    """
    previous = threading.excepthook

    def _hook(args: threading.ExceptHookArgs) -> None:
        if args.exc_type is not None and issubclass(args.exc_type, (BrokenPipeError, EOFError)):
            return
        previous(args)

    threading.excepthook = _hook
    try:
        yield
    finally:
        threading.excepthook = previous


def _smooth_all_feature_sets(
    *,
    combined_bed: Path,
    feature_sets: list[str],
    features: Features,
    indices: dict[str, HierarchyIndex],
    presmoothed_paths: dict[str, Path],
    smoothed_paths: dict[str, Path],
    threads: int,
    chunk_size: int = 50000,
    preserve_input_order: bool = True,
    is_reads_input: bool = False,
) -> None:
    """Stream ``combined_bed`` through ONE pool that handles every FS.

    Replaces the legacy per-feature-set pool loop. Each worker is
    initialised with the full ``{fs: HierarchyIndex}`` /
    ``{fs: FeaturesForWorker}`` state up front, parses each chunk
    once, and runs the smoothing pipeline for every requested
    feature set on each sequence in the chunk.

    Wins over the per-FS loop:

    * One pool spawn instead of one per feature set (saves ~30 s x
      ``N_fs-1`` of process startup + initargs unpickle on whole-
      genome inputs).
    * One pass over ``combined_bed`` instead of one per feature set
      (saves ``N_fs-1`` x I/O of a multi-GB combined BED).
    * No inter-feature-set idle gaps (previously ~2 min per
      transition on HG002 while the old pool drained, temp files
      were concatenated, and the next pool initialised).

    ``presmoothed_paths`` and ``smoothed_paths`` are each
    ``{feature_set: path}``; either may be empty (when
    ``--no-keep-presmoothed`` / ``--no-smooth``) but not both. The
    caller is responsible for that check.

    Three codepaths, selected by ``(preserve_input_order, is_reads_input)``:

    * **(True, False)** -- *assembly + preserve*. Workers write each
      ``(fs, seq)`` pair's output directly to a per-feature-set temp
      directory; the main process concatenates the per-sequence
      files in input-FASTA order at the end. Output is
      byte-equivalent to the legacy per-FS code. Per-sequence temp
      files are essential here because a single chunk (e.g. chr1)
      can produce hundreds of MB of output BED lines that would
      overload the IPC pipe.
    * **(True, True)** -- *reads + preserve*. Streaming dispatch via
      :meth:`Pool.imap` (ordered): workers return BED lines via IPC,
      main writes in input order. No per-sequence temp files, because
      per-read temp files would scale catastrophically (millions of
      tiny files). Safe because read chunks are uniform-size, so the
      ordered iterator doesn't stall waiting for a single slow chunk
      and per-chunk IPC payloads stay small.
    * **(False, *)** -- *no preserve*. Streaming dispatch via
      :meth:`Pool.imap_unordered`: sequences appear in
      worker-completion order. Fastest path when order doesn't
      matter downstream. Used for reads when output order is
      irrelevant, or for long-read FASTA where the user opted out of
      the temp-file machinery.
    """
    if not presmoothed_paths and not smoothed_paths:
        raise KaryoscopeError(
            "internal error: _smooth_all_feature_sets called with no output paths"
        )

    pool_size = _cpus.resolve_threads(threads)
    features_by_fs = {fs: make_features_for_worker(features, fs) for fs in feature_sets}

    ctx = mp.get_context("spawn")
    # Inherit the main process's root log level so worker INFO lines
    # (per-sequence smoothing progress) appear when the user passed
    # -v, but stay silent at the default WARNING level.
    worker_log_level = logging.getLogger().level

    if preserve_input_order and not is_reads_input:
        # Assembly + preserve: per-sequence temp files (chunks can be
        # multi-GB; can't return via IPC without blowing memory).
        _smooth_with_per_sequence_tempfiles(
            combined_bed=combined_bed,
            feature_sets=feature_sets,
            indices=indices,
            features_by_fs=features_by_fs,
            presmoothed_paths=presmoothed_paths,
            smoothed_paths=smoothed_paths,
            ctx=ctx,
            pool_size=pool_size,
            worker_log_level=worker_log_level,
            chunk_size=chunk_size,
        )
    else:
        # Reads (preserve or not) or assembly + no-preserve:
        # streaming. The ordered flag picks pool.imap vs
        # pool.imap_unordered.
        _smooth_streaming(
            combined_bed=combined_bed,
            feature_sets=feature_sets,
            indices=indices,
            features_by_fs=features_by_fs,
            presmoothed_paths=presmoothed_paths,
            smoothed_paths=smoothed_paths,
            ctx=ctx,
            pool_size=pool_size,
            worker_log_level=worker_log_level,
            chunk_size=chunk_size,
            ordered=preserve_input_order,
        )


def _smooth_streaming(
    *,
    combined_bed: Path,
    feature_sets: list[str],
    indices: dict[str, HierarchyIndex],
    features_by_fs: dict,
    presmoothed_paths: dict[str, Path],
    smoothed_paths: dict[str, Path],
    ctx,
    pool_size: int,
    worker_log_level: int,
    chunk_size: int,
    ordered: bool,
) -> None:
    """Streaming smoothing: workers return BED lines via IPC, main writes.

    No per-(fs, seq) temp files -- workers return their BED lines as
    dicts via the IPC pipe, and the main process writes them straight
    to the per-feature-set output BEDs. Appropriate for inputs where
    each chunk's returned data is small (read-level FASTQ/BAM with
    ~50k reads per chunk = a few MB per chunk per FS), as opposed to
    assembly inputs where a chunk can be a whole chromosome's worth
    of intervals (hundreds of MB per FS, which would overload IPC and
    motivated the per-sequence temp-files codepath).

    Two ordering modes, selected by ``ordered``:

    * **``ordered=True``** -- uses :meth:`Pool.imap`. Chunk results
      are yielded in *input* order regardless of worker-completion
      order; the output BED preserves input-FASTQ/BAM order. Safe for
      reads because all chunks are roughly the same size, so the
      ordered iterator doesn't stall waiting for a single slow chunk.
      Used when ``preserve_input_order=True`` on a FASTQ/BAM input.

    * **``ordered=False``** -- uses :meth:`Pool.imap_unordered`.
      Sequences appear in worker-completion order, which is faster
      because workers don't wait on the head-of-queue chunk but loses
      input order. Used when ``preserve_input_order=False``.
    """
    pre_handles = {fs: presmoothed_paths[fs].open("w") for fs in presmoothed_paths}
    smo_handles = {fs: smoothed_paths[fs].open("w") for fs in smoothed_paths}
    try:
        with ctx.Pool(
            processes=pool_size,
            initializer=worker_initializer,
            initargs=(indices, features_by_fs, feature_sets, None, worker_log_level),
        ) as pool:
            stop_event = _spawn_pool_watchdog(pool)
            try:
                # pool.imap preserves input order; pool.imap_unordered
                # yields whichever result is ready first. Both consume
                # the same chunk iterator.
                map_method = pool.imap if ordered else pool.imap_unordered
                for chunk_result in map_method(
                    process_seq_chunk,
                    chunked_seq_reader(combined_bed, chunk_size),
                ):
                    # chunk_result: {fs: {seq: (smoothed_lines, presmoothed_lines)}}
                    for fs, per_seq in chunk_result.items():
                        pre_h = pre_handles.get(fs)
                        smo_h = smo_handles.get(fs)
                        for _seq, (smo_lines, pre_lines) in per_seq.items():
                            if pre_h is not None:
                                pre_h.writelines(pre_lines)
                            if smo_h is not None:
                                smo_h.writelines(smo_lines)
            finally:
                stop_event.set()
    finally:
        for h in pre_handles.values():
            h.close()
        for h in smo_handles.values():
            h.close()


def _smooth_with_per_sequence_tempfiles(
    *,
    combined_bed: Path,
    feature_sets: list[str],
    indices: dict[str, HierarchyIndex],
    features_by_fs: dict,
    presmoothed_paths: dict[str, Path],
    smoothed_paths: dict[str, Path],
    ctx,
    pool_size: int,
    worker_log_level: int,
    chunk_size: int,
) -> None:
    """Assembly-mode smoothing: workers write per-(fs, seq) temp files.

    Workers write each ``(feature_set, sequence)`` pair's smoothed /
    presmoothed BED lines directly to
    ``{tmpdir}/{fs}/{safe_seq_name}.smo`` / ``.pre`` and return only
    line-count metadata, avoiding the multi-GB IPC transfer that
    whole-chromosome chunks would otherwise produce. At the end the
    main process concatenates the per-sequence files in input-FASTA
    order to produce each feature set's final BED.

    Input sequence order is captured at dispatch time: we wrap
    :func:`chunked_seq_reader` with a peek-and-yield generator that
    records each chunk's sequence names before yielding the chunk
    to the pool. Collisions between distinct seq names that sanitise
    to the same filename component are detected here (real genomes
    don't trigger this).
    """
    import tempfile

    # Temp dir lives next to the output so cleanup happens on the
    # same filesystem (no cross-device shutil.copy cost) and so the
    # temp dir is on /scratch when the output is.
    anchor = (
        next(iter(smoothed_paths.values()))
        if smoothed_paths
        else next(iter(presmoothed_paths.values()))
    )
    out_dir = anchor.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    input_seq_order: list[str] = []
    seen_seqs: set[str] = set()
    safe_to_orig: dict[str, str] = {}

    def _record_and_yield():
        """Peek each chunk for sequence names, record in input order, yield."""
        for chunk in chunked_seq_reader(combined_bed, chunk_size):
            for raw in chunk:
                seq = raw.partition("\t")[0]
                if seq and seq not in seen_seqs:
                    seen_seqs.add(seq)
                    input_seq_order.append(seq)
                    safe = safe_filename(seq)
                    existing = safe_to_orig.get(safe)
                    if existing is not None and existing != seq:
                        raise KaryoscopeError(
                            f"sequence names {existing!r} and {seq!r} both "
                            f"sanitise to the same temp-file path component "
                            f"{safe!r}. Rename one of the inputs."
                        )
                    safe_to_orig[safe] = seq
            yield chunk

    with tempfile.TemporaryDirectory(prefix="ks_smooth_", dir=out_dir) as tmpdir_str:
        tmpdir = Path(tmpdir_str)
        tmpdir_by_fs = {fs: tmpdir / fs for fs in feature_sets}
        for fs_dir in tmpdir_by_fs.values():
            fs_dir.mkdir(parents=True, exist_ok=True)

        with ctx.Pool(
            processes=pool_size,
            initializer=worker_initializer,
            initargs=(
                indices,
                features_by_fs,
                feature_sets,
                tmpdir_by_fs,
                worker_log_level,
            ),
        ) as pool:
            stop_event = _spawn_pool_watchdog(pool)
            try:
                # Drain the iterator -- workers have already written
                # the per-(fs, seq) files; we just need to consume the
                # metadata returns to drive imap_unordered to completion.
                for _chunk_meta in pool.imap_unordered(
                    process_seq_chunk,
                    _record_and_yield(),
                ):
                    pass
            finally:
                stop_event.set()

        # Concatenate per-(fs, seq) temp files in input sequence order
        # to produce each feature set's final BEDs. Each (fs, kind)
        # target is independent (writes a distinct output file from a
        # distinct source temp dir), so we parallelise across the 6 x 2
        # = 12 typical-case targets with a thread pool. Pure I/O work,
        # so threads are fine -- the GIL is released across the
        # blocking read/write inside shutil.copyfileobj, and we get
        # genuine parallel disk throughput on /scratch SSD.
        concat_tasks: list[tuple[Path, Path, str]] = []
        for fs in feature_sets:
            fs_dir = tmpdir_by_fs[fs]
            if fs in smoothed_paths:
                concat_tasks.append((smoothed_paths[fs], fs_dir, ".smo"))
            if fs in presmoothed_paths:
                concat_tasks.append((presmoothed_paths[fs], fs_dir, ".pre"))
        concat_threads = min(pool_size, len(concat_tasks)) if concat_tasks else 1
        logger.info(
            "concat pass: %d temp-file group(s) (threads=%d)",
            len(concat_tasks),
            concat_threads,
        )
        t_concat_start = time.perf_counter()
        with ThreadPoolExecutor(max_workers=concat_threads) as exe:
            # list() forces .result() on each future, so any worker
            # exception is re-raised here rather than silently swallowed.
            list(
                exe.map(
                    lambda task: _concat_per_sequence_to_bed(
                        task[0], task[1], task[2], input_seq_order
                    ),
                    concat_tasks,
                )
            )
        logger.info("concat pass complete in %.1fs", time.perf_counter() - t_concat_start)
        # tmpdir cleaned up automatically on context manager exit


def _concat_per_sequence_to_bed(
    dest: Path,
    fs_dir: Path,
    ext: str,
    input_seq_order: list[str],
) -> None:
    """Concatenate per-(fs, seq) temp files into one BED in input order.

    Streams ``{fs_dir}/{safe_filename(seq)}{ext}`` for each ``seq`` in
    ``input_seq_order`` via :func:`shutil.copyfileobj` (64 KB block
    streaming -- never holds a whole file in memory). Missing temp
    files are skipped silently: a sequence may have produced no
    output for this feature set if every interval merged away.

    Pulled out so the bgzip-pass / concat-pass parallelisation in
    :func:`_smooth_with_per_sequence_tempfiles` can call it from a
    thread pool -- each (feature_set, kind) target reads a distinct
    temp dir and writes a distinct output file, so there is no
    cross-task contention to worry about.
    """
    with dest.open("wb") as out_h:
        for seq in input_seq_order:
            p = fs_dir / f"{safe_filename(seq)}{ext}"
            if not p.exists():
                continue
            with p.open("rb") as src:
                shutil.copyfileobj(src, out_h)


def _human_bytes(n: int) -> str:
    """Render a byte count as ``"X.YY GB"`` (>=1 GB) or ``"X.Y MB"``."""
    if n >= 1024**3:
        return f"{n / 1024**3:.2f} GB"
    return f"{n / 1024**2:.1f} MB"


def _peak_child_rss_bytes() -> int | None:
    """Peak resident set size across every child process reaped so far.

    This is the number that decides how much machine an ``annotate`` run
    needs. KaryoScope's own footprint is small next to the ``hks``
    invocations it waits on, and those are not visible in this process's
    own ``ru_maxrss``.

    It is a high-water mark over all reaped children, so it only ever rises
    — report it as a running peak, never as one step's cost. Returns None
    where ``resource`` is unavailable.
    """
    try:
        import resource
    except ImportError:  # pragma: no cover — not a platform we ship for
        return None
    raw = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    # ru_maxrss is kilobytes on Linux and bytes on macOS.
    return raw if sys.platform == "darwin" else raw * 1024


def _rate(nbytes: int, seconds: float) -> str:
    """Render a throughput, or ``"-"`` when the elapsed time is unusable."""
    if seconds <= 0:
        return "-"
    return f"{nbytes / 1024**3 / seconds:.2f} GB/s"


# --- HKS backend ------------------------------------------------------


def _run_hks_backend(
    *,
    manifest,
    db_dir: Path,
    input_path: Path,
    prefix: str,
    output_dir: Path,
    requested: list[str],
    smooth: bool,
    keep_presmoothed: bool,
    presmoothed_paths: dict[str, Path],
    smoothed_paths: dict[str, Path],
    threads: int,
    k: int,
    progress: Progress = SILENT,
) -> None:
    """Run the HKS lookup and optional smoothing for every requested feature set.

    Unlike the KMC backend (one combined query, then translate integer feature
    ids), HKS queries one ``.hksf`` per feature set and reads label names
    directly. Each feature set is processed independently:

    1. ``hks lookup`` against ``<basename>.<fs>.hksf`` -> the presmoothed BED.
    2. If ``smooth``: ``hks smooth`` (using ``<basename>.<fs>.hierarchy.txt``)
       reads that same file -> the smoothed BED.

    There is no conversion step between them. ``hks`` is told the output shape
    KaryoScope wants -- headerless, ``novel`` for misses -- so the lookup output
    *is* the presmoothed BED, and smooth reads it in place. When the caller does
    not want the presmoothed BED kept, the lookup still has to write somewhere
    for smooth to read, so it goes to a per-feature-set temp file deleted after
    each set.

    ``k`` is the query k-mer length (``manifest.kmer.size`` unless overridden
    for a variable-k index).
    """
    base_path = db_dir / (manifest.index.basename + ".hksb")

    # Reads emit integer query ranks instead of names: HKS otherwise loads
    # every read name into memory (~10 GB at hundreds of millions of reads),
    # and read names carry no downstream meaning (unlike assembly contig names,
    # which map to karyotype chromosomes).
    is_reads = _is_reads_input(input_path)

    # A BAM is converted to FASTA once, up front: run_hks_lookup would
    # otherwise re-run samtools fasta for every feature set.
    query_path = input_path
    tmp_fasta: Path | None = None
    if input_path.suffix.lower() == ".bam":
        logger.info(
            "converting BAM %s to FASTA once for %d feature set(s)",
            input_path.name,
            len(requested),
        )
        tmp_fasta = convert_bam_to_fasta(input_path, output_dir, capture=True)
        query_path = tmp_fasta

    t_hks_start = time.perf_counter()
    tracker = progress.track(requested)
    try:
        for fs in requested:
            t_fs = time.perf_counter()
            fs_file = db_dir / f"{manifest.index.basename}.{fs}.hksf"
            hierarchy_file = db_dir / f"{manifest.index.basename}.{fs}.hierarchy.txt"
            # The lookup output is already the presmoothed BED, so when it is being
            # kept it is written straight to its final home rather than copied
            # there. Otherwise smooth still needs it on disk to read, and it is a
            # temp file we drop afterwards.
            lookup_out = (
                presmoothed_paths[fs]
                if keep_presmoothed
                else output_dir / f"{prefix}.{fs}.lookup_raw.tmp.bed"
            )

            logger.info(
                "running hks lookup for feature set %r on %s (threads=%d)",
                fs,
                input_path.name,
                threads,
            )
            t_lookup = time.perf_counter()
            run_hks_lookup(
                base_path=base_path,
                feature_set_file=fs_file,
                k=k,
                input_path=query_path,
                output_path=lookup_out,
                threads=threads,
                report_query_names=not is_reads,
                capture=True,
            )
            dt_lookup = time.perf_counter() - t_lookup
            if not lookup_out.is_file():
                raise KaryoscopeError(f"hks lookup did not produce expected output at {lookup_out}")

            # Timed and sized per phase rather than per feature set. The two do
            # very different work -- the lookup is the parallel k-mer query, the
            # smooth a largely serial pass over what it wrote -- so a single
            # per-feature-set number cannot say which of them a change moved.
            lookup_bytes = lookup_out.stat().st_size
            # Deliberately no throughput here. A lookup's time is dominated by
            # loading the index and querying the input, not by writing its
            # output, so output-bytes-per-second would be a rate of nothing --
            # it read 0.02 GB/s on a real run purely because the index load is
            # large and the BED is small. `hks -vv` reports the phases that do
            # have meaningful rates.
            logger.info(
                "hks lookup for %r wrote %s in %.1fs",
                fs,
                _human_bytes(lookup_bytes),
                dt_lookup,
            )

            try:
                if smooth:
                    t_smo = time.perf_counter()
                    logger.info("running hks smooth for feature set %r", fs)
                    run_hks_smooth(
                        hierarchy_file=hierarchy_file,
                        input_path=lookup_out,
                        output_path=smoothed_paths[fs],
                        threads=threads,
                        capture=True,
                    )
                    dt_smo = time.perf_counter() - t_smo
                    smoothed_bytes = smoothed_paths[fs].stat().st_size
                    logger.info(
                        "hks smooth for %r wrote %s in %.1fs (read %s at %s)",
                        fs,
                        _human_bytes(smoothed_bytes),
                        dt_smo,
                        _human_bytes(lookup_bytes),
                        _rate(lookup_bytes, dt_smo),
                    )
            finally:
                if not keep_presmoothed:
                    try:
                        lookup_out.unlink()
                    except OSError as exc:
                        logger.warning(
                            "could not remove temp lookup output %s: %s", lookup_out, exc
                        )
            peak = _peak_child_rss_bytes()
            if peak is not None:
                logger.info("peak hks memory so far: %s", _human_bytes(peak))
            # Reported here rather than per sub-step: one line per feature set
            # is the granularity a user waiting on the run actually needs, and
            # HKS processes them strictly in sequence so the counter is honest.
            tracker.step(fs, time.perf_counter() - t_fs)
    finally:
        if tmp_fasta is not None:
            tmp_fasta.unlink(missing_ok=True)

    peak = _peak_child_rss_bytes()
    logger.info(
        "hks backend complete in %.1fs (%d feature set(s)%s)",
        time.perf_counter() - t_hks_start,
        len(requested),
        f", peak hks memory {_human_bytes(peak)}" if peak is not None else "",
    )


# --- main entry point -------------------------------------------------


def annotate(
    *,
    input_path: Path,
    output_dir: Path,
    db_root: Path,
    db_id: str | None = None,
    feature_sets: list[str] | None = None,
    threads: int = 0,
    smooth: bool = True,
    keep_presmoothed: bool = True,
    keep_intermediates: bool = False,
    bgzip: bool = True,
    preserve_input_order: bool = True,
    force: bool = False,
    k: int | None = None,
    check_space: bool = True,
    progress: Progress = SILENT,
) -> AnnotateResult:
    """Run the full annotate pipeline for one input FASTA.

    Parameters
    ----------
    input_path
        FASTA (plain or gzipped) to annotate.
    output_dir
        Directory to write output BEDs to. Created if absent.
    db_root
        KaryoScope database root (typically ``ensure_db_root(...)``).
    db_id
        Specific database id to use, or ``None`` to pick the unique
        installed database (errors if there are zero or many).
    feature_sets
        Restrict output to these feature sets. ``None`` means "all sets
        declared in the database's manifest".
    threads
        Threads for both the C++ k-mer query and the smoothing pool.
        ``0`` means auto (the CPUs this process may actually use --
        a SLURM allocation or CPU affinity, not the machine's core count).
    smooth
        Produce the hierarchy-smoothed BED. Default: ``True``.
    keep_presmoothed
        Keep the presmoothed BED. Default: ``True``. When ``False`` and
        ``smooth=True``, only the smoothed output is written.
    keep_intermediates
        Keep the combined ``.combined.presmoothed.featureIDs.bed`` from
        the C++ step. Default: delete after processing.
    bgzip
        bgzip the per-feature-set output BEDs. Default: yes.
    force
        Regenerate the combined intermediate even when a complete one
        already exists on disk. Default: ``False`` -- a rerun reuses a
        verified combined BED from a previous (possibly crashed) run and
        skips the expensive ``get_featureIDs`` step. "Verified" means the
        BED is present and its ``.done`` completion marker matches its
        size/mtime; a partial file left by a killed run has no matching
        marker and is regenerated regardless of this flag.
    check_space
        Estimate the output footprint and refuse to start if ``output_dir``
        can't hold it. Default: ``True``. The estimate is derived from the
        input size, so pass ``False`` if it misjudges an unusual input.
    progress
        Milestone reporter for stdout. Defaults to silence, so importing
        KaryoScope as a library never prints; the CLI passes an enabled one.

    Raises
    ------
    InsufficientDiskSpaceError
        If ``check_space`` is True and ``output_dir`` looks too small, or if
        the filesystem fills up during the run regardless.
    MissingDependencyError
        If an external tool this run needs isn't installed.
    KaryoscopeError
        If ``smooth=False`` and ``keep_presmoothed=False`` (no output
        would be produced), if the input file doesn't exist, if the
        requested feature set isn't in the manifest, or if the
        hierarchy fails validation (only when smoothing is enabled).
    DatabaseNotFoundError
        If the database can't be resolved or its layout is broken.
    ToolNotFoundError
        If the k-mer query binary (``get_featureIDs`` for KMC, ``hks`` for HKS)
        or ``bgzip`` (when requested) isn't found.
    ExternalToolError
        If a subprocess exits non-zero.
    """
    if not smooth and not keep_presmoothed:
        raise KaryoscopeError(
            "no output would be produced: cannot combine --no-smooth with "
            "--no-keep-presmoothed. Choose at least one output type."
        )
    if not input_path.is_file():
        raise KaryoscopeError(f"input file not found: {input_path}")

    db_id_resolved, db_dir = resolve_database(db_root, db_id)
    manifest = validate_database_layout(db_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolve the query k-mer length. Defaults to the manifest's size; an
    # explicit override is only honoured on a variable-k HKS index (kmer.type
    # == "variable"), which can answer any k <= max_size. A fixed-k index (KMC,
    # or an HKS index built without --variable-k) can only be queried at its s.
    query_k = _resolve_query_k(manifest, k, db_id_resolved)

    # Pick feature sets, validating the user's choices against the manifest.
    available = list(manifest.feature_sets)
    if feature_sets is None:
        requested = available
    else:
        unknown = [fs for fs in feature_sets if fs not in available]
        if unknown:
            raise KaryoscopeError(
                f"feature set(s) {unknown!r} not declared in {db_id_resolved}'s "
                f"manifest. Available: {available!r}"
            )
        requested = list(feature_sets)
    logger.info("annotating %s against %s, sets=%s", input_path, db_id_resolved, requested)

    # Preflight. Both checks are cheap and both catch failures that would
    # otherwise land many minutes in, after the k-mer query has run: a
    # missing bgzip only bites at the very last step, and a full disk only
    # once gigabytes have been written.
    preflight.require(
        _annotate_dependencies(index_type=manifest.index.type, input_path=input_path, bgzip=bgzip),
        context=f"annotate against {db_id_resolved}",
    )
    _cpus.warn_if_oversubscribed(threads, what=f"annotate against {db_id_resolved}")
    input_bases = estimate_input_bases(input_path)
    needed_bytes = estimate_output_bytes(
        input_bases=input_bases,
        n_feature_sets=len(requested),
        keep_presmoothed=keep_presmoothed,
        smooth=smooth,
        index_type=manifest.index.type,
    )
    logger.info(
        "estimated output footprint: %s for %d feature set(s) over ~%.2f Gbp of input",
        _human_bytes(needed_bytes),
        len(requested),
        input_bases / 1e9,
    )
    diskspace.require_free_space(
        output_dir,
        needed_bytes,
        what=f"annotating {input_path.name} ({len(requested)} feature set(s))",
        estimated=True,
        hint=(
            "Note that --bgzip shrinks the final output but not the peak: every "
            "BED is written in full before the compression pass starts.\n"
            "Options:\n"
            "  - write to a larger filesystem with --outdir\n"
            "  - annotate fewer feature sets at a time with --feature-set\n"
            "  - drop one of the two outputs with --no-keep-presmoothed or --no-smooth\n"
            "  - pass --no-space-check if this estimate looks wrong for your input"
        ),
        skip=not check_space,
    )

    # Announce the run before the first expensive step. Everything below
    # this point can take twenty minutes, and until now the terminal stayed
    # blank for all of it.
    n_threads = _cpus.resolve_threads(threads)
    progress.start(
        f"Annotating {input_path.name} against {db_id_resolved}",
        f"{len(requested)} feature set(s), {n_threads} thread(s), "
        f"~{_human_bytes(needed_bytes)} estimated output",
    )

    t_annotate_start = time.perf_counter()

    # Parse features.tsv up front so we fail fast if it's malformed. It maps
    # integer feature ids to names for the KMC backend; the HKS backend reads
    # names from the index and omits features.tsv entirely (manifest.features is
    # then None), so parse it only when present.
    features: Features | None = None
    if manifest.features is not None:
        features = parse_features(db_dir / manifest.features)
        missing_from_table = [fs for fs in requested if fs not in features.feature_sets]
        if missing_from_table:
            raise KaryoscopeError(
                f"feature set(s) {missing_from_table!r} declared in manifest but "
                f"missing from features.tsv columns ({features.feature_sets!r})"
            )

    # Parse + validate the hierarchy if smoothing is enabled. The
    # validation is hard-fail here (in contrast to `info`, where it's
    # a warning) because producing smoothed output against a broken
    # hierarchy would silently give wrong results.
    hierarchy: Hierarchy | None = None
    indices: dict[str, HierarchyIndex] = {}
    if smooth:
        hierarchy = parse_hierarchy(db_dir / manifest.hierarchy)
        # Build the features-columns map for the cross-validation check. Only
        # available for the KMC backend (features.tsv); for HKS there is no
        # feature table, so skip that particular check (pass None).
        feature_columns: dict[str, set[str]] | None = (
            {fs: {row[fs] for row in features.table.values()} for fs in requested}
            if features is not None
            else None
        )
        issues = validate_hierarchy(hierarchy, feature_columns=feature_columns)
        if issues:
            raise KaryoscopeError(
                "hierarchy.tsv failed validation; refusing to produce "
                "smoothed output:\n  - " + "\n  - ".join(issues)
            )
        for fs in requested:
            if fs not in hierarchy.feature_sets():
                raise KaryoscopeError(
                    f"feature set {fs!r} has no rows in hierarchy.tsv; "
                    "smoothing this set is not possible. Re-run with "
                    "--no-smooth or restrict --feature-set."
                )
            indices[fs] = HierarchyIndex.from_hierarchy(hierarchy, fs)

    input_basename = _derive_input_basename(input_path)
    prefix = f"{input_basename}.{db_id_resolved}"
    # When an explicit k is used, tag the outputs with it so a k-sweep into one
    # directory doesn't overwrite itself (and the default run stays unchanged).
    if k is not None:
        prefix = f"{prefix}.k{query_k}"

    # Compute output paths (uncompressed names; bgzip later if requested).
    presmoothed_paths: dict[str, Path] = (
        {fs: output_dir / f"{prefix}.{fs}.presmoothed.bed" for fs in requested}
        if keep_presmoothed
        else {}
    )
    smoothed_paths: dict[str, Path] = (
        {fs: output_dir / f"{prefix}.{fs}.smoothed.bed" for fs in requested} if smooth else {}
    )

    # --- Backend dispatch -------------------------------------------------
    # The combined-BED intermediate is a KMC-only artifact; the HKS backend
    # writes per-feature-set BEDs directly and has nothing to keep here.
    combined_kept: Path | None = None

    if manifest.index.type == "hks":
        _run_hks_backend(
            manifest=manifest,
            db_dir=db_dir,
            input_path=input_path,
            prefix=prefix,
            output_dir=output_dir,
            requested=requested,
            smooth=smooth,
            keep_presmoothed=keep_presmoothed,
            presmoothed_paths=presmoothed_paths,
            smoothed_paths=smoothed_paths,
            threads=threads,
            k=query_k,
            progress=progress,
        )
    else:  # "kmc" -- the only other supported type (guaranteed by parse_manifest)
        # Run the C++ helper -- unless a complete combined BED from a prior
        # run is already on disk. get_featureIDs is the most expensive and
        # most memory-hungry step; reusing a verified result lets a user who
        # was OOM-killed during smoothing simply rerun (e.g. with fewer
        # --threads) and resume straight into the smoothing pass, instead of
        # paying for -- and risking another OOM in -- the k-mer query again.
        # The combined BED is feature-set- and thread-count-agnostic, so
        # reuse stays correct even if those args changed between runs.
        kmc_db_basename = db_dir / manifest.index.basename
        combined_bed = combined_bed_path(output_dir, prefix)

        if not force and combined_bed_is_complete(combined_bed):
            logger.info(
                "reusing existing combined BED from a previous run (%s): %s "
                "-- skipping get_featureIDs. Pass --force to regenerate.",
                _human_bytes(combined_bed.stat().st_size),
                combined_bed,
            )
            progress.note(
                f"reusing the combined BED from a previous run "
                f"({_human_bytes(combined_bed.stat().st_size)}); skipping the k-mer query"
            )
        else:
            logger.info(
                "running get_featureIDs on %s (threads=%d); this may take several minutes",
                input_path.name,
                threads,
            )
            t_kmc_start = time.perf_counter()
            combined_bed = run_get_featureids(
                db_path=kmc_db_basename,
                input_path=input_path,
                output_dir=output_dir,
                threads=threads,
                prefix=prefix,
                capture=True,
            )
            if not combined_bed.is_file():
                raise KaryoscopeError(
                    f"get_featureIDs did not produce expected output at {combined_bed}"
                )
            logger.info(
                "ran get_featureIDs in %.1fs (combined BED: %s)",
                time.perf_counter() - t_kmc_start,
                _human_bytes(combined_bed.stat().st_size),
            )
            # Named stages rather than [i/N]: the KMC backend runs one
            # combined query and then one streaming smoothing pass over
            # every feature set at once, so there is no per-feature-set
            # completion moment to count.
            progress.stage("k-mer query", time.perf_counter() - t_kmc_start)
        logger.debug("combined BED at %s", combined_bed)

        # Run the smoothing pass. One pool initialised with every
        # requested feature set's state; each chunk is processed for all
        # feature sets in one worker invocation. See
        # :func:`_smooth_all_feature_sets` for the architectural rationale.
        # When smoothing is off we use the simpler in-process splitter
        # (no need to fork workers for a one-line translation).
        if smooth:
            is_reads = _is_reads_input(input_path)
            logger.info(
                "smoothing pass: %d feature set(s), threads=%d",
                len(requested),
                threads,
            )
            logger.debug(
                "smoothing pass with threads=%d, preserve_input_order=%s, "
                "is_reads_input=%s, feature_sets=%s",
                threads,
                preserve_input_order,
                is_reads,
                requested,
            )
            t_smooth_start = time.perf_counter()
            # Quiet the benign BrokenPipe/EOF tracebacks the pool's daemon
            # threads emit if a worker is OOM-killed, so the watchdog's
            # FATAL message isn't buried under a wall of noise.
            with _quiet_worker_pipe_errors():
                _smooth_all_feature_sets(
                    combined_bed=combined_bed,
                    feature_sets=requested,
                    features=features,
                    indices=indices,
                    presmoothed_paths=presmoothed_paths,
                    smoothed_paths=smoothed_paths,
                    threads=threads,
                    preserve_input_order=preserve_input_order,
                    is_reads_input=is_reads,
                )
            logger.info("smoothing pass complete in %.1fs", time.perf_counter() - t_smooth_start)
            progress.stage(
                f"smoothing {len(requested)} feature set(s)",
                time.perf_counter() - t_smooth_start,
            )
        else:
            # Only presmoothed output, no smoothing.
            logger.info("splitting combined BED into %d per-feature-set BED(s)", len(requested))
            t_split_start = time.perf_counter()
            _split_combined_bed(combined_bed, requested, features, presmoothed_paths)
            logger.info("split complete in %.1fs", time.perf_counter() - t_split_start)
            progress.stage(
                f"splitting into {len(requested)} feature set(s)",
                time.perf_counter() - t_split_start,
            )

        # Tidy up the combined intermediate unless asked to keep it. Remove
        # its completion marker alongside it so no dangling marker is left
        # pointing at a deleted file.
        if not keep_intermediates:
            try:
                combined_bed.unlink()
                clear_combined_marker(combined_bed)
                logger.debug("removed combined intermediate %s", combined_bed)
            except OSError as e:
                logger.warning("could not remove intermediate %s: %s", combined_bed, e)
                combined_kept = combined_bed
        else:
            combined_kept = combined_bed

    # bgzip (or not).
    if bgzip:
        n_to_bgzip = sum(1 for fs in requested if fs in presmoothed_paths) + sum(
            1 for fs in requested if fs in smoothed_paths
        )
        logger.info("bgzip pass: %d BED(s) (threads=%d each)", n_to_bgzip, threads)
        t_bgzip_start = time.perf_counter()
        for fs in requested:
            if fs in presmoothed_paths:
                presmoothed_paths[fs] = bgzip_file(presmoothed_paths[fs], threads=threads)
            if fs in smoothed_paths:
                smoothed_paths[fs] = bgzip_file(smoothed_paths[fs], threads=threads)
        logger.info("bgzip pass complete in %.1fs", time.perf_counter() - t_bgzip_start)
        # Worth its own line: compressing 12 BEDs of a human diploid run
        # takes minutes, and it happens after the last feature-set line, so
        # without this the run looks finished-but-hung right at the end.
        progress.stage(f"bgzip ({n_to_bgzip} file(s))", time.perf_counter() - t_bgzip_start)

    n_outputs = len(presmoothed_paths) + len(smoothed_paths)
    logger.info(
        "annotate complete in %.1fs (%d output BED(s))",
        time.perf_counter() - t_annotate_start,
        n_outputs,
    )

    return AnnotateResult(
        presmoothed_paths=presmoothed_paths,
        smoothed_paths=smoothed_paths,
        combined_intermediate=combined_kept,
    )


# --- multi-input batch (HKS backend) ----------------------------------


def _run_hks_backend_batch(
    *,
    manifest,
    db_dir: Path,
    input_paths: list[Path],
    prefixes: dict[Path, str],
    output_dir: Path,
    requested: list[str],
    smooth: bool,
    keep_presmoothed: bool,
    presmoothed_by_input: dict[Path, dict[str, Path]],
    smoothed_by_input: dict[Path, dict[str, Path]],
    threads: int,
    k: int,
    progress: Progress = SILENT,
) -> None:
    """HKS lookup+smoothing for MANY inputs, one index load per feature set.

    The key difference from :func:`_run_hks_backend` (single input) is loop
    order: the feature set is the OUTER loop, and for each set every input is
    queried in a single ``hks lookup`` invocation (see
    :func:`karyoscope.core.io.hks.run_hks_lookup_batch`). The (~6 GB base + ~3 GB
    feature-set) index is therefore loaded once per feature set for the whole
    cohort, instead of once per (input, feature set) pair.

    ``--report-query-names`` is a single per-invocation ``hks`` flag, so inputs
    are grouped by :func:`_is_reads_input` (reads emit query ranks, assemblies
    emit names) and each group gets its own batched call — at most two per
    feature set, still far fewer than one-per-input. After the queries, each
    input's raw TSV is converted (presmoothed) and/or smoothed independently,
    exactly as in the single-input backend.
    """
    base_path = db_dir / (manifest.index.basename + ".hksb")

    t_hks_start = time.perf_counter()
    tracker = progress.track(requested)
    for fs in requested:
        t_fs = time.perf_counter()
        fs_file = db_dir / f"{manifest.index.basename}.{fs}.hksf"
        hierarchy_file = db_dir / f"{manifest.index.basename}.{fs}.hierarchy.txt"
        # As in the single-input path: the lookup output already IS the
        # presmoothed BED, so when it is being kept it goes straight to its
        # final home and no copy is ever made.
        lookup_by_input = {
            p: (
                presmoothed_by_input[p][fs]
                if keep_presmoothed
                else output_dir / f"{prefixes[p]}.{fs}.lookup_raw.tmp.bed"
            )
            for p in input_paths
        }

        # One batched lookup per report-query-names group (reads vs assemblies).
        t_lookup = time.perf_counter()
        for is_reads_group in (False, True):
            group = [p for p in input_paths if _is_reads_input(p) is is_reads_group]
            if not group:
                continue
            logger.info(
                "hks lookup: feature set %r over %d input(s) (reads=%s, threads=%d)",
                fs,
                len(group),
                is_reads_group,
                threads,
            )
            run_hks_lookup_batch(
                base_path=base_path,
                feature_set_file=fs_file,
                k=k,
                io_pairs=[(p, lookup_by_input[p]) for p in group],
                threads=threads,
                report_query_names=not is_reads_group,
                capture=True,
            )

        logger.info(
            "hks lookup for feature set %r over %d input(s) took %.1fs",
            fs,
            len(input_paths),
            time.perf_counter() - t_lookup,
        )

        # Per-input smoothing. Timed separately from the lookup so a
        # batch-vs-per-input comparison can attribute any difference to the
        # query or to this tail, rather than only seeing one wall-clock
        # number for the whole feature set.
        t_smooth = 0.0
        for p in input_paths:
            lookup_out = lookup_by_input[p]
            if not lookup_out.is_file():
                raise KaryoscopeError(f"hks lookup did not produce expected output at {lookup_out}")
            try:
                if smooth:
                    t0 = time.perf_counter()
                    run_hks_smooth(
                        hierarchy_file=hierarchy_file,
                        input_path=lookup_out,
                        output_path=smoothed_by_input[p][fs],
                        threads=threads,
                        capture=True,
                    )
                    t_smooth += time.perf_counter() - t0
            finally:
                if not keep_presmoothed:
                    try:
                        lookup_out.unlink()
                    except OSError as exc:
                        logger.warning(
                            "could not remove temp lookup output %s: %s", lookup_out, exc
                        )

        logger.info(
            "feature set %r: smooth %.1fs (summed over %d input(s))",
            fs,
            t_smooth,
            len(input_paths),
        )
        peak = _peak_child_rss_bytes()
        if peak is not None:
            logger.info("peak hks memory so far: %s", _human_bytes(peak))
        tracker.step(fs, time.perf_counter() - t_fs)

    logger.info(
        "hks batch backend complete in %.1fs (%d feature set(s) x %d input(s))",
        time.perf_counter() - t_hks_start,
        len(requested),
        len(input_paths),
    )


def annotate_batch(
    *,
    input_paths: list[Path],
    output_dir: Path,
    db_root: Path,
    db_id: str | None = None,
    feature_sets: list[str] | None = None,
    threads: int = 0,
    smooth: bool = True,
    keep_presmoothed: bool = True,
    keep_intermediates: bool = False,
    bgzip: bool = True,
    preserve_input_order: bool = True,
    force: bool = False,
    k: int | None = None,
    check_space: bool = True,
    progress: Progress = SILENT,
) -> dict[Path, AnnotateResult]:
    """Annotate several inputs, loading the index once per feature set (HKS).

    For the **HKS** backend, all inputs are queried against each feature set in a
    single ``hks lookup`` (one index load per feature set for the whole cohort,
    versus one per input), then each input's TSVs are converted/smoothed/bgzipped
    independently. All outputs go into the shared ``output_dir``; filenames are
    prefixed by the input basename (as in single-input mode) so they never
    collide.

    For the **KMC** backend there is no batch query primitive, so inputs are
    annotated one at a time via :func:`annotate` (no regression, no speed-up).

    Returns ``{input_path: AnnotateResult}``. ``preserve_input_order`` is accepted
    for signature parity with :func:`annotate` but has no effect on the HKS
    backend (its per-sequence output order is fixed by the query file).
    """
    if not input_paths:
        return {}

    # A single input has nothing to batch — delegate to the exact, battle-tested
    # single-input path (covers both HKS and KMC backends). This makes
    # annotate_batch a safe universal entry point for callers that may pass one
    # or many inputs (CLI, karyotype, scaffold).
    if len(input_paths) == 1:
        p = input_paths[0]
        return {
            p: annotate(
                input_path=p,
                output_dir=output_dir,
                db_root=db_root,
                db_id=db_id,
                feature_sets=feature_sets,
                threads=threads,
                smooth=smooth,
                keep_presmoothed=keep_presmoothed,
                keep_intermediates=keep_intermediates,
                bgzip=bgzip,
                preserve_input_order=preserve_input_order,
                force=force,
                k=k,
                check_space=check_space,
                progress=progress,
            )
        }

    if not smooth and not keep_presmoothed:
        raise KaryoscopeError(
            "no output would be produced: cannot combine --no-smooth with "
            "--no-keep-presmoothed. Choose at least one output type."
        )
    for p in input_paths:
        if not p.is_file():
            raise KaryoscopeError(f"input file not found: {p}")

    db_id_resolved, db_dir = resolve_database(db_root, db_id)
    manifest = validate_database_layout(db_dir)

    # KMC backend: no batch query primitive — annotate inputs one at a time.
    if manifest.index.type != "hks":
        return {
            p: annotate(
                input_path=p,
                output_dir=output_dir,
                db_root=db_root,
                db_id=db_id,
                feature_sets=feature_sets,
                threads=threads,
                smooth=smooth,
                keep_presmoothed=keep_presmoothed,
                keep_intermediates=keep_intermediates,
                bgzip=bgzip,
                preserve_input_order=preserve_input_order,
                force=force,
                k=k,
                check_space=check_space,
                progress=progress,
            )
            for p in input_paths
        }

    output_dir.mkdir(parents=True, exist_ok=True)
    query_k = _resolve_query_k(manifest, k, db_id_resolved)

    # Feature-set selection (identical rules to annotate()).
    available = list(manifest.feature_sets)
    if feature_sets is None:
        requested = available
    else:
        unknown = [fs for fs in feature_sets if fs not in available]
        if unknown:
            raise KaryoscopeError(
                f"feature set(s) {unknown!r} not declared in {db_id_resolved}'s "
                f"manifest. Available: {available!r}"
            )
        requested = list(feature_sets)

    # Validate the hierarchy up front when smoothing (HKS smooths via `hks smooth`
    # against the per-feature-set hierarchy.txt, but we still fail fast on a broken
    # hierarchy.tsv, matching the single-input path).
    if smooth:
        hierarchy = parse_hierarchy(db_dir / manifest.hierarchy)
        issues = validate_hierarchy(hierarchy, feature_columns=None)
        if issues:
            raise KaryoscopeError(
                "hierarchy.tsv failed validation; refusing to produce "
                "smoothed output:\n  - " + "\n  - ".join(issues)
            )
        for fs in requested:
            if fs not in hierarchy.feature_sets():
                raise KaryoscopeError(
                    f"feature set {fs!r} has no rows in hierarchy.tsv; "
                    "smoothing this set is not possible. Re-run with "
                    "--no-smooth or restrict --feature-set."
                )

    logger.info(
        "annotating %d inputs against %s, sets=%s",
        len(input_paths),
        db_id_resolved,
        requested,
    )

    # Preflight, as in the single-input path — but once for the whole cohort.
    # The inputs share one output directory, so the disk-space estimate has to
    # be the sum: checking each input separately would pass N times over on a
    # filesystem that only has room for one.
    for p in input_paths:
        preflight.require(
            _annotate_dependencies(index_type=manifest.index.type, input_path=p, bgzip=bgzip),
            context=f"annotate against {db_id_resolved}",
        )
    _cpus.warn_if_oversubscribed(threads, what=f"annotate against {db_id_resolved}")
    total_bases = sum(estimate_input_bases(p) for p in input_paths)
    needed_bytes = estimate_output_bytes(
        input_bases=total_bases,
        n_feature_sets=len(requested),
        keep_presmoothed=keep_presmoothed,
        smooth=smooth,
    )
    logger.info(
        "estimated output footprint: %s for %d feature set(s) over ~%.2f Gbp "
        "of input across %d file(s)",
        _human_bytes(needed_bytes),
        len(requested),
        total_bases / 1e9,
        len(input_paths),
    )
    diskspace.require_free_space(
        output_dir,
        needed_bytes,
        what=(f"annotating {len(input_paths)} input(s) ({len(requested)} feature set(s))"),
        estimated=True,
        hint=(
            "Note that --bgzip shrinks the final output but not the peak: every "
            "BED is written in full before the compression pass starts.\n"
            "Options:\n"
            "  - write to a larger filesystem with --outdir\n"
            "  - annotate fewer inputs or feature sets at a time\n"
            "  - drop one of the two outputs with --no-keep-presmoothed or --no-smooth\n"
            "  - pass --no-space-check if this estimate looks wrong for your input"
        ),
        skip=not check_space,
    )

    n_threads = _cpus.resolve_threads(threads)
    progress.start(
        f"Annotating {len(input_paths)} inputs against {db_id_resolved}",
        f"{len(requested)} feature set(s), {n_threads} thread(s), "
        f"~{_human_bytes(needed_bytes)} estimated output",
        "one index load per feature set for the whole batch",
    )

    t_start = time.perf_counter()

    # Per-input prefixes and output-path dicts (same naming as single-input).
    prefixes: dict[Path, str] = {}
    presmoothed_by_input: dict[Path, dict[str, Path]] = {}
    smoothed_by_input: dict[Path, dict[str, Path]] = {}
    for p in input_paths:
        prefix = f"{_derive_input_basename(p)}.{db_id_resolved}"
        if k is not None:
            prefix = f"{prefix}.k{query_k}"
        prefixes[p] = prefix
        presmoothed_by_input[p] = (
            {fs: output_dir / f"{prefix}.{fs}.presmoothed.bed" for fs in requested}
            if keep_presmoothed
            else {}
        )
        smoothed_by_input[p] = (
            {fs: output_dir / f"{prefix}.{fs}.smoothed.bed" for fs in requested} if smooth else {}
        )

    _run_hks_backend_batch(
        manifest=manifest,
        db_dir=db_dir,
        input_paths=input_paths,
        prefixes=prefixes,
        output_dir=output_dir,
        requested=requested,
        smooth=smooth,
        keep_presmoothed=keep_presmoothed,
        presmoothed_by_input=presmoothed_by_input,
        smoothed_by_input=smoothed_by_input,
        threads=threads,
        k=query_k,
        progress=progress,
    )

    # bgzip per input (mutating each input's path dicts in place), then build results.
    results: dict[Path, AnnotateResult] = {}
    for p in input_paths:
        pre = presmoothed_by_input[p]
        smo = smoothed_by_input[p]
        if bgzip:
            for fs in requested:
                if fs in pre:
                    pre[fs] = bgzip_file(pre[fs], threads=threads)
                if fs in smo:
                    smo[fs] = bgzip_file(smo[fs], threads=threads)
        results[p] = AnnotateResult(
            presmoothed_paths=pre, smoothed_paths=smo, combined_intermediate=None
        )

    logger.info(
        "annotate batch complete in %.1fs (%d input(s))",
        time.perf_counter() - t_start,
        len(input_paths),
    )
    return results
