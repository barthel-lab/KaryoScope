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

from karyoscope import installed as _installed
from karyoscope.core.external import require_tool, run_tool
from karyoscope.core.io.features import Features, parse_features, render_feature
from karyoscope.core.io.hierarchy import (
    Hierarchy,
    parse_hierarchy,
    validate_hierarchy,
)
from karyoscope.core.io.hks import (
    convert_hks_tsv_to_bed,
    run_hks_lookup,
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
                "Install one with `karyoscope download`, or pass --db <ID>."
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
                for fs in feature_sets:
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

    pool_size = threads if threads > 0 else (os.cpu_count() or 1)
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


def _bgzip_file(path: Path, threads: int = 1) -> Path:
    """Compress ``path`` in-place with ``bgzip``, returning the new path.

    ``bgzip`` removes the source file by default (matches gzip's behaviour).
    Returns ``Path(str(path) + ".gz")``. Logs per-file start + completion
    at INFO so a long bgzip pass (12 files for a 6-feature-set human
    database) doesn't look like the pipeline has hung.

    ``threads`` is forwarded as ``bgzip -@``; the htslib bgzip compresses
    a single file in parallel when given more than one thread. We
    process files sequentially within the bgzip pass, so passing the
    user's full ``--threads`` here is the right call (no contention
    with concurrent file compressions). ``threads=1`` (the default)
    omits ``-@`` entirely for cleanest subprocess invocation.
    """
    bgzip = require_tool(
        "bgzip",
        install_hint="Install htslib (`conda install -c bioconda htslib`), "
        "or rerun with --no-bgzip to skip compression.",
    )
    orig_size = path.stat().st_size
    logger.info("bgzipping %s (%s, threads=%d)", path.name, _human_bytes(orig_size), threads)
    t0 = time.perf_counter()
    cmd = [bgzip, "-f"]
    if threads > 1:
        cmd.extend(["-@", str(threads)])
    cmd.append(str(path))
    run_tool(cmd)
    out_path = Path(str(path) + ".gz")
    out_size = out_path.stat().st_size if out_path.is_file() else 0
    dt = time.perf_counter() - t0
    logger.info(
        "bgzipped %s (%s -> %s) in %.1fs",
        out_path.name,
        _human_bytes(orig_size),
        _human_bytes(out_size),
        dt,
    )
    return out_path


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
) -> None:
    """Run the HKS lookup and optional smoothing for every requested feature set.

    Unlike the KMC backend (one combined query, then translate integer feature
    ids), HKS queries one ``.hksf`` per feature set and reads label names
    directly. Each feature set is processed independently:

    1. ``hks lookup`` against ``<basename>.<fs>.hksf`` -> a raw TSV.
    2. If ``keep_presmoothed``: convert the raw TSV to the presmoothed BED.
    3. If ``smooth``: ``hks smooth`` (using ``<basename>.<fs>.hierarchy.txt``)
       -> the smoothed BED.

    The raw TSV is a per-feature-set temp file, deleted after each set.
    """
    base_path = db_dir / (manifest.index.basename + ".hksb")
    k = manifest.kmer.size

    t_hks_start = time.perf_counter()
    for fs in requested:
        fs_file = db_dir / f"{manifest.index.basename}.{fs}.hksf"
        hierarchy_file = db_dir / f"{manifest.index.basename}.{fs}.hierarchy.txt"
        raw_tsv = output_dir / f"{prefix}.{fs}.lookup_raw.tmp.tsv"

        logger.info(
            "running hks lookup for feature set %r on %s (threads=%d)",
            fs,
            input_path.name,
            threads,
        )
        run_hks_lookup(
            base_path=base_path,
            feature_set_file=fs_file,
            k=k,
            input_path=input_path,
            output_path=raw_tsv,
            threads=threads,
            capture=True,
        )
        if not raw_tsv.is_file():
            raise KaryoscopeError(f"hks lookup did not produce expected output at {raw_tsv}")

        try:
            if keep_presmoothed:
                convert_hks_tsv_to_bed(raw_tsv, presmoothed_paths[fs])

            if smooth:
                t_smo = time.perf_counter()
                logger.info("running hks smooth for feature set %r", fs)
                run_hks_smooth(
                    hierarchy_file=hierarchy_file,
                    input_path=raw_tsv,
                    output_path=smoothed_paths[fs],
                    threads=threads,
                    capture=True,
                )
                logger.info("smoothed feature set %r in %.1fs", fs, time.perf_counter() - t_smo)
        finally:
            try:
                raw_tsv.unlink()
            except OSError as exc:
                logger.warning("could not remove temp lookup TSV %s: %s", raw_tsv, exc)

    logger.info(
        "hks backend complete in %.1fs (%d feature set(s))",
        time.perf_counter() - t_hks_start,
        len(requested),
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
        ``0`` means auto (``os.cpu_count()``).
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

    Raises
    ------
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
    t_annotate_start = time.perf_counter()

    # Parse features.tsv up front so we fail fast if it's malformed.
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
        # Build the features-columns map for the cross-validation check.
        feature_columns: dict[str, set[str]] = {
            fs: {row[fs] for row in features.table.values()} for fs in requested
        }
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
        else:
            # Only presmoothed output, no smoothing.
            logger.info("splitting combined BED into %d per-feature-set BED(s)", len(requested))
            t_split_start = time.perf_counter()
            _split_combined_bed(combined_bed, requested, features, presmoothed_paths)
            logger.info("split complete in %.1fs", time.perf_counter() - t_split_start)

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
                presmoothed_paths[fs] = _bgzip_file(presmoothed_paths[fs], threads=threads)
            if fs in smoothed_paths:
                smoothed_paths[fs] = _bgzip_file(smoothed_paths[fs], threads=threads)
        logger.info("bgzip pass complete in %.1fs", time.perf_counter() - t_bgzip_start)

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
