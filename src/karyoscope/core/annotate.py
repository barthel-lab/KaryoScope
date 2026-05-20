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

import logging
import multiprocessing as mp
import os
import sys
import threading
from collections.abc import Iterator
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
from karyoscope.core.io.kmc import run_get_featureids
from karyoscope.core.smooth import (
    HierarchyIndex,
    chunked_seq_reader,
    make_features_for_worker,
    process_seq_chunk,
    worker_initializer,
)
from karyoscope.exceptions import (
    DatabaseNotFoundError,
    KaryoscopeError,
)
from karyoscope.manifest import validate_database_layout

logger = logging.getLogger(__name__)


#: FASTA extensions we recognise when deriving the output basename.
#: Order matters — longer extensions first so they win over shorter ones.
_FASTA_EXTENSIONS: tuple[str, ...] = (
    ".fasta.gz",
    ".fa.gz",
    ".fna.gz",
    ".fasta",
    ".fa",
    ".fna",
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
    """Strip recognised FASTA extensions from a path's filename.

    ``my_assembly.fa.gz`` -> ``my_assembly``. Falls back to the raw
    stem if no known extension matches.
    """
    name = input_path.name
    name_lower = name.lower()
    for ext in _FASTA_EXTENSIONS:
        if name_lower.endswith(ext):
            return name[: -len(ext)]
    return input_path.stem


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


def _smooth_one_feature_set(
    *,
    combined_bed: Path,
    feature_set: str,
    features: Features,
    index: HierarchyIndex,
    presmoothed_path: Path | None,
    smoothed_path: Path | None,
    threads: int,
    chunk_size: int = 50000,
    preserve_input_order: bool = True,
) -> None:
    """Stream ``combined_bed`` through a worker pool, writing per-set BEDs.

    ``presmoothed_path`` and ``smoothed_path`` are independently
    optional -- pass ``None`` for either to skip writing it. At least
    one must be set; the caller is responsible for not calling this
    function with both ``None``.

    The worker pool is sized at ``threads`` (or auto via
    ``os.cpu_count()`` when ``threads <= 0``). Workers are initialised
    once with the shared :class:`HierarchyIndex` and a per-set
    :class:`FeaturesForWorker` projection; each chunk arrives as a
    list of combined-BED lines that the worker translates, smooths,
    and returns keyed per-sequence.

    Two output codepaths, selected by ``preserve_input_order``:

    * **True (default)** -- *assembly mode*. Each worker's
      per-sequence results stream into a per-sequence temp file in a
      :class:`tempfile.TemporaryDirectory`. After the pool drains,
      the temp files are concatenated in input-sequence order to
      produce the final BED. ``pool.imap_unordered`` is used so the
      slowest chunk (chr1, typically) doesn't queue up GB-scale
      results from finished chunks in the IPC pipe. Output is
      byte-identical to the old ``pool.imap`` codepath; memory peak
      drops to ~ workers x one chunk's worth.

    * **False** -- *reads mode*. ``pool.imap_unordered`` straight
      through to the output file, no temp files. Use for inputs
      with millions of small sequences (long-read FASTA, future
      FASTQ/BAM) where per-sequence temp files aren't feasible
      (file-descriptor limits, syscall overhead) and the user
      doesn't care about input order in the output BED.
    """
    if presmoothed_path is None and smoothed_path is None:
        raise KaryoscopeError("internal error: _smooth_one_feature_set called with no output paths")

    pool_size = threads if threads > 0 else (os.cpu_count() or 1)
    features_for_worker = make_features_for_worker(features, feature_set)

    ctx = mp.get_context("spawn")
    # Inherit the main process's root log level so worker INFO lines
    # (per-sequence smoothing progress) appear when the user passed
    # -v, but stay silent at the default WARNING level.
    worker_log_level = logging.getLogger().level

    if preserve_input_order:
        _smooth_with_per_sequence_tempfiles(
            combined_bed=combined_bed,
            presmoothed_path=presmoothed_path,
            smoothed_path=smoothed_path,
            ctx=ctx,
            pool_size=pool_size,
            index=index,
            features_for_worker=features_for_worker,
            worker_log_level=worker_log_level,
            chunk_size=chunk_size,
        )
    else:
        _smooth_streaming_unordered(
            combined_bed=combined_bed,
            presmoothed_path=presmoothed_path,
            smoothed_path=smoothed_path,
            ctx=ctx,
            pool_size=pool_size,
            index=index,
            features_for_worker=features_for_worker,
            worker_log_level=worker_log_level,
            chunk_size=chunk_size,
        )


def _smooth_streaming_unordered(
    *,
    combined_bed: Path,
    presmoothed_path: Path | None,
    smoothed_path: Path | None,
    ctx,
    pool_size: int,
    index: HierarchyIndex,
    features_for_worker,
    worker_log_level: int,
    chunk_size: int,
) -> None:
    """Reads-mode smoothing: ``imap_unordered`` direct to output files.

    Sequences appear in the output BED in the order their workers
    finished, not input order. Acceptable when the caller has set
    ``preserve_input_order=False`` (i.e. input is sequencing reads,
    not an assembly, and order is irrelevant downstream).
    """
    pre_handle = presmoothed_path.open("w") if presmoothed_path is not None else None
    smo_handle = smoothed_path.open("w") if smoothed_path is not None else None
    try:
        with ctx.Pool(
            processes=pool_size,
            initializer=worker_initializer,
            initargs=(index, features_for_worker, worker_log_level),
        ) as pool:
            stop_event = _spawn_pool_watchdog(pool)
            try:
                for chunk_result in pool.imap_unordered(
                    process_seq_chunk,
                    chunked_seq_reader(combined_bed, chunk_size),
                ):
                    # chunk_result: {seq_name: (smoothed_lines, presmoothed_lines)}
                    for _seq, (smo_lines, pre_lines) in chunk_result.items():
                        if pre_handle is not None:
                            pre_handle.writelines(pre_lines)
                        if smo_handle is not None:
                            smo_handle.writelines(smo_lines)
            finally:
                stop_event.set()
    finally:
        if pre_handle is not None:
            pre_handle.close()
        if smo_handle is not None:
            smo_handle.close()


def _smooth_with_per_sequence_tempfiles(
    *,
    combined_bed: Path,
    presmoothed_path: Path | None,
    smoothed_path: Path | None,
    ctx,
    pool_size: int,
    index: HierarchyIndex,
    features_for_worker,
    worker_log_level: int,
    chunk_size: int,
) -> None:
    """Assembly-mode smoothing: per-sequence temp files + concat in input order.

    ``pool.imap_unordered`` keeps memory bounded (no GB-scale result
    queue waiting on the slowest chunk to finish). Per-sequence
    output streams to its own temp file; at the end we concatenate
    the temp files in input-sequence order. Output is byte-identical
    to the legacy ``pool.imap`` codepath.

    Input sequence order is captured at dispatch time: we wrap
    :func:`chunked_seq_reader` with a peek-and-yield generator that
    records each chunk's sequence names before yielding the chunk
    to the pool.
    """
    import shutil
    import tempfile

    # Temp dir lives next to the output so cleanup happens on the
    # same filesystem (no cross-device shutil.copy cost) and so the
    # temp dir is on /scratch when the output is.
    anchor = smoothed_path if smoothed_path is not None else presmoothed_path
    assert anchor is not None  # guarded by the caller check above
    out_dir = anchor.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    input_seq_order: list[str] = []
    seen_seqs: set[str] = set()

    def _record_and_yield():
        """Peek each chunk for sequence names, record in input order, yield."""
        for chunk in chunked_seq_reader(combined_bed, chunk_size):
            for raw in chunk:
                seq = raw.partition("\t")[0]
                if seq and seq not in seen_seqs:
                    seen_seqs.add(seq)
                    input_seq_order.append(seq)
            yield chunk

    # Per-sequence temp file paths, opened lazily on first write.
    with tempfile.TemporaryDirectory(prefix="ks_smooth_", dir=out_dir) as tmpdir_str:
        tmpdir = Path(tmpdir_str)
        smo_temps: dict[str, Path] = {}
        pre_temps: dict[str, Path] = {}

        with ctx.Pool(
            processes=pool_size,
            initializer=worker_initializer,
            initargs=(index, features_for_worker, worker_log_level),
        ) as pool:
            stop_event = _spawn_pool_watchdog(pool)
            try:
                for chunk_result in pool.imap_unordered(
                    process_seq_chunk,
                    _record_and_yield(),
                ):
                    for seq, (smo_lines, pre_lines) in chunk_result.items():
                        if smoothed_path is not None and smo_lines:
                            if seq not in smo_temps:
                                smo_temps[seq] = tmpdir / f"{len(smo_temps):08d}.smo"
                            with smo_temps[seq].open("a") as f:
                                f.writelines(smo_lines)
                        if presmoothed_path is not None and pre_lines:
                            if seq not in pre_temps:
                                pre_temps[seq] = tmpdir / f"{len(pre_temps):08d}.pre"
                            with pre_temps[seq].open("a") as f:
                                f.writelines(pre_lines)
            finally:
                stop_event.set()

        # Concatenate temp files in input sequence order.
        if smoothed_path is not None:
            with smoothed_path.open("wb") as out_h:
                for seq in input_seq_order:
                    p = smo_temps.get(seq)
                    if p is None:
                        continue
                    with p.open("rb") as src:
                        shutil.copyfileobj(src, out_h)
        if presmoothed_path is not None:
            with presmoothed_path.open("wb") as out_h:
                for seq in input_seq_order:
                    p = pre_temps.get(seq)
                    if p is None:
                        continue
                    with p.open("rb") as src:
                        shutil.copyfileobj(src, out_h)
        # tmpdir cleaned up automatically on context manager exit


def _bgzip_file(path: Path) -> Path:
    """Compress ``path`` in-place with ``bgzip``, returning the new path.

    ``bgzip`` removes the source file by default (matches gzip's behaviour).
    Returns ``Path(str(path) + ".gz")``.
    """
    bgzip = require_tool(
        "bgzip",
        install_hint="Install htslib (`conda install -c bioconda htslib`), "
        "or rerun with --no-bgzip to skip compression.",
    )
    run_tool([bgzip, "-f", str(path)])
    return Path(str(path) + ".gz")


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
        If ``get_featureIDs`` (or ``bgzip``, when requested) isn't found.
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

    # Run the C++ helper.
    input_basename = _derive_input_basename(input_path)
    prefix = f"{input_basename}.{db_id_resolved}"
    kmc_db_basename = db_dir / manifest.index.basename
    combined_bed = run_get_featureids(
        db_path=kmc_db_basename,
        input_path=input_path,
        output_dir=output_dir,
        threads=threads,
        prefix=prefix,
        capture=True,
    )
    if not combined_bed.is_file():
        raise KaryoscopeError(f"get_featureIDs did not produce expected output at {combined_bed}")
    logger.debug("combined BED at %s", combined_bed)

    # Compute output paths (uncompressed names; bgzip later if requested).
    presmoothed_paths: dict[str, Path] = (
        {fs: output_dir / f"{prefix}.{fs}.presmoothed.bed" for fs in requested}
        if keep_presmoothed
        else {}
    )
    smoothed_paths: dict[str, Path] = (
        {fs: output_dir / f"{prefix}.{fs}.smoothed.bed" for fs in requested} if smooth else {}
    )

    # Run the per-feature-set pass. When smoothing is on we go through
    # the worker pool; when smoothing is off we use the simpler
    # in-process splitter (no need to fork workers for a one-line
    # translation).
    if smooth:
        logger.debug(
            "smoothing pass with threads=%d, preserve_input_order=%s",
            threads,
            preserve_input_order,
        )
        for fs in requested:
            _smooth_one_feature_set(
                combined_bed=combined_bed,
                feature_set=fs,
                features=features,
                index=indices[fs],
                presmoothed_path=presmoothed_paths.get(fs),
                smoothed_path=smoothed_paths.get(fs),
                threads=threads,
                preserve_input_order=preserve_input_order,
            )
    else:
        # Only presmoothed output, no smoothing.
        _split_combined_bed(combined_bed, requested, features, presmoothed_paths)

    # bgzip (or not).
    if bgzip:
        for fs in requested:
            if fs in presmoothed_paths:
                presmoothed_paths[fs] = _bgzip_file(presmoothed_paths[fs])
            if fs in smoothed_paths:
                smoothed_paths[fs] = _bgzip_file(smoothed_paths[fs])

    # Tidy up the combined intermediate unless asked to keep it.
    if not keep_intermediates:
        try:
            combined_bed.unlink()
            combined_kept: Path | None = None
            logger.debug("removed combined intermediate %s", combined_bed)
        except OSError as e:
            logger.warning("could not remove intermediate %s: %s", combined_bed, e)
            combined_kept = combined_bed
    else:
        combined_kept = combined_bed

    return AnnotateResult(
        presmoothed_paths=presmoothed_paths,
        smoothed_paths=smoothed_paths,
        combined_intermediate=combined_kept,
    )
