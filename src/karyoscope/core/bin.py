"""Aggregate base-pair annotation BED records into fixed-size bins.

This module is the ``karyoscope bin`` building block. It takes a sorted
per-base-pair (or run-length-encoded) BED with one feature label per
interval and emits a coarser-grained BED whose intervals are
``bin_size``-bp windows, each labelled with the feature that wins the
within-bin selection rule below.

Selection rule (faithful port of the archive's ``bin_features.py``):

1. **Sum overlap per feature**: for each bin, accumulate the bp of
   overlap contributed by every feature label seen in that bin.
2. **Leaf prioritisation (optional)**: if a leaf-feature set was
   provided and any leaf feature has overlap in the bin, only leaf
   features compete. Falls back to all features when no leaf is
   present. Leaves come from the database's :mod:`hierarchy.tsv` for
   the relevant feature set: any node that is a child of some other
   node but is itself never a parent.
3. **Novel deprioritisation**: ``novel`` only wins if its overlap is
   strictly greater than half the bin. Otherwise the best non-novel
   feature wins. This stops bins with fragmented real coverage from
   being incorrectly called novel just because the unannotated bases
   form a plurality.
4. **Tie-break**: when two features have equal overlap, the
   alphabetically lower feature wins, with non-novel beating novel.

Bins are written in order. Adjacent bins that share a chromosome and
label are merged into a single output interval. The final output is
therefore coordinate-sorted within each chromosome and contains no two
adjacent same-label rows.

Differences from the archive port:

* The legacy ``--specific`` suffix path is dropped. v0.1's
  ``hierarchy.tsv`` makes leaf detection structural (a child that is
  never a parent), so the old "anything ending in ``_specific``" hack
  is no longer needed. Callers that want leaf prioritisation pass a
  ``leaf_set`` explicitly (typically built via :func:`leaves_for`).
* The module exposes :func:`bin_features` as a clean library entry
  point so :mod:`karyoscope.core.scaffold` can call it in-process
  rather than shelling out.
"""

from __future__ import annotations

import gzip
import logging
import multiprocessing as mp
import os
from collections import defaultdict
from collections.abc import Callable, Iterable, Iterator
from pathlib import Path
from typing import IO

from karyoscope.core.io.features import NOVEL_NAME
from karyoscope.core.io.hierarchy import Hierarchy
from karyoscope.core.smooth import chunked_seq_reader
from karyoscope.exceptions import BinError

logger = logging.getLogger(__name__)

#: Default lines-per-chunk for the worker pool. Matches the smoothing
#: module's default; chunks always end at a sequence boundary, so the
#: actual chunk size is "at least N lines and ends at the next
#: sequence boundary".
DEFAULT_CHUNK_SIZE = 50_000


#: Features that lose ties and that win only when they cover more than
#: half the bin. Currently just the novel sentinel; centralised so it
#: stays in sync with :data:`karyoscope.core.io.features.NOVEL_NAME`.
_DEPRIORITIZED: frozenset[str] = frozenset({NOVEL_NAME})


# --- leaf-set helper ------------------------------------------------


def leaves_for(hierarchy: Hierarchy, feature_set: str) -> set[str]:
    """Return the leaf nodes of ``feature_set`` in ``hierarchy``.

    A leaf is any child that never appears as a parent within the same
    feature set. Empty for vacuous feature sets; callers may want to
    treat an empty result the same as ``None`` (no leaf prioritisation).
    """
    children: set[str] = set()
    parents: set[str] = set()
    for row in hierarchy.rows_in(feature_set):
        children.add(row.child)
        parents.add(row.parent)
    return children - parents


# --- selection rule -------------------------------------------------


def _pick_winner(candidates: dict[str, int], bin_size: int) -> str:
    """Pick the winning feature from a ``{feature: overlap_bp}`` map.

    Encodes the novel-deprioritisation and tie-break rules. Callers
    must pass a non-empty ``candidates``; ``KeyError`` otherwise.
    """
    best_overlap = max(candidates.values())
    leaders = [f for f, v in candidates.items() if v == best_overlap]
    # Non-deprioritized wins ties; then alphabetical order.
    leaders.sort(key=lambda f: (f in _DEPRIORITIZED, f))
    winner = leaders[0]

    if winner in _DEPRIORITIZED and candidates[winner] <= bin_size / 2:
        non_depr = {f: v for f, v in candidates.items() if f not in _DEPRIORITIZED}
        if non_depr:
            best_v = max(non_depr.values())
            return sorted(f for f, v in non_depr.items() if v == best_v)[0]
    return winner


def _best_feature(counts: dict[str, int], leaf_set: set[str] | None, bin_size: int) -> str:
    """Pick the best feature for a bin, applying leaf prioritisation."""
    if leaf_set is None:
        return _pick_winner(counts, bin_size)
    leaves = {k: v for k, v in counts.items() if k in leaf_set}
    if leaves:
        return _pick_winner(leaves, bin_size)
    return _pick_winner(counts, bin_size)


# --- merger ---------------------------------------------------------


class _Coalescer:
    """Forward (chrom,start,end,feature) records to a sink, merging adjacent same-label rows.

    Adjacent here means abutting (``prev.end == next.start``) on the
    same chromosome with the same feature label.
    """

    __slots__ = ("_chrom", "_end", "_feature", "_sink", "_start")

    def __init__(self, sink: Callable[[str, int, int, str], None]) -> None:
        self._sink = sink
        self._chrom: str | None = None
        self._start = -1
        self._end = -1
        self._feature: str | None = None

    def add(self, chrom: str, start: int, end: int, feature: str) -> None:
        if chrom == self._chrom and start == self._end and feature == self._feature:
            self._end = end
            return
        self.flush()
        self._chrom = chrom
        self._start = start
        self._end = end
        self._feature = feature

    def flush(self) -> None:
        if self._chrom is None:
            return
        self._sink(self._chrom, self._start, self._end, self._feature)  # type: ignore[arg-type]
        self._chrom = None


# --- core driver ----------------------------------------------------


def _drive(
    records: Iterable[tuple[str, int, int, str]],
    bin_size: int,
    leaf_set: set[str] | None,
    sink: Callable[[str, int, int, str], None],
) -> None:
    """Stream binning loop. Calls ``sink`` once per output row."""
    coalescer = _Coalescer(sink)
    active: dict[int, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    current_chrom: str | None = None
    current_max_end = 0

    def _flush_below(limit_idx: int) -> None:
        for b_idx in sorted(k for k in active if k < limit_idx):
            best = _best_feature(active[b_idx], leaf_set, bin_size)
            start = b_idx * bin_size
            end = min(start + bin_size, current_max_end)
            if end > start:
                coalescer.add(current_chrom, start, end, best)  # type: ignore[arg-type]
            del active[b_idx]

    def _flush_all() -> None:
        if not active or current_chrom is None:
            return
        for b_idx in sorted(active):
            best = _best_feature(active[b_idx], leaf_set, bin_size)
            start = b_idx * bin_size
            end = min(start + bin_size, current_max_end)
            if end > start:
                coalescer.add(current_chrom, start, end, best)
        active.clear()

    for chrom, start, end, feature in records:
        if chrom != current_chrom:
            _flush_all()
            coalescer.flush()
            current_chrom = chrom
            current_max_end = 0

        current_max_end = max(current_max_end, end)

        start_bin = start // bin_size
        end_bin = (end - 1) // bin_size

        _flush_below(start_bin)

        for b_idx in range(start_bin, end_bin + 1):
            bin_start = b_idx * bin_size
            bin_end = bin_start + bin_size
            ov = min(end, bin_end) - max(start, bin_start)
            if ov > 0:
                active[b_idx][feature] += ov

    _flush_all()
    coalescer.flush()


# --- public APIs ----------------------------------------------------


def bin_records(
    records: Iterable[tuple[str, int, int, str]],
    *,
    bin_size: int,
    leaf_set: set[str] | None = None,
) -> Iterator[tuple[str, int, int, str]]:
    """Bin an iterable of BED records and yield binned, merged records.

    Pure function — no I/O. Used directly by unit tests and by
    :func:`bin_features` (which wraps it with file handling).

    ``records`` must be grouped by chromosome and sorted by start
    position within each chromosome (the order produced by
    :mod:`karyoscope.core.annotate`). Bins are emitted in coordinate
    order per chromosome, and adjacent bins with the same winning
    feature are coalesced.
    """
    if bin_size < 1:
        raise BinError(f"bin_size must be a positive integer, got {bin_size}")

    out: list[tuple[str, int, int, str]] = []
    _drive(records, bin_size, leaf_set, lambda c, s, e, f: out.append((c, s, e, f)))
    yield from out


def _open_in(path: Path) -> IO[str]:
    """Open ``path`` for text reading.

    Returns ``sys.stdin`` when ``path`` is ``-``. Callers must not
    close the returned handle in that case; use :func:`_close_if_owned`
    or check ``handle is sys.stdin`` before closing.
    """
    if str(path) == "-":
        import sys

        return sys.stdin
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def _open_out(path: Path, *, gzip_out: bool) -> IO[str]:
    """Open ``path`` for text writing.

    Returns ``sys.stdout`` when ``path`` is ``-``. Same close-semantics
    caveat as :func:`_open_in`.
    """
    if str(path) == "-":
        import sys

        return sys.stdout
    if gzip_out:
        return gzip.open(path, "wt")
    return path.open("w")


def _close_if_owned(handle: IO[str]) -> None:
    """Close ``handle`` unless it's stdin/stdout (which we don't own)."""
    import sys

    if handle is sys.stdin or handle is sys.stdout:
        return
    handle.close()


def _iter_bed(handle: IO[str]) -> Iterator[tuple[str, int, int, str]]:
    """Yield ``(chrom, start, end, feature)`` from a 4-col BED stream.

    Skips blank lines. Raises :class:`BinError` with the line number on
    malformed rows so the user can find the bad input.
    """
    for i, raw in enumerate(handle, start=1):
        line = raw.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            raise BinError(
                f"line {i}: expected 4+ tab-separated columns, got {len(parts)}: {raw!r}"
            )
        try:
            start = int(parts[1])
            end = int(parts[2])
        except ValueError as e:
            raise BinError(f"line {i}: non-integer coordinates: {raw!r}") from e
        yield parts[0], start, end, parts[3]


def _resolve_pool_size(threads: int) -> int:
    """Translate the user-facing ``threads`` arg into an actual pool size.

    ``threads <= 0`` means "auto" (``os.cpu_count()``). ``threads == 1``
    is single-threaded (caller short-circuits the pool entirely).
    """
    if threads <= 0:
        return os.cpu_count() or 1
    return threads


# --- worker-pool dispatch ------------------------------------------

# Per-worker globals, set once via :func:`worker_initializer`. Pickling
# cost for ``bin_size`` and ``leaf_set`` is paid once per worker
# (when the Pool spawns it), not once per chunk.
_worker_bin_size: int | None = None
_worker_leaf_set: frozenset[str] | None = None


def worker_initializer(
    bin_size: int,
    leaf_set: frozenset[str] | None,
) -> None:
    """Initialise a :class:`multiprocessing.Pool` worker.

    Public so it can be the ``initializer=`` callback. Not part of the
    user-facing API; tests may reference it.
    """
    global _worker_bin_size, _worker_leaf_set
    _worker_bin_size = bin_size
    _worker_leaf_set = leaf_set


def process_seq_chunk(chunk: list[str]) -> list[str]:
    """Bin one chunk of BED lines into a list of output lines.

    The chunk must contain complete sequences (no fragment ever
    crosses a chunk boundary -- :func:`karyoscope.core.smooth.chunked_seq_reader`
    enforces this). Returns newline-terminated output lines ready for
    the caller to ``writelines`` into the output file.

    Output order within the chunk matches input order: contigs appear
    in their input order, and within each contig the bins are
    coordinate-sorted.
    """
    if _worker_bin_size is None:
        raise BinError("worker_initializer was not called; cannot process chunk")
    bin_size = _worker_bin_size
    leaf_set: set[str] | None = set(_worker_leaf_set) if _worker_leaf_set else None

    out_lines: list[str] = []

    def _sink(c: str, s: int, e: int, f: str) -> None:
        out_lines.append(f"{c}\t{s}\t{e}\t{f}\n")

    def _records() -> Iterator[tuple[str, int, int, str]]:
        for i, raw in enumerate(chunk, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                raise BinError(
                    f"chunk line {i}: expected 4+ tab-separated columns, got {len(parts)}: {raw!r}"
                )
            try:
                yield parts[0], int(parts[1]), int(parts[2]), parts[3]
            except ValueError as e:
                raise BinError(f"chunk line {i}: non-integer coordinates: {raw!r}") from e

    _drive(_records(), bin_size, leaf_set, _sink)
    return out_lines


def bin_features(
    input_path: Path,
    output_path: Path,
    *,
    bin_size: int,
    leaf_set: set[str] | None = None,
    gzip_out: bool | None = None,
    threads: int = 1,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> None:
    """Bin a BED file and write the result.

    ``input_path`` may be plain or gzipped (``.gz``). When
    ``gzip_out`` is ``None`` (the default) the output is gzipped iff
    ``output_path`` ends in ``.gz``. Use ``Path("-")`` for either
    path to mean stdin/stdout (always single-threaded in that case --
    the pool path requires a real on-disk input).

    Parallelism: ``threads`` (default 1) controls the worker-pool
    size. ``threads == 0`` means "auto" (``os.cpu_count()``).
    ``threads == 1`` short-circuits the pool entirely and runs the
    in-process binner directly -- output is byte-for-byte identical
    to the multi-threaded path.

    The pool divides the input into chunks of at least ``chunk_size``
    lines, never splitting a sequence across chunks. Each chunk goes
    to one worker. ``pool.imap`` preserves input order, so the output
    BED is grouped per sequence in source order with bins
    coordinate-sorted within each sequence (the same property the
    single-threaded path produces).
    """
    if gzip_out is None:
        gzip_out = str(output_path).endswith(".gz") and str(output_path) != "-"

    is_stdio = str(input_path) == "-" or str(output_path) == "-"
    pool_size = _resolve_pool_size(threads)
    use_pool = pool_size > 1 and not is_stdio

    logger.info(
        "binning %s -> %s (bin_size=%d, leaf_set=%s, threads=%s)",
        input_path,
        output_path,
        bin_size,
        "yes" if leaf_set else "no",
        pool_size if use_pool else 1,
    )

    if not use_pool:
        # Single-threaded path: stream directly through the binner.
        in_h = _open_in(input_path)
        out_h = _open_out(output_path, gzip_out=gzip_out)
        try:

            def _emit(c: str, s: int, e: int, f: str) -> None:
                out_h.write(f"{c}\t{s}\t{e}\t{f}\n")

            _drive(_iter_bed(in_h), bin_size, leaf_set, _emit)
        finally:
            _close_if_owned(in_h)
            _close_if_owned(out_h)
        return

    # Multi-threaded path: chunked-by-sequence-boundary dispatch
    # through a process pool. ``chunked_seq_reader`` lives in
    # :mod:`karyoscope.core.smooth`; the sequence-boundary guarantee
    # is the same one smoothing relies on, so per-sequence binning
    # semantics are preserved.
    leaf_frozen = frozenset(leaf_set) if leaf_set else None
    ctx = mp.get_context("spawn")
    out_h = _open_out(output_path, gzip_out=gzip_out)
    try:
        with ctx.Pool(
            processes=pool_size,
            initializer=worker_initializer,
            initargs=(bin_size, leaf_frozen),
        ) as pool:
            for out_lines in pool.imap(
                process_seq_chunk,
                chunked_seq_reader(input_path, chunk_size),
            ):
                out_h.writelines(out_lines)
    finally:
        _close_if_owned(out_h)
