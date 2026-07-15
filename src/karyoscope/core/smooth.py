"""Hierarchy-aware smoothing of per-feature-set BED tracks.

This module implements the "smoothing" step of the annotate pipeline.
The intuition: when a sequence is annotated by querying every k-mer
against the database, the result is a per-position label that can be
noisy — short novel runs (k-mers missing from the index) appear inside
otherwise-uniform feature blocks, and labels can flicker between
sibling features when the k-mer-level evidence is mixed. Smoothing
uses the database's feature hierarchy to clean this up: when two
flanking intervals share a common ancestor in the hierarchy, the noisy
intervals between them are *promoted* to that ancestor, generalising
from "I don't know exactly what this is" to "I know it's at least
this".

Algorithm summary (faithful port of ``smooth_features.py`` with the
``phylogeny → hierarchy`` rename and the ``_specific``-suffix logic
removed):

1. For each sequence (independently), build a list of intervals from
   its rows in the presmoothed BED.
2. Repeatedly scan the list with a moving window. At each position,
   look forward and backward to find the largest window of intervals
   whose features are mutually ancestor-compatible (each is an
   ancestor of the previous, with bounded gap distance), and identify
   any unrelated interval that bounds the window.
3. Compute the LCA of the window's left and right features. Promote
   any interval within the window whose feature is an ancestor of the
   LCA (i.e., less specific) up to the LCA itself.
4. Repeat the entire pass until no interval is promoted (fixed point).

Two sentinels at output time:

* Feature ids of 0 are rendered as ``"novel"`` upstream of this module
  (in :func:`karyoscope.core.io.features.render_feature`); their
  hierarchy label internally is the root (:data:`REQUIRED_ROOT`,
  ``"categorized"``). When emitting BED records we map any interval
  still labelled with the root *and* flagged as novel back to
  ``"novel"``, so users never see the internal sentinel.
* Other intervals can also legitimately end up promoted to the root by
  smoothing (LCA of two top-level siblings is the root). We keep
  those as ``"categorized"`` in the output — the user explicitly
  asked for the smoothed track and the root is a real feature name in
  this database, not an error.

Concurrency: this module exposes :func:`process_seq_chunk` as the
per-worker entry point for ``multiprocessing.Pool``. Workers are
initialised once via :func:`worker_initializer` with the full set of
``(HierarchyIndex, FeaturesForWorker)`` pairs -- one per requested
feature set -- so a single pool can produce all output BEDs in one
streaming pass over the combined BED rather than re-spawning per
feature set. The chunk-flushing boundary is the sequence id (the
BED's first column): a chunk always contains complete sequences,
never a fragment, because smoothing needs the full flanking context
for each sequence.

In *assembly mode* (the default; ``preserve_input_order=True`` in
:func:`karyoscope.core.annotate.annotate`) workers also receive a
per-feature-set temp directory and write their smoothed /
presmoothed output for each ``(feature_set, sequence)`` pair
directly to disk -- the worker returns only small line-count
metadata, avoiding the multi-GB IPC transfer that would otherwise
occur for whole-chromosome chunks. In *reads mode*
(``preserve_input_order=False``, intended for FASTQ/long-read input
with millions of small sequences) workers return the output lines
in a dict and the main process writes them straight to a single
per-feature-set BED.
"""

from __future__ import annotations

import gzip
import logging
import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO

from karyoscope.core.io.features import NOVEL_NAME, Features
from karyoscope.core.io.hierarchy import REQUIRED_ROOT, Hierarchy
from karyoscope.exceptions import KaryoscopeError

logger = logging.getLogger(__name__)


class SmoothError(KaryoscopeError):
    """Problems running the smoothing pass."""


#: Maximum distance (in BED coordinates) between two intervals for them
#: to be considered "flanking" for smoothing purposes. Matches the
#: archive's hard-coded constant; intervals separated by more than this
#: are treated as independent regions.
DEFAULT_MAX_GAP = 1000

#: FIX(asat): default density floor for :func:`density_filter`. 0.0 disables it
#: (original behaviour). Sparse-probe DBs (e.g. the divergent alpha-satellite
#: set) set this > 0 so a same-feature run that is mostly bridged-novel gets
#: reverted to ``novel`` -- kills the sparse-paralogous-chain inflation that the
#: root-as-gap guard alone does not (adjacent cross-mapped hits < max_gap apart
#: still chain). A genuine dense array (novel fraction near 0) is untouched.
DEFAULT_MIN_DENSITY = 0.0


# --- HierarchyIndex --------------------------------------------------


@dataclass
class HierarchyIndex:
    """Ancestor / LCA queries for one feature set's hierarchy.

    Wraps the parent map from :class:`karyoscope.core.io.hierarchy.Hierarchy`
    with memoised :meth:`get_ancestors` and :meth:`get_lca` queries. One
    instance per feature set; cheap to construct (microseconds for the
    real CHM13 database).

    Attributes
    ----------
    feature_set
        Name of the feature set this index is for.
    parent_of
        ``{child: parent}`` map. The root is absent.
    root
        Name of this feature set's root node (always
        :data:`REQUIRED_ROOT` for v0.1).
    """

    feature_set: str
    parent_of: dict[str, str]
    root: str
    _ancestor_cache: dict[str, list[str]] = field(default_factory=dict)
    _lca_cache: dict[tuple[str, str], str | None] = field(default_factory=dict)

    @classmethod
    def from_hierarchy(cls, hierarchy: Hierarchy, feature_set: str) -> HierarchyIndex:
        """Build an index for one feature set from a parsed hierarchy.

        Trusts that the caller already ran :func:`validate_hierarchy` —
        will raise :class:`SmoothError` on the obvious shape problems
        (missing root, non-categorized root) but doesn't re-run the
        full validation pass.
        """
        parent_of = hierarchy.parent_map(feature_set)
        all_nodes = set(parent_of.keys()) | set(parent_of.values())
        roots = all_nodes - set(parent_of.keys())
        if len(roots) != 1:
            raise SmoothError(
                f"feature set {feature_set!r} does not have exactly one root: {sorted(roots)!r}"
            )
        root = next(iter(roots))
        if root != REQUIRED_ROOT:
            raise SmoothError(
                f"feature set {feature_set!r}: root is {root!r}, must be {REQUIRED_ROOT!r}"
            )
        return cls(feature_set=feature_set, parent_of=parent_of, root=root)

    def get_ancestors(self, node: str) -> list[str]:
        """Return ``[node, parent, ..., root]`` for ``node``.

        For an unknown node (one that's not in the hierarchy at all),
        returns ``[node]`` — a single-element list. This matches the
        archive's behaviour and is the right default for the novel
        sentinel and for any unknown name a caller might pass.
        """
        cached = self._ancestor_cache.get(node)
        if cached is not None:
            return cached
        out = [node]
        current = node
        while current in self.parent_of:
            current = self.parent_of[current]
            out.append(current)
        self._ancestor_cache[node] = out
        return out

    def is_ancestor(self, node: str, potential_ancestor: str) -> bool:
        """Return True if ``potential_ancestor`` is on the path from ``node`` to root.

        A node is considered its own ancestor (matches the archive).
        """
        return potential_ancestor in self.get_ancestors(node)

    def get_lca(self, node1: str, node2: str) -> str | None:
        """Return the lowest common ancestor of ``node1`` and ``node2``.

        Returns ``None`` if the two nodes have no common ancestor
        (impossible for well-formed inputs but defensive).
        """
        key = (node1, node2) if node1 <= node2 else (node2, node1)
        if key in self._lca_cache:
            return self._lca_cache[key]
        ancestors1 = self.get_ancestors(node1)
        ancestors2 = set(self.get_ancestors(node2))
        for a in ancestors1:
            if a in ancestors2:
                self._lca_cache[key] = a
                return a
        self._lca_cache[key] = None
        return None


# --- Interval representation ----------------------------------------


@dataclass
class Interval:
    """A single BED-style interval being smoothed.

    Mutable on purpose: :func:`smooth_intervals` updates the
    ``feature`` field in place during the fixed-point loop.
    """

    seq_name: str
    start: int
    end: int
    feature: str
    is_novel: bool


# --- The algorithm itself -------------------------------------------


def smooth_intervals(
    intervals: list[Interval],
    index: HierarchyIndex,
    *,
    max_gap: int = DEFAULT_MAX_GAP,
    out_stats: dict[str, int] | None = None,
) -> list[Interval]:
    """Promote noisy intermediate intervals to their LCA with flankers.

    This is the faithful port of the archive's ``smooth_chunk``.
    Operates on a list of :class:`Interval` objects belonging to a
    single sequence; returns a *new* list with feature labels updated
    where smoothing promoted them. Repeats passes until no further
    changes happen (fixed point).

    Parameters
    ----------
    intervals
        Per-sequence intervals in left-to-right order.
    index
        The :class:`HierarchyIndex` for the feature set being smoothed.
    max_gap
        Maximum BED-coordinate gap between two intervals for them to
        be considered "flanking" for smoothing.
    out_stats
        Optional mutable dict the function populates with diagnostic
        counters. When provided, on return it carries:
        ``{"passes": int}`` -- the number of fixed-point iterations
        actually executed (always >= 1; converges when an iteration
        makes no changes). Useful for debugging perf, e.g. confirming
        the algorithm isn't pathologically iterating.
    """
    if not intervals:
        return []

    # Work on a copy so the caller's list isn't mutated.
    work = [
        Interval(
            seq_name=iv.seq_name,
            start=iv.start,
            end=iv.end,
            feature=iv.feature,
            is_novel=iv.is_novel,
        )
        for iv in intervals
    ]

    pass_count = 0
    while True:
        pass_count += 1
        changes_made = False
        i = 0
        was_related = [False] * len(work)

        while i < len(work):
            window_start = i
            related_indices = [i]
            first_unrelated_idx: int | None = None

            # Forward scan
            disallowed_feats: set[str] = set()
            anchor_feat = work[i].feature
            last_related_idx = i
            last_related_feat = anchor_feat
            last_related_end = work[i].end
            j = i
            while j + 1 < len(work):
                next_feat = work[j + 1].feature
                next_start = work[j + 1].start
                if next_start > last_related_end + max_gap:
                    break
                # FIX(asat): a novel run carries the ROOT label, and the root is an
                # ancestor of every feature, so is_ancestor() alone would treat a
                # novel interval as a 'related' continuation and advance
                # last_related_end across its WHOLE length -- bypassing max_gap
                # (which only measures gaps BETWEEN intervals, not the LENGTH of
                # one; annotate emits one contiguous interval per novel run). Keep
                # the intended short-gap smoothing (a novel run <= max_gap is still
                # bridged and promoted to the flankers' LCA) but stop a novel run
                # LONGER than max_gap from being bridged -- that is the 21 Mb-arm
                # artifact. Sparse chains (many <= max_gap novel gaps) are handled
                # separately by density_filter.
                if index.is_ancestor(last_related_feat, next_feat) and not (
                    work[j + 1].is_novel and work[j + 1].end - work[j + 1].start > max_gap
                ):
                    if next_feat in disallowed_feats:
                        break
                    related_indices.append(j + 1)
                    last_related_idx = j + 1
                    last_related_feat = next_feat
                    last_related_end = work[j + 1].end
                    j += 1
                elif not index.is_ancestor(next_feat, last_related_feat):
                    disallowed_ancestors = index.get_ancestors(next_feat)
                    disallowed_feats.update(disallowed_ancestors)
                    if first_unrelated_idx is None and not was_related[j + 1]:
                        first_unrelated_idx = j + 1
                    j += 1
                else:
                    break
            window_end = last_related_idx

            # Backward scan from the current peak
            k = window_end
            peak_feat = work[k].feature
            while k + 1 < len(work):
                next_feat = work[k + 1].feature
                next_start = work[k + 1].start
                if next_start > last_related_end + max_gap:
                    break
                # FIX(asat): mirror the forward scan's length-bounded novel guard
                # (a novel run longer than max_gap is not a bridgeable flanker).
                if index.is_ancestor(next_feat, last_related_feat) and not (
                    work[k + 1].is_novel and work[k + 1].end - work[k + 1].start > max_gap
                ):
                    related_indices.append(k + 1)
                    last_related_idx = k + 1
                    last_related_feat = next_feat
                    last_related_end = work[k + 1].end
                    k += 1
                elif not index.is_ancestor(last_related_feat, next_feat):
                    if peak_feat in index.get_ancestors(next_feat):
                        break
                    if first_unrelated_idx is None and not was_related[k + 1]:
                        first_unrelated_idx = k + 1
                    k += 1
                else:
                    break
            window_end = last_related_idx

            # Archive quirk: drop the last "related" index from the
            # was_related bookkeeping (preserved for fidelity).
            if len(related_indices) >= 2:
                related_indices.pop()
            for idx in set(related_indices):
                was_related[idx] = True

            # LCA-based promotion within the window
            if window_end > window_start:
                left_feat = work[window_start].feature
                right_feat = work[window_end].feature
                lca = index.get_lca(left_feat, right_feat)
                if lca is not None:
                    for w_idx in range(window_start + 1, window_end):
                        original = work[w_idx].feature
                        if index.is_ancestor(lca, original) and original != lca:
                            work[w_idx].feature = lca
                            changes_made = True

            # Advance i (faithful to the archive)
            if first_unrelated_idx is not None:
                if first_unrelated_idx > window_end > window_start:
                    i = window_end
                else:
                    i = first_unrelated_idx
            elif window_end > i:
                i = window_end
            else:
                while i < len(work) and was_related[i]:
                    i += 1

        if not changes_made:
            break

    if out_stats is not None:
        out_stats["passes"] = pass_count
    return work


def density_filter(
    intervals: list[Interval],
    root: str,
    *,
    min_density: float = DEFAULT_MIN_DENSITY,
) -> list[Interval]:
    """Revert bridged-novel bp in low-density same-feature runs to ``novel``.

    FIX(asat): the smoother bridges a same-leaf hit every < ``max_gap`` amid
    novel into one contiguous run even when the real (non-novel) support is a
    tiny fraction of the span (a paralogous k-mer sprinkle, e.g. chr3's D3Z1
    monomer cross-mapping across the chr10 centromere -> a 1.8 Mb ``D3Z1``
    block from 0.18 Mb of real hits). The root-as-gap guard does not catch this
    because each individual gap is < ``max_gap``.

    For each maximal run of consecutive same-``feature`` intervals whose feature
    is not the root, compute ``real_bp / span`` where ``real_bp`` is the summed
    length of the non-novel intervals. If that is below ``min_density``, the run
    is a bridging artifact: its *promoted-novel* intervals (``is_novel``) are
    reverted to ``(root, novel)`` while the genuine hits keep their label. A
    dense real array (density ~1) is left untouched. ``min_density == 0``
    disables the pass (returns the input unchanged).

    Operates on the pre-:func:`merge_adjacent` list so per-interval ``is_novel``
    is still available; returns a new list.
    """
    if min_density <= 0.0 or not intervals:
        return intervals
    out = [
        Interval(iv.seq_name, iv.start, iv.end, iv.feature, iv.is_novel)
        for iv in intervals
    ]
    i = 0
    n = len(out)
    while i < n:
        j = i
        while j + 1 < n and out[j + 1].feature == out[i].feature:
            j += 1
        run = out[i : j + 1]
        feat = run[0].feature
        if feat != root:
            span = sum(iv.end - iv.start for iv in run)
            real_bp = sum(iv.end - iv.start for iv in run if not iv.is_novel)
            if span > 0 and real_bp / span < min_density:
                for iv in run:
                    if iv.is_novel:
                        iv.feature = root
        i = j + 1
    return out


def merge_adjacent(
    intervals: list[Interval],
    root: str,
) -> list[Interval]:
    """Coalesce contiguous intervals sharing the same feature label.

    Matches the archive's ``merge_intervals``: adjacent intervals with
    identical features merge if they're contiguous (``next.start ==
    current.end``). Intervals at the root level additionally require
    matching ``is_novel`` — this keeps the novel boundary visible in
    the output even after smoothing has projected several novel-and-
    nonnovel runs up to the root.
    """
    if not intervals:
        return []

    merged: list[Interval] = []
    current = Interval(
        seq_name=intervals[0].seq_name,
        start=intervals[0].start,
        end=intervals[0].end,
        feature=intervals[0].feature,
        is_novel=intervals[0].is_novel,
    )

    for nxt in intervals[1:]:
        same_seq = nxt.seq_name == current.seq_name
        same_feature = nxt.feature == current.feature
        contiguous = nxt.start == current.end
        compatible_novelty = current.is_novel == nxt.is_novel if current.feature == root else True
        if same_seq and same_feature and contiguous and compatible_novelty:
            current.end = nxt.end
        else:
            merged.append(current)
            current = Interval(
                seq_name=nxt.seq_name,
                start=nxt.start,
                end=nxt.end,
                feature=nxt.feature,
                is_novel=nxt.is_novel,
            )

    merged.append(current)
    return merged


# --- BED I/O for the worker ------------------------------------------


def _render_for_output(iv: Interval, root: str) -> str:
    """Render one interval as a BED line.

    Maps ``feature == root`` with ``is_novel`` back to ``"novel"`` so
    the internal sentinel never leaks to users.
    """
    label = NOVEL_NAME if iv.feature == root and iv.is_novel else iv.feature
    return f"{iv.seq_name}\t{iv.start}\t{iv.end}\t{label}\n"


# --- Worker entry point for multiprocessing --------------------------

# These globals are populated by worker_initializer on each Pool worker
# process. Workers can't pickle the HierarchyIndex caches through every
# imap call cheaply, so we share constants once per worker. All
# requested feature sets live in the worker simultaneously -- the
# single-pool architecture (see module docstring) processes each chunk
# for every feature set in one invocation.
_worker_indices: dict[str, HierarchyIndex] | None = None
_worker_features_by_fs: dict[str, FeaturesForWorker] | None = None
_worker_feature_sets: list[str] | None = None
# Set in assembly mode (per-(fs, seq) temp files); ``None`` in reads
# mode (worker returns output lines via IPC).
_worker_tmpdir_by_fs: dict[str, Path] | None = None
# FIX(asat): density-filter threshold (see density_filter). Set per-run by
# worker_initializer from the DB manifest; 0.0 disables (original behaviour).
_worker_min_density: float = DEFAULT_MIN_DENSITY


# A trimmed worker-friendly view of the features table. We hold it as a
# plain dict to keep the pickle payload tiny when forking workers; the
# full Features object is fine but the only data the worker needs is
# the id->name map for one feature set.
@dataclass
class FeaturesForWorker:
    """Per-feature-set id-to-name mapping for the smoothing worker.

    Built once on the main process from the full
    :class:`karyoscope.core.io.features.Features` and sent to the pool
    as part of ``initargs``.
    """

    feature_set: str
    id_to_name: dict[int, str]
    # The name used for ids absent from the table. For 5b we made
    # missing ids a hard error in the main process; the worker mirrors
    # that with SmoothError if it ever sees one (shouldn't happen
    # because the C++ binary only emits ids in the index).
    novel_label: str


def make_features_for_worker(features: Features, feature_set: str) -> FeaturesForWorker:
    """Project a full :class:`Features` onto one feature set for a worker."""
    id_to_name = {fid: row[feature_set] for fid, row in features.table.items()}
    return FeaturesForWorker(
        feature_set=feature_set,
        id_to_name=id_to_name,
        novel_label=NOVEL_NAME,
    )


_SAFE_FILENAME_RE = re.compile(r"[^A-Za-z0-9_.-]")


def safe_filename(seq_name: str) -> str:
    """Sanitise a FASTA sequence name for use as a path component.

    Replaces any character outside ``[A-Za-z0-9_.-]`` with ``_``.
    Real FASTA sequence names from human assemblies are already
    filesystem-safe (``chr1_MATERNAL`` etc.); the regex defends
    against the occasional ``/``, space, or pipe character that
    would otherwise break the per-(feature_set, sequence) temp file
    paths the assembly-mode worker writes to.

    Two distinct sequence names sanitising to the same filename
    would silently concatenate their output. Real genomes don't hit
    this; the main process raises if it detects a collision while
    enumerating input sequences.
    """
    return _SAFE_FILENAME_RE.sub("_", seq_name)


def worker_initializer(
    indices: dict[str, HierarchyIndex],
    features_by_fs: dict[str, FeaturesForWorker],
    feature_sets: list[str],
    tmpdir_by_fs: dict[str, Path] | None,
    log_level: int = logging.WARNING,
) -> None:
    """Initialise a multiprocessing-pool worker for the single-pool path.

    Workers process every requested feature set per chunk, so all of
    the per-feature-set state (``HierarchyIndex`` + ``FeaturesForWorker``)
    arrives once via the ``initargs`` rather than per chunk. The
    ``feature_sets`` argument fixes the iteration order so the worker
    visits feature sets in a deterministic order across chunks.

    ``tmpdir_by_fs`` selects the output mode:

    * **dict**: assembly mode. The worker writes each
      ``(feature_set, sequence)`` pair's smoothed / presmoothed BED
      lines directly to ``{tmpdir_by_fs[fs]}/{safe_seq_name}.smo`` /
      ``.pre``, and returns only ``{fs: {seq: (n_pre, n_smo)}}``
      line counts. Avoids the multi-GB IPC transfer that whole-chr
      chunks would otherwise produce.
    * **None**: reads mode. The worker returns
      ``{fs: {seq: (smo_lines, pre_lines)}}`` for the main process
      to write. Reads produce small per-chunk outputs so IPC is
      cheap; the absence of input-order tracking is acceptable
      because read order isn't meaningful in FASTQ/long-read input.

    Also installs a stderr log handler in the worker. With the
    ``spawn`` context (which :mod:`karyoscope.core.annotate` uses),
    the worker is a fresh Python process with no logging
    configuration; without this step ``logger.info`` calls vanish.
    """
    import sys

    handler = logging.StreamHandler(stream=sys.stderr)
    handler.setFormatter(
        logging.Formatter(
            "%(asctime)s %(levelname)s [worker %(process)d] %(message)s",
            datefmt="%H:%M:%S",
        )
    )
    root = logging.getLogger()
    root.handlers = [handler]
    root.setLevel(log_level)

    global _worker_indices, _worker_features_by_fs
    global _worker_feature_sets, _worker_tmpdir_by_fs, _worker_min_density
    _worker_indices = indices
    _worker_features_by_fs = features_by_fs
    _worker_feature_sets = list(feature_sets)
    # FIX(asat): density-filter threshold via env (inherited by spawn workers);
    # annotate() sets it from the DB manifest's smoothing.min_density. Default
    # 0.0 disables the pass -> original behaviour for every other DB.
    import os as _os
    try:
        _worker_min_density = float(_os.environ.get("KARYOSCOPE_MIN_DENSITY", DEFAULT_MIN_DENSITY))
    except ValueError:
        _worker_min_density = DEFAULT_MIN_DENSITY
    _worker_tmpdir_by_fs = tmpdir_by_fs


def _smooth_one_seq_for_fs(
    seq_name: str,
    intervals_raw: list[tuple[int, int, int]],
    feature_set: str,
    index: HierarchyIndex,
    features: FeaturesForWorker,
) -> tuple[list[str], list[str]]:
    """Smooth one sequence for one feature set; return (smo, pre) lines.

    Translates raw ``(start, end, feature_id)`` tuples to FS-specific
    :class:`Interval` objects, runs the merge / smooth / merge pipeline,
    renders to BED lines, and logs per-sequence progress at ``INFO``.

    Pulled out of :func:`process_seq_chunk` so the worker can iterate
    over feature sets while letting Python free each FS's intermediate
    state (Interval list + smoothing scratch) between iterations.
    """
    import time as _time

    t0 = _time.perf_counter()
    root = index.root
    id_to_name = features.id_to_name

    intervals: list[Interval] = []
    for start, end, fid in intervals_raw:
        if fid == 0:
            internal_label = root
            is_novel = True
        else:
            name = id_to_name.get(fid)
            if name is None:
                raise SmoothError(
                    f"feature id {fid} not in features.tsv for feature set {feature_set!r}"
                )
            internal_label = name
            is_novel = False
        intervals.append(
            Interval(
                seq_name=seq_name,
                start=start,
                end=end,
                feature=internal_label,
                is_novel=is_novel,
            )
        )

    n_in = len(intervals)
    merged_pre = merge_adjacent(intervals, root)
    pre_lines = [_render_for_output(iv, root) for iv in merged_pre]
    stats: dict[str, int] = {}
    smoothed = smooth_intervals(intervals, index, out_stats=stats)
    smoothed = density_filter(smoothed, root, min_density=_worker_min_density)
    merged_post = merge_adjacent(smoothed, root)
    smo_lines = [_render_for_output(iv, root) for iv in merged_post]

    dt = _time.perf_counter() - t0
    passes = stats.get("passes", 0)
    logger.info(
        "smoothed %s/%s: %d intervals -> %d presmoothed, %d smoothed in %.1fs (%d pass%s)",
        feature_set,
        seq_name,
        n_in,
        len(pre_lines),
        len(smo_lines),
        dt,
        passes,
        "" if passes == 1 else "es",
    )
    return smo_lines, pre_lines


def process_seq_chunk(
    chunk: list[str],
) -> dict[str, dict[str, tuple[int, int] | tuple[list[str], list[str]]]]:
    """Smooth one chunk of combined-BED lines for *every* feature set.

    The chunk must consist of complete sequences (no fragments).
    Input lines are 4-column combined BED records with the integer
    feature id in column 4 (the C++ binary's output format).

    The chunk is parsed *once* into per-sequence raw
    ``(start, end, feature_id)`` tuples; then for each sequence the
    function iterates over every requested feature set, building an
    FS-specific :class:`Interval` list and running the smoothing
    pipeline. Letting each FS's intermediate state go out of scope
    between iterations keeps peak worker RAM bounded by one FS's
    working set rather than ``N_fs x`` that.

    Return type depends on the mode set in :func:`worker_initializer`:

    * **Assembly mode** (``_worker_tmpdir_by_fs`` is a dict): for each
      ``(fs, seq)`` the BED lines are written directly to
      ``{tmpdir}/{safe_seq_name}.smo`` / ``.pre``. The return value
      is ``{fs: {seq: (n_pre_lines, n_smo_lines)}}`` -- counts only,
      so IPC stays microscopic regardless of chr1-sized output.
    * **Reads mode** (``_worker_tmpdir_by_fs`` is ``None``): returns
      ``{fs: {seq: (smo_lines, pre_lines)}}``. Reads produce small
      per-chunk outputs so this is cheap to transfer.

    Novel feature ids (id 0) are mapped to the hierarchy root
    sentinel internally during smoothing; :func:`_render_for_output`
    rewrites them back to :data:`karyoscope.core.io.features.NOVEL_NAME`
    so users never see the internal sentinel.
    """
    if _worker_indices is None or _worker_features_by_fs is None or _worker_feature_sets is None:
        raise SmoothError("worker_initializer was not called; cannot process chunk")
    indices = _worker_indices
    features_by_fs = _worker_features_by_fs
    feature_sets = _worker_feature_sets
    tmpdir_by_fs = _worker_tmpdir_by_fs

    # Phase 1: parse the chunk into per-sequence raw intervals.
    # Holding tuples (24 bytes each + python overhead) instead of
    # full Interval objects (~150 bytes each) keeps the parsed-but-
    # not-yet-smoothed state small. The per-FS Interval list is
    # built lazily inside the smoothing loop.
    by_seq: dict[str, list[tuple[int, int, int]]] = {}
    seq_order: list[str] = []
    for raw in chunk:
        line = raw.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            raise SmoothError(f"malformed line in chunk (need at least 4 columns): {raw!r}")
        try:
            start = int(parts[1])
            end = int(parts[2])
            feature_id = int(parts[3])
        except ValueError as e:
            raise SmoothError(f"non-integer field in chunk line: {raw!r}") from e
        seq_name = parts[0]
        bucket = by_seq.get(seq_name)
        if bucket is None:
            bucket = []
            by_seq[seq_name] = bucket
            seq_order.append(seq_name)
        bucket.append((start, end, feature_id))

    # Phase 2: per-sequence, per-feature-set smoothing.
    out: dict[str, dict[str, tuple[int, int] | tuple[list[str], list[str]]]] = {
        fs: {} for fs in feature_sets
    }
    for seq in seq_order:
        intervals_raw = by_seq.pop(seq)
        for fs in feature_sets:
            smo_lines, pre_lines = _smooth_one_seq_for_fs(
                seq_name=seq,
                intervals_raw=intervals_raw,
                feature_set=fs,
                index=indices[fs],
                features=features_by_fs[fs],
            )
            if tmpdir_by_fs is not None:
                # Assembly mode: write directly, return counts only.
                fs_dir = tmpdir_by_fs[fs]
                safe = safe_filename(seq)
                if smo_lines:
                    with (fs_dir / f"{safe}.smo").open("a") as f:
                        f.writelines(smo_lines)
                if pre_lines:
                    with (fs_dir / f"{safe}.pre").open("a") as f:
                        f.writelines(pre_lines)
                out[fs][seq] = (len(pre_lines), len(smo_lines))
            else:
                out[fs][seq] = (smo_lines, pre_lines)
    return out


# --- Chunked reading respecting sequence boundaries ------------------


def chunked_seq_reader(
    path: Path,
    chunk_size: int = 50000,
) -> Iterator[list[str]]:
    """Yield chunks of BED lines, flushing only at sequence boundaries.

    A "chunk" is a list of newline-terminated strings. The reader
    holds back the current sequence's lines until either the buffer
    has at least ``chunk_size`` lines AND the next line is from a new
    sequence, or EOF is reached. This guarantees that ``smooth_intervals``
    never sees a partial sequence.

    Supports plain text and gzip-compressed inputs based on extension.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:  # type: ignore[operator]
        yield from _chunked_seq_reader_from_handle(fh, chunk_size)


def _chunked_seq_reader_from_handle(
    fh: IO[str],
    chunk_size: int,
) -> Iterator[list[str]]:
    """Implementation factored out for testability with in-memory streams."""
    chunk: list[str] = []
    last_seq_id: str | None = None

    for line in fh:
        seq_id, _, _ = line.partition("\t")
        if seq_id != last_seq_id and last_seq_id is not None and len(chunk) >= chunk_size:
            yield chunk
            chunk = []
        chunk.append(line)
        last_seq_id = seq_id

    if chunk:
        yield chunk
