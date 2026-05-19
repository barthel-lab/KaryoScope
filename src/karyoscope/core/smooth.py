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
initialised with the shared :class:`HierarchyIndex` and id-to-name
table via :func:`worker_initializer`. The chunk-flushing boundary is
the sequence id (the BED's first column): a chunk always contains
complete sequences, never a fragment, because smoothing needs the full
flanking context for each sequence.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO

from karyoscope.core.io.features import NOVEL_NAME, Features
from karyoscope.core.io.hierarchy import REQUIRED_ROOT, Hierarchy
from karyoscope.exceptions import KaryoscopeError


class SmoothError(KaryoscopeError):
    """Problems running the smoothing pass."""


#: Maximum distance (in BED coordinates) between two intervals for them
#: to be considered "flanking" for smoothing purposes. Matches the
#: archive's hard-coded constant; intervals separated by more than this
#: are treated as independent regions.
DEFAULT_MAX_GAP = 1000


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

    while True:
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
                if index.is_ancestor(last_related_feat, next_feat):
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
                if index.is_ancestor(next_feat, last_related_feat):
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

    return work


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
# process. Workers can't pickle the HierarchyIndex's caches efficiently
# through every imap call, so we share the constants once per worker.
_worker_index: HierarchyIndex | None = None
_worker_features: FeaturesForWorker | None = None
_worker_feature_set: str | None = None


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


def worker_initializer(
    index: HierarchyIndex,
    features_for_worker: FeaturesForWorker,
) -> None:
    """Initialise a multiprocessing-pool worker.

    Sets the per-worker globals used by :func:`process_seq_chunk`.
    Pickling cost is paid once per worker, not once per chunk.
    """
    global _worker_index, _worker_features, _worker_feature_set
    _worker_index = index
    _worker_features = features_for_worker
    _worker_feature_set = features_for_worker.feature_set


def process_seq_chunk(
    chunk: list[str],
) -> tuple[list[str], list[str]]:
    """Smooth one chunk of combined-BED lines for one feature set.

    The chunk must consist of complete sequences (no fragments).
    Input lines are 4-column combined BED records with the integer
    feature id in column 4 (the C++ binary's output format).

    Returns ``(smoothed_lines, presmoothed_lines)`` — two lists of
    newline-terminated BED records ready for the caller to write to
    files. The caller chooses which (or both) to write based on the
    user's ``--smooth`` and ``--keep-presmoothed`` flags.

    Each line is processed once: id translation happens here, the
    presmoothed track is the result of merging adjacent same-name
    runs, and the smoothed track adds the hierarchy-aware promotion
    pass on top. Running both passes in one worker avoids a re-read
    of the (potentially large) presmoothed BED.

    The archive smoothing treats novel intervals (id 0) as the root
    sentinel during the algorithm (``id_map['0'] = root``) and renames
    ``root + is_novel`` back to ``"novel"`` at output. We do the same:
    id 0 becomes the root internally and is unprojected to
    :data:`karyoscope.core.io.features.NOVEL_NAME` in
    :func:`_render_for_output`.
    """
    if _worker_index is None or _worker_features is None:
        raise SmoothError("worker_initializer was not called; cannot process chunk")
    index = _worker_index
    features = _worker_features
    root = index.root

    smoothed_out: list[str] = []
    presmoothed_out: list[str] = []

    current_seq: str | None = None
    current_intervals: list[Interval] = []

    def flush() -> None:
        if not current_intervals:
            return
        # Presmoothed: just merge adjacent same-feature runs.
        merged_pre = merge_adjacent(current_intervals, root)
        for iv in merged_pre:
            presmoothed_out.append(_render_for_output(iv, root))
        # Smoothed: run the algorithm, then re-merge.
        smoothed = smooth_intervals(current_intervals, index)
        merged_post = merge_adjacent(smoothed, root)
        for iv in merged_post:
            smoothed_out.append(_render_for_output(iv, root))

    for raw in chunk:
        line = raw.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            raise SmoothError(f"malformed line in chunk (need at least 4 columns): {raw!r}")
        seq_name = parts[0]
        try:
            start = int(parts[1])
            end = int(parts[2])
            feature_id = int(parts[3])
        except ValueError as e:
            raise SmoothError(f"non-integer field in chunk line: {raw!r}") from e

        if seq_name != current_seq:
            flush()
            current_seq = seq_name
            current_intervals = []

        # Translate the feature id to a hierarchy node label. id 0 is
        # the novel sentinel; ids >0 should be in the table (the
        # main-process annotate() function verifies the features file
        # is consistent with the index before kicking off workers, but
        # we still raise rather than silently mislabel if an id
        # slipped through).
        if feature_id == 0:
            internal_label = root
            is_novel = True
        else:
            name = features.id_to_name.get(feature_id)
            if name is None:
                raise SmoothError(
                    f"feature id {feature_id} not in features.tsv for "
                    f"feature set {features.feature_set!r}"
                )
            internal_label = name
            is_novel = False

        current_intervals.append(
            Interval(
                seq_name=seq_name,
                start=start,
                end=end,
                feature=internal_label,
                is_novel=is_novel,
            )
        )

    flush()
    return smoothed_out, presmoothed_out


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
