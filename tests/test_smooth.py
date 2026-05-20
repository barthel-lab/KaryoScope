"""Unit tests for :mod:`karyoscope.core.smooth`."""

from __future__ import annotations

import io

import pytest

from karyoscope.core.io.hierarchy import REQUIRED_ROOT, Hierarchy, HierarchyRow
from karyoscope.core.smooth import (
    DEFAULT_MAX_GAP,
    HierarchyIndex,
    Interval,
    SmoothError,
    _chunked_seq_reader_from_handle,
    merge_adjacent,
    smooth_intervals,
)

# --- Test fixture: a small but non-trivial hierarchy ----------------


@pytest.fixture
def region_hierarchy() -> Hierarchy:
    """Three-level region hierarchy matching the dummy db.

    Tree:
                categorized
                     |
                centromeric
                /         \\
             aSat         HSat
            /    \\         |
           rA     rB        rC
    """
    return Hierarchy(
        rows=[
            HierarchyRow("region", "centromeric", "categorized"),
            HierarchyRow("region", "aSat", "centromeric"),
            HierarchyRow("region", "HSat", "centromeric"),
            HierarchyRow("region", "rA", "aSat"),
            HierarchyRow("region", "rB", "aSat"),
            HierarchyRow("region", "rC", "HSat"),
        ]
    )


@pytest.fixture
def region_index(region_hierarchy: Hierarchy) -> HierarchyIndex:
    return HierarchyIndex.from_hierarchy(region_hierarchy, "region")


# --- HierarchyIndex -------------------------------------------------


def test_index_ancestors_root_to_leaf(region_index: HierarchyIndex) -> None:
    assert region_index.get_ancestors("rA") == ["rA", "aSat", "centromeric", "categorized"]


def test_index_ancestors_of_root(region_index: HierarchyIndex) -> None:
    assert region_index.get_ancestors("categorized") == ["categorized"]


def test_index_ancestors_of_unknown_node(region_index: HierarchyIndex) -> None:
    """Unknown node returns a single-element list (archive behaviour)."""
    assert region_index.get_ancestors("mystery_feature") == ["mystery_feature"]


def test_index_is_ancestor_reflexive(region_index: HierarchyIndex) -> None:
    """A node is its own ancestor (matches the archive)."""
    assert region_index.is_ancestor("rA", "rA")


def test_index_is_ancestor_true_cases(region_index: HierarchyIndex) -> None:
    assert region_index.is_ancestor("rA", "aSat")
    assert region_index.is_ancestor("rA", "centromeric")
    assert region_index.is_ancestor("rA", "categorized")


def test_index_is_ancestor_false_cases(region_index: HierarchyIndex) -> None:
    assert not region_index.is_ancestor("rA", "rB")  # siblings
    assert not region_index.is_ancestor("rA", "HSat")  # different subtree
    assert not region_index.is_ancestor("aSat", "rA")  # parent of, not ancestor of


def test_index_lca_siblings(region_index: HierarchyIndex) -> None:
    assert region_index.get_lca("rA", "rB") == "aSat"


def test_index_lca_across_subtrees(region_index: HierarchyIndex) -> None:
    assert region_index.get_lca("rA", "rC") == "centromeric"
    assert region_index.get_lca("rA", "HSat") == "centromeric"


def test_index_lca_self(region_index: HierarchyIndex) -> None:
    assert region_index.get_lca("rA", "rA") == "rA"


def test_index_lca_includes_root(region_index: HierarchyIndex) -> None:
    """Two top-level siblings have the root as their LCA."""
    # 'centromeric' is the only top-level child of categorized in this
    # fixture, so use an unknown node + centromeric instead.
    # Actually: unknown node returns just itself as ancestors → no LCA.
    # Set up a fresh hierarchy with two top-level siblings:
    h = Hierarchy(
        rows=[
            HierarchyRow("x", "A", "categorized"),
            HierarchyRow("x", "B", "categorized"),
        ]
    )
    idx = HierarchyIndex.from_hierarchy(h, "x")
    assert idx.get_lca("A", "B") == "categorized"


def test_index_lca_none_for_disjoint_nodes(region_index: HierarchyIndex) -> None:
    """Two unknown nodes have no common ancestor."""
    assert region_index.get_lca("unknown_a", "unknown_b") is None


def test_index_caches_results(region_index: HierarchyIndex) -> None:
    """get_ancestors / get_lca memoise — second call hits the cache."""
    region_index.get_ancestors("rA")
    region_index.get_lca("rA", "rB")
    # Caches are populated
    assert "rA" in region_index._ancestor_cache
    assert any("rA" in k for k in region_index._lca_cache)


def test_index_rejects_multiple_roots() -> None:
    h = Hierarchy(
        rows=[
            HierarchyRow("x", "a", "rootA"),
            HierarchyRow("x", "b", "rootB"),
        ]
    )
    with pytest.raises(SmoothError, match="exactly one root"):
        HierarchyIndex.from_hierarchy(h, "x")


def test_index_rejects_non_categorized_root() -> None:
    h = Hierarchy(
        rows=[
            HierarchyRow("x", "a", "wrongroot"),
        ]
    )
    with pytest.raises(SmoothError, match="must be 'categorized'"):
        HierarchyIndex.from_hierarchy(h, "x")


# --- merge_adjacent ----------------------------------------------------


def test_merge_adjacent_collapses_same_feature() -> None:
    intervals = [
        Interval("seqA", 0, 100, "rA", False),
        Interval("seqA", 100, 200, "rA", False),
        Interval("seqA", 200, 300, "rA", False),
    ]
    out = merge_adjacent(intervals, REQUIRED_ROOT)
    assert len(out) == 1
    assert (out[0].start, out[0].end) == (0, 300)


def test_merge_adjacent_keeps_distinct_features() -> None:
    intervals = [
        Interval("seqA", 0, 100, "rA", False),
        Interval("seqA", 100, 200, "rB", False),
    ]
    out = merge_adjacent(intervals, REQUIRED_ROOT)
    assert len(out) == 2


def test_merge_adjacent_doesnt_bridge_gaps() -> None:
    """Non-contiguous same-feature intervals stay separate."""
    intervals = [
        Interval("seqA", 0, 100, "rA", False),
        Interval("seqA", 150, 200, "rA", False),  # gap from 100 to 150
    ]
    out = merge_adjacent(intervals, REQUIRED_ROOT)
    assert len(out) == 2


def test_merge_adjacent_root_keeps_novel_boundary() -> None:
    """At root level, is_novel boundaries are preserved."""
    intervals = [
        Interval("seqA", 0, 100, REQUIRED_ROOT, True),  # novel
        Interval("seqA", 100, 200, REQUIRED_ROOT, False),  # not novel, same feature
    ]
    out = merge_adjacent(intervals, REQUIRED_ROOT)
    # Stays separate because is_novel differs at root
    assert len(out) == 2


def test_merge_adjacent_nonroot_ignores_novelty() -> None:
    """Non-root features merge regardless of is_novel."""
    intervals = [
        Interval("seqA", 0, 100, "rA", True),
        Interval("seqA", 100, 200, "rA", False),
    ]
    out = merge_adjacent(intervals, REQUIRED_ROOT)
    assert len(out) == 1


def test_merge_adjacent_empty() -> None:
    assert merge_adjacent([], REQUIRED_ROOT) == []


# --- smooth_intervals --------------------------------------------------


def test_smooth_promotes_novel_between_siblings(region_index: HierarchyIndex) -> None:
    """rA → novel → rB should promote the novel to aSat (LCA of rA, rB)."""
    intervals = [
        Interval("seqA", 0, 100, "rA", False),
        Interval("seqA", 100, 200, REQUIRED_ROOT, True),  # novel
        Interval("seqA", 200, 300, "rB", False),
    ]
    out = smooth_intervals(intervals, region_index)
    # The middle interval should be promoted to aSat
    assert out[1].feature == "aSat"
    # Flanking intervals unchanged
    assert out[0].feature == "rA"
    assert out[2].feature == "rB"


def test_smooth_promotes_across_deeper_lca(region_index: HierarchyIndex) -> None:
    """rA → novel → rC: LCA is centromeric, two levels up."""
    intervals = [
        Interval("seqA", 0, 100, "rA", False),
        Interval("seqA", 100, 200, REQUIRED_ROOT, True),
        Interval("seqA", 200, 300, "rC", False),
    ]
    out = smooth_intervals(intervals, region_index)
    assert out[1].feature == "centromeric"


def test_smooth_no_op_when_no_flankers(region_index: HierarchyIndex) -> None:
    """A single interval with no neighbours stays put."""
    intervals = [Interval("seqA", 0, 100, "rA", False)]
    out = smooth_intervals(intervals, region_index)
    assert out[0].feature == "rA"


def test_smooth_respects_max_gap(region_index: HierarchyIndex) -> None:
    """An interval too far from the flanker is not smoothed."""
    intervals = [
        Interval("seqA", 0, 100, "rA", False),
        # Gap from 100 to 100 + DEFAULT_MAX_GAP + 1 = beyond max_gap
        Interval(
            "seqA", 100 + DEFAULT_MAX_GAP + 1, 100 + DEFAULT_MAX_GAP + 100, REQUIRED_ROOT, True
        ),
        Interval("seqA", 100 + DEFAULT_MAX_GAP + 100, 100 + DEFAULT_MAX_GAP + 200, "rB", False),
    ]
    out = smooth_intervals(intervals, region_index)
    # The novel interval is too far from rA to be considered part of the
    # rA window; its label should be unchanged.
    assert out[1].feature == REQUIRED_ROOT


def test_smooth_empty(region_index: HierarchyIndex) -> None:
    assert smooth_intervals([], region_index) == []


def test_smooth_doesnt_demote(region_index: HierarchyIndex) -> None:
    """A specific feature shouldn't be replaced by a less-specific one.

    Algorithm only promotes ancestors of the LCA; an interval already
    labelled with the LCA or below shouldn't move.
    """
    intervals = [
        Interval("seqA", 0, 100, "rA", False),
        Interval("seqA", 100, 200, "rA", False),  # already specific
        Interval("seqA", 200, 300, "rB", False),
    ]
    out = smooth_intervals(intervals, region_index)
    # Middle interval was already "rA", LCA is "aSat". rA is more
    # specific than aSat, so it should NOT be promoted.
    assert out[1].feature == "rA"


# --- pass-count diagnostic (out_stats) -------------------------------


def test_smooth_out_stats_records_pass_count(region_index: HierarchyIndex) -> None:
    """``out_stats['passes']`` reports the number of fixed-point iterations.

    For a productive smoothing run the count is >=2: one pass that
    promotes the middle interval, plus a final no-op pass that
    confirms convergence and breaks the loop. Used by the smoothing
    worker to report per-sequence diagnostics in INFO logs.
    """
    intervals = [
        Interval("seqA", 0, 100, "rA", False),
        Interval("seqA", 100, 200, "centromeric", False),
        Interval("seqA", 200, 300, "rB", False),
    ]
    stats: dict[str, int] = {}
    smooth_intervals(intervals, region_index, out_stats=stats)
    assert stats["passes"] >= 2  # at least one productive + one confirmatory


def test_smooth_out_stats_no_op_is_one_pass(region_index: HierarchyIndex) -> None:
    """An input with nothing to smooth still records ``passes >= 1``.

    The fixed-point loop always runs the first pass; the loop exits
    when that pass made no changes. So an input that needs no
    smoothing reports exactly 1 pass.
    """
    intervals = [
        Interval("seqA", 0, 100, "rA", False),
        Interval("seqA", 100, 200, "rA", False),
    ]
    stats: dict[str, int] = {}
    smooth_intervals(intervals, region_index, out_stats=stats)
    assert stats["passes"] == 1


def test_smooth_out_stats_optional(region_index: HierarchyIndex) -> None:
    """``out_stats`` is optional -- callers that don't care don't pass it."""
    intervals = [Interval("seqA", 0, 100, "rA", False)]
    # No stats dict; must still work.
    result = smooth_intervals(intervals, region_index)
    assert len(result) == 1


# --- chunked_seq_reader -----------------------------------------------


def test_chunked_reader_yields_at_sequence_boundary() -> None:
    """A chunk should never contain a partial sequence."""
    # Three sequences, lots of records each
    lines: list[str] = []
    for seq_name in ("seqA", "seqB", "seqC"):
        for i in range(100):
            lines.append(f"{seq_name}\t{i * 10}\t{(i + 1) * 10}\t1\n")
    fh = io.StringIO("".join(lines))

    chunks = list(_chunked_seq_reader_from_handle(fh, chunk_size=50))
    # Each chunk must contain complete sequences only
    for chunk in chunks:
        seqs_in_chunk = {line.split("\t", 1)[0] for line in chunk}
        # Every seq present in a chunk should have all 100 of its lines
        for seq in seqs_in_chunk:
            count = sum(1 for line in chunk if line.startswith(seq + "\t"))
            assert count == 100, (
                f"sequence {seq!r} was split across chunks ({count} of 100 lines in one chunk)"
            )


def test_chunked_reader_small_input_one_chunk() -> None:
    """Input smaller than chunk_size produces a single chunk."""
    fh = io.StringIO("seqA\t0\t10\t1\nseqB\t0\t10\t2\n")
    chunks = list(_chunked_seq_reader_from_handle(fh, chunk_size=1000))
    assert len(chunks) == 1
    assert len(chunks[0]) == 2


def test_chunked_reader_empty_input() -> None:
    fh = io.StringIO("")
    chunks = list(_chunked_seq_reader_from_handle(fh, chunk_size=100))
    assert chunks == []
