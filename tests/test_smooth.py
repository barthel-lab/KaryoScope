"""Unit tests for :mod:`karyoscope.core.smooth`."""

from __future__ import annotations

import io

import pytest

from karyoscope.core.io.hierarchy import REQUIRED_ROOT, Hierarchy, HierarchyRow
from karyoscope.core.smooth import (
    DEFAULT_MAX_GAP,
    FeaturesForWorker,
    HierarchyIndex,
    Interval,
    SmoothError,
    _chunked_seq_reader_from_handle,
    _smooth_one_seq_for_fs,
    merge_adjacent,
    process_seq_chunk,
    smooth_intervals,
    worker_initializer,
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
    """get_ancestors memoises: repeated calls return the same object."""
    assert region_index.get_ancestors("rA") is region_index.get_ancestors("rA")
    assert region_index.get_lca("rA", "rB") == region_index.get_lca("rA", "rB")


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


# --- the multiprocessing pool path -----------------------------------
#
# worker_initializer + process_seq_chunk are how annotate's KMC-backend
# smoothing actually runs; until now only the functions they call had
# direct tests. The oracle is _smooth_one_seq_for_fs applied per
# sequence in-process; the pool path must produce byte-identical lines
# through chunk parsing, worker initialisation in a fresh interpreter
# (spawn), and both output modes.

_POOL_CHUNK_LINES = [
    "seq1\t0\t100\t1\n",
    "seq1\t100\t104\t0\n",  # short novel run flanked by rA: smoothing fodder
    "seq1\t104\t200\t1\n",
    "seq1\t200\t260\t2\n",
    "seq2\t0\t50\t1\n",
    "seq2\t50\t55\t3\n",  # rC island inside rA: promotes toward the LCA
    "seq2\t55\t120\t1\n",
    "seq3\t0\t40\t3\n",
]


def _region_worker_state(region_hierarchy: Hierarchy):
    index = HierarchyIndex.from_hierarchy(region_hierarchy, "region")
    feats = FeaturesForWorker(
        feature_set="region",
        id_to_name={1: "rA", 2: "rB", 3: "rC"},
        novel_label="novel",
    )
    return {"region": index}, {"region": feats}, ["region"]


def _expected_by_seq(region_hierarchy: Hierarchy):
    """The oracle: per-sequence direct smoothing, no chunking, no pool."""
    indices, feats, _ = _region_worker_state(region_hierarchy)
    by_seq: dict[str, list[tuple[int, int, int]]] = {}
    for line in _POOL_CHUNK_LINES:
        seq, start, end, fid = line.split("\t")
        by_seq.setdefault(seq, []).append((int(start), int(end), int(fid)))
    return {
        seq: _smooth_one_seq_for_fs(
            seq_name=seq,
            intervals_raw=raw,
            feature_set="region",
            index=indices["region"],
            features=feats["region"],
        )
        for seq, raw in by_seq.items()
    }


def test_process_seq_chunk_matches_direct_smoothing(
    region_hierarchy: Hierarchy, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Chunk parsing and per-FS orchestration reproduce the direct path,
    and a chunk boundary between sequences changes nothing."""
    import karyoscope.core.smooth as sm

    indices, feats, fs_list = _region_worker_state(region_hierarchy)
    monkeypatch.setattr(sm, "_worker_indices", indices)
    monkeypatch.setattr(sm, "_worker_features_by_fs", feats)
    monkeypatch.setattr(sm, "_worker_feature_sets", fs_list)
    monkeypatch.setattr(sm, "_worker_tmpdir_by_fs", None)

    expected = _expected_by_seq(region_hierarchy)

    one_chunk = process_seq_chunk(list(_POOL_CHUNK_LINES))["region"]
    assert set(one_chunk) == set(expected)
    for seq, (smo, pre) in expected.items():
        assert one_chunk[seq] == (smo, pre)

    split = [
        [line for line in _POOL_CHUNK_LINES if line.startswith("seq1")],
        [line for line in _POOL_CHUNK_LINES if not line.startswith("seq1")],
    ]
    merged: dict[str, object] = {}
    for chunk in split:
        merged.update(process_seq_chunk(chunk)["region"])
    assert merged == one_chunk


def test_process_seq_chunk_assembly_mode_writes_the_same_lines(
    region_hierarchy: Hierarchy, monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    """Assembly mode writes per-sequence files holding exactly the lines
    reads mode returns, and reports their line counts."""
    import karyoscope.core.smooth as sm

    indices, feats, fs_list = _region_worker_state(region_hierarchy)
    fs_dir = tmp_path / "region"
    fs_dir.mkdir()
    monkeypatch.setattr(sm, "_worker_indices", indices)
    monkeypatch.setattr(sm, "_worker_features_by_fs", feats)
    monkeypatch.setattr(sm, "_worker_feature_sets", fs_list)
    monkeypatch.setattr(sm, "_worker_tmpdir_by_fs", {"region": fs_dir})

    counts = process_seq_chunk(list(_POOL_CHUNK_LINES))["region"]
    expected = _expected_by_seq(region_hierarchy)

    for seq, (smo, pre) in expected.items():
        n_pre, n_smo = counts[seq]
        assert (n_pre, n_smo) == (len(pre), len(smo))
        assert (fs_dir / f"{seq}.smo").read_text() == "".join(smo)
        assert (fs_dir / f"{seq}.pre").read_text() == "".join(pre)


@pytest.mark.slow
def test_pool_path_matches_single_process(region_hierarchy: Hierarchy) -> None:
    """The real spawn pool — worker_initializer in a fresh interpreter,
    initargs pickled, chunks distributed — matches the direct path."""
    import multiprocessing as mp

    indices, feats, fs_list = _region_worker_state(region_hierarchy)
    chunks = [
        [line for line in _POOL_CHUNK_LINES if line.startswith(seq)]
        for seq in ("seq1", "seq2", "seq3")
    ]

    ctx = mp.get_context("spawn")
    with ctx.Pool(
        processes=2,
        initializer=worker_initializer,
        initargs=(indices, feats, fs_list, None),
    ) as pool:
        results = pool.map(process_seq_chunk, chunks)

    merged: dict[str, object] = {}
    for r in results:
        merged.update(r["region"])
    expected = _expected_by_seq(region_hierarchy)
    assert merged == {seq: (smo, pre) for seq, (smo, pre) in expected.items()}
