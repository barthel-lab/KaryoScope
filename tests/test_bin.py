"""Unit + integration tests for :mod:`karyoscope.core.bin` and ``karyoscope bin``."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest
from click.testing import CliRunner

from karyoscope.cli import main
from karyoscope.core.bin import (
    _best_feature,
    _pick_winner,
    bin_features,
    bin_records,
    leaves_for,
)
from karyoscope.core.io.hierarchy import Hierarchy, HierarchyRow
from karyoscope.exceptions import BinError

# --- _pick_winner ---------------------------------------------------


class TestPickWinner:
    def test_single_feature_wins(self) -> None:
        assert _pick_winner({"A": 50}, bin_size=100) == "A"

    def test_largest_overlap_wins(self) -> None:
        assert _pick_winner({"A": 60, "B": 40}, bin_size=100) == "A"

    def test_novel_majority_wins(self) -> None:
        # 60% novel beats 40% real combined.
        assert _pick_winner({"A": 20, "B": 20, "novel": 60}, bin_size=100) == "novel"

    def test_novel_minority_loses_to_best_real(self) -> None:
        # Novel has the largest overlap but not strict majority — falls back.
        assert _pick_winner({"A": 30, "B": 30, "novel": 40}, bin_size=100) == "A"

    def test_novel_exactly_half_loses(self) -> None:
        # Strict majority threshold: equal-to-half does not count.
        assert _pick_winner({"A": 50, "novel": 50}, bin_size=100) == "A"

    def test_alphabetic_tiebreak_among_reals(self) -> None:
        assert _pick_winner({"B": 50, "A": 50}, bin_size=100) == "A"

    def test_novel_loses_tiebreak_to_real(self) -> None:
        # Equal overlap; novel is deprioritised regardless of name.
        assert _pick_winner({"A": 30, "novel": 30}, bin_size=100) == "A"

    def test_all_novel_returns_novel(self) -> None:
        # Even though novel is deprioritised, with no alternative we keep it.
        assert _pick_winner({"novel": 25}, bin_size=100) == "novel"


# --- _best_feature with leaf set ------------------------------------


class TestBestFeatureLeafPrioritisation:
    def test_no_leaf_set_means_pick_winner(self) -> None:
        # Internal node "centromeric" outvotes leaf "rA" without leaf
        # prioritisation; with leaf_set the leaf still wins.
        counts = {"centromeric": 80, "rA": 20}
        leaves = {"rA", "rB"}
        assert _best_feature(counts, None, 100) == "centromeric"
        assert _best_feature(counts, leaves, 100) == "rA"

    def test_falls_back_to_all_when_no_leaf_present(self) -> None:
        counts = {"centromeric": 80, "novel": 20}
        leaves = {"rA", "rB"}  # neither present in counts
        assert _best_feature(counts, leaves, 100) == "centromeric"

    def test_leaves_compete_under_normal_rules(self) -> None:
        # Two leaves present — leaf-only competition picks the larger.
        counts = {"rA": 60, "rB": 40, "centromeric": 80}
        leaves = {"rA", "rB"}
        assert _best_feature(counts, leaves, 100) == "rA"


# --- leaves_for -----------------------------------------------------


class TestLeavesFor:
    @pytest.fixture
    def region_hierarchy(self) -> Hierarchy:
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

    def test_leaves(self, region_hierarchy: Hierarchy) -> None:
        assert leaves_for(region_hierarchy, "region") == {"rA", "rB", "rC"}

    def test_empty_feature_set(self, region_hierarchy: Hierarchy) -> None:
        assert leaves_for(region_hierarchy, "does_not_exist") == set()


# --- bin_records: pure binning behaviour ----------------------------


class TestBinRecords:
    def test_single_record_one_bin(self) -> None:
        recs = [("chr1", 0, 100, "A")]
        assert list(bin_records(recs, bin_size=100)) == [("chr1", 0, 100, "A")]

    def test_coalesces_adjacent_same_label(self) -> None:
        # Three 100-bp bins all dominated by A should collapse to one row.
        recs = [
            ("chr1", 0, 100, "A"),
            ("chr1", 100, 200, "A"),
            ("chr1", 200, 300, "A"),
        ]
        assert list(bin_records(recs, bin_size=100)) == [("chr1", 0, 300, "A")]

    def test_does_not_coalesce_across_chroms(self) -> None:
        recs = [
            ("chr1", 0, 100, "A"),
            ("chr2", 0, 100, "A"),
        ]
        assert list(bin_records(recs, bin_size=100)) == [
            ("chr1", 0, 100, "A"),
            ("chr2", 0, 100, "A"),
        ]

    def test_truncates_partial_trailing_bin(self) -> None:
        # 0-150 of an A: bin 0 (0-100) is full A, bin 1 (100-200) only
        # covers 100-150 → output should not extend past 150.
        recs = [("chr1", 0, 150, "A")]
        assert list(bin_records(recs, bin_size=100)) == [("chr1", 0, 150, "A")]


class TestRuntTrailingBin:
    """A trailing remainder must not cast a label vote of its own.

    The genome view sets ``bin_size = round(longest_seq / 250)``, so the
    longest sequence is left with a remainder of order tens of bp every
    time. On CHM13 that surfaced as a lone 48 bp ``categorized`` row at the
    end of chr1 -- one that also earned a full karyotype legend entry.
    """

    def test_runt_absorbed_into_previous_row(self) -> None:
        # 10 bp of B trailing a full 100 bp bin of A must not become its own row.
        recs = [("chr1", 0, 100, "A"), ("chr1", 100, 110, "B")]
        assert list(bin_records(recs, bin_size=100)) == [("chr1", 0, 110, "A")]

    def test_runt_absorption_loses_no_bases(self) -> None:
        recs = [("chr1", 0, 100, "A"), ("chr1", 100, 110, "B")]
        out = list(bin_records(recs, bin_size=100))
        assert out[-1][2] == 110

    def test_chm13_chr1_case(self) -> None:
        # The exact observed case: 250 full bins plus a 48 bp remainder whose
        # k-mers stayed ambiguous up to the hierarchy root.
        bin_size = 993_549
        recs = [
            ("chr1", 0, 248_387_250, "chr1"),
            ("chr1", 248_387_250, 248_387_298, "categorized"),
        ]
        assert list(bin_records(recs, bin_size=bin_size)) == [("chr1", 0, 248_387_298, "chr1")]

    def test_trailing_bin_above_cutoff_keeps_its_label(self) -> None:
        # Over half a bin is real signal, not a remainder -- keep it.
        recs = [("chr1", 0, 100, "A"), ("chr1", 100, 180, "B")]
        assert list(bin_records(recs, bin_size=100)) == [
            ("chr1", 0, 100, "A"),
            ("chr1", 100, 180, "B"),
        ]

    def test_lone_runt_still_emitted(self) -> None:
        # No preceding row to absorb into: emit it rather than drop the contig.
        assert list(bin_records([("tig", 0, 40, "A")], bin_size=100)) == [("tig", 0, 40, "A")]

    def test_runt_does_not_bleed_across_chromosomes(self) -> None:
        recs = [
            ("chr1", 0, 100, "A"),
            ("chr1", 100, 110, "Z"),
            ("chr2", 0, 100, "B"),
        ]
        assert list(bin_records(recs, bin_size=100)) == [
            ("chr1", 0, 110, "A"),
            ("chr2", 0, 100, "B"),
        ]

    def test_winning_feature_per_bin(self) -> None:
        recs = [
            ("chr1", 0, 60, "A"),
            ("chr1", 60, 100, "B"),
        ]
        # Bin 0: A=60, B=40 → A.
        assert list(bin_records(recs, bin_size=100)) == [("chr1", 0, 100, "A")]

    def test_leaf_set_prioritises_leaves(self) -> None:
        # Internal node dominates the bin by overlap, but the leaf wins.
        recs = [
            ("chr1", 0, 80, "internal"),
            ("chr1", 80, 100, "leafX"),
        ]
        assert list(bin_records(recs, bin_size=100, leaf_set={"leafX"})) == [
            ("chr1", 0, 100, "leafX")
        ]

    def test_invalid_bin_size_raises(self) -> None:
        with pytest.raises(BinError):
            list(bin_records([("chr1", 0, 1, "A")], bin_size=0))


# --- bin_features: file-to-file -------------------------------------


class TestBinFeaturesFile:
    def _write_bed(
        self, p: Path, rows: list[tuple[str, int, int, str]], *, gz: bool = False
    ) -> None:
        text = "".join(f"{c}\t{s}\t{e}\t{f}\n" for c, s, e, f in rows)
        if gz:
            with gzip.open(p, "wt") as h:
                h.write(text)
        else:
            p.write_text(text)

    def _read_bed(self, p: Path) -> list[tuple[str, int, int, str]]:
        opener = gzip.open if p.suffix == ".gz" else open
        out = []
        with opener(p, "rt") as h:
            for line in h:
                if not line.strip():
                    continue
                c, s, e, f = line.rstrip("\n").split("\t")
                out.append((c, int(s), int(e), f))
        return out

    def test_plain_in_plain_out(self, tmp_path: Path) -> None:
        src = tmp_path / "in.bed"
        dst = tmp_path / "out.bed"
        self._write_bed(src, [("chr1", 0, 300, "A")])
        bin_features(src, dst, bin_size=100)
        assert self._read_bed(dst) == [("chr1", 0, 300, "A")]

    def test_gzip_in_gzip_out(self, tmp_path: Path) -> None:
        src = tmp_path / "in.bed.gz"
        dst = tmp_path / "out.bed.gz"
        self._write_bed(
            src,
            [("chr1", 0, 1000, "X"), ("chr1", 1000, 2000, "novel"), ("chr1", 2000, 3000, "Y")],
            gz=True,
        )
        bin_features(src, dst, bin_size=1000)
        # X dominates bin 0 fully; novel dominates bin 1 fully (covers
        # the full 1000 > 500 majority); Y dominates bin 2 fully.
        assert self._read_bed(dst) == [
            ("chr1", 0, 1000, "X"),
            ("chr1", 1000, 2000, "novel"),
            ("chr1", 2000, 3000, "Y"),
        ]

    def test_bad_input_row_raises(self, tmp_path: Path) -> None:
        src = tmp_path / "bad.bed"
        dst = tmp_path / "out.bed"
        src.write_text("chr1\tnot_a_number\t100\tA\n")
        with pytest.raises(BinError):
            bin_features(src, dst, bin_size=100)


# --- CLI ------------------------------------------------------------


class TestCLI:
    def test_bare_binning_via_stdio(self, cli_runner: CliRunner, tmp_path: Path) -> None:
        src = tmp_path / "in.bed"
        dst = tmp_path / "out.bed"
        src.write_text("chr1\t0\t300\tA\n")
        result = cli_runner.invoke(main, ["bin", "-i", str(src), "-o", str(dst), "-b", "100"])
        assert result.exit_code == 0, result.output
        assert dst.read_text() == "chr1\t0\t300\tA\n"

    def test_db_requires_feature_set(self, cli_runner: CliRunner, tmp_path: Path) -> None:
        src = tmp_path / "in.bed"
        src.write_text("chr1\t0\t100\tA\n")
        result = cli_runner.invoke(
            main, ["bin", "-i", str(src), "-o", "-", "-b", "100", "--db", "X"]
        )
        assert result.exit_code != 0
        assert "feature-set" in result.output.lower()

    def test_bin_size_must_be_positive(self, cli_runner: CliRunner, tmp_path: Path) -> None:
        src = tmp_path / "in.bed"
        src.write_text("chr1\t0\t100\tA\n")
        result = cli_runner.invoke(main, ["bin", "-i", str(src), "-o", "-", "-b", "0"])
        assert result.exit_code != 0

    def test_leaf_prioritisation_via_dummy_db(
        self, cli_runner: CliRunner, populated_db_root: Path, tmp_path: Path
    ) -> None:
        # Build a small input where an internal node has more bp than
        # the leaf, then confirm leaf prioritisation flips the winner.
        # The dummy db's region hierarchy has rA/rB/rC as leaves under
        # aSat/HSat under centromeric.
        src = tmp_path / "in.bed"
        src.write_text("chr1\t0\t80\tcentromeric\nchr1\t80\t100\trA\n")
        dst = tmp_path / "out.bed"
        result = cli_runner.invoke(
            main,
            [
                "bin",
                "-i",
                str(src),
                "-o",
                str(dst),
                "-b",
                "100",
                "--db",
                "KS_dummy_test_v1",
                "--feature-set",
                "region",
            ],
        )
        assert result.exit_code == 0, result.output
        assert dst.read_text() == "chr1\t0\t100\trA\n"


# --- threaded path produces byte-identical output -------------------


@pytest.mark.slow
class TestBinFeaturesThreaded:
    """``bin_features(threads=N)`` should byte-for-byte match
    ``bin_features(threads=1)`` regardless of N.

    Pool-based execution path: chunked-by-sequence-boundary dispatch
    via :class:`multiprocessing.Pool`. The output ordering is
    preserved by ``pool.imap`` (input-chunk-order) and within each
    chunk by the in-worker ``_drive`` loop, so byte equality is the
    right invariant to test.
    """

    def _write_bed(self, p: Path, rows: list[tuple[str, int, int, str]]) -> None:
        p.write_text("".join(f"{c}\t{s}\t{e}\t{f}\n" for c, s, e, f in rows))

    def test_single_sequence_threads_match(self, tmp_path: Path) -> None:
        src = tmp_path / "in.bed"
        self._write_bed(
            src,
            [("chr1", i * 1000, (i + 1) * 1000, "A" if i % 2 == 0 else "B") for i in range(100)],
        )
        dst1 = tmp_path / "out1.bed"
        dst2 = tmp_path / "out2.bed"
        bin_features(src, dst1, bin_size=2000, threads=1)
        bin_features(src, dst2, bin_size=2000, threads=2, chunk_size=20)
        assert dst1.read_text() == dst2.read_text()

    def test_multi_sequence_threads_match(self, tmp_path: Path) -> None:
        # Several sequences forces chunking across sequence boundaries.
        rows: list[tuple[str, int, int, str]] = []
        for chrom in ("chr1", "chr2", "chr3", "chr4"):
            for i in range(50):
                rows.append((chrom, i * 1000, (i + 1) * 1000, "A" if i % 3 else "B"))
        src = tmp_path / "in.bed"
        self._write_bed(src, rows)
        dst1 = tmp_path / "out1.bed"
        dst2 = tmp_path / "out2.bed"
        bin_features(src, dst1, bin_size=10000, threads=1)
        # Force lots of chunks: 4 sequences x 50 lines = 200 lines,
        # chunk_size=30 makes at least 4 chunks (one per sequence).
        bin_features(src, dst2, bin_size=10000, threads=4, chunk_size=30)
        assert dst1.read_text() == dst2.read_text()

    def test_with_leaf_set_threads_match(self, tmp_path: Path) -> None:
        # Worker must receive and apply leaf_set correctly.
        rows = [
            ("chr1", 0, 80, "centromeric"),
            ("chr1", 80, 100, "rA"),
            ("chr2", 0, 100, "rB"),
        ]
        src = tmp_path / "in.bed"
        self._write_bed(src, rows)
        dst1 = tmp_path / "out1.bed"
        dst2 = tmp_path / "out2.bed"
        leaf = {"rA", "rB", "rC"}
        bin_features(src, dst1, bin_size=100, leaf_set=leaf, threads=1)
        bin_features(src, dst2, bin_size=100, leaf_set=leaf, threads=2, chunk_size=1)
        assert dst1.read_text() == dst2.read_text()

    def test_gzip_io_threads_match(self, tmp_path: Path) -> None:
        # Confirm the gzipped-output path works through the pool too.
        rows = [("chr1", i * 100, (i + 1) * 100, "X") for i in range(20)]
        src = tmp_path / "in.bed"
        self._write_bed(src, rows)
        dst1 = tmp_path / "out1.bed.gz"
        dst2 = tmp_path / "out2.bed.gz"
        bin_features(src, dst1, bin_size=200, threads=1)
        bin_features(src, dst2, bin_size=200, threads=2, chunk_size=5)
        # Compare decompressed text -- compressed-byte equality is
        # not guaranteed across runs (gzip metadata).
        with gzip.open(dst1, "rt") as h1, gzip.open(dst2, "rt") as h2:
            assert h1.read() == h2.read()

    def test_threads_zero_auto_detects(self, tmp_path: Path) -> None:
        # threads=0 means auto -- should not raise.
        src = tmp_path / "in.bed"
        self._write_bed(src, [("chr1", 0, 1000, "A")])
        dst = tmp_path / "out.bed"
        bin_features(src, dst, bin_size=100, threads=0)
        assert dst.read_text() == "chr1\t0\t1000\tA\n"

    def test_stdio_forces_single_threaded(self, cli_runner: CliRunner, tmp_path: Path) -> None:
        # The pool path requires a real input path; stdin/stdout
        # invocations should still work with --threads > 1 (we
        # silently fall back to single-threaded).
        src = tmp_path / "in.bed"
        src.write_text("chr1\t0\t300\tA\n")
        result = cli_runner.invoke(main, ["bin", "-i", str(src), "-o", "-", "-b", "100", "-t", "4"])
        assert result.exit_code == 0, result.output
        assert "chr1\t0\t300\tA\n" in result.output
