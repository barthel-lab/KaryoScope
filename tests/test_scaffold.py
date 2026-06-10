"""Unit tests for :mod:`karyoscope.core.scaffold`."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.core.io.fasta import read_fasta_records
from karyoscope.core.io.scaffold_map import MapRow
from karyoscope.core.io.telo import TeloFlags
from karyoscope.core.scaffold import (
    ContigInput,
    Interval,
    _to_agp_objects,
    assign_main_chromosome,
    category_index,
    chromosome_sort_key,
    classify_and_orient,
    find_largest_contiguous_region,
    flip_bins,
    get_simple_region,
    half_region_totals,
    need_to_flip,
    plan_combined_layout,
    rewrite_bed,
    rewrite_bed_combined,
    rewrite_fasta,
    scaffold_region_majority,
    write_combined_fasta,
)

# --- pure helpers ---------------------------------------------------


class TestGetSimpleRegion:
    @pytest.mark.parametrize(
        ("name", "expected"),
        [
            ("p_arm", "p_arm"),
            ("chr1_p_arm", "p_arm"),
            ("q_arm", "q_arm"),
            ("q_arm_long", "q_arm"),
            ("arm", "arm"),
            ("telomere", "telomere"),
            ("tel_repeat", "telomere"),
            ("novel", "novel"),
            ("centromeric", "centromere"),
            ("alpha_sat", "centromere"),
            ("HSat1A", "centromere"),
        ],
    )
    def test_substring_matching(self, name: str, expected: str) -> None:
        assert get_simple_region(name) == expected


class TestChromosomeSortKey:
    def test_numeric_sorts_first(self) -> None:
        keys = sorted(["chrX", "chr1", "chr10", "chr2"], key=chromosome_sort_key)
        assert keys == ["chr1", "chr2", "chr10", "chrX"]

    def test_handles_no_prefix(self) -> None:
        # Defensive: contigs without "chr" prefix still get a stable order.
        keys = sorted(["1", "2", "10", "X"], key=chromosome_sort_key)
        assert keys == ["1", "2", "10", "X"]


class TestAssignMainChromosome:
    def test_picks_largest_overlap_leaf(self) -> None:
        bins: list[Interval] = [
            (0, 1000, "chr1"),
            (1000, 2000, "chr1"),
            (2000, 2500, "chr2"),
        ]
        assert assign_main_chromosome(bins, {"chr1", "chr2"}) == "chr1"

    def test_ignores_internal_nodes(self) -> None:
        # Internal node "autosome" has huge overlap but isn't a leaf.
        bins: list[Interval] = [
            (0, 5000, "autosome"),
            (5000, 5500, "chr2"),
        ]
        assert assign_main_chromosome(bins, {"chr1", "chr2"}) == "chr2"

    def test_no_leaves_returns_none(self) -> None:
        bins: list[Interval] = [(0, 1000, "novel"), (1000, 2000, "categorized")]
        assert assign_main_chromosome(bins, {"chr1"}) is None


class TestFindLargestContiguousRegion:
    def test_full_match(self) -> None:
        bins = [(0, 1000, "chr1"), (1000, 2000, "chr1")]
        assert find_largest_contiguous_region(bins, "chr1", {"chr1", "chr2"}) == (0, 2000)

    def test_picks_longest_compatible_block(self) -> None:
        # Two compatible blocks separated by an incompatible chr2 block.
        # Compatible runs: chr1 (0-1000) then chr1+chr1 (2000-4000).
        # Second one is longer → it wins.
        bins = [
            (0, 1000, "chr1"),
            (1000, 2000, "chr2"),  # incompatible
            (2000, 3000, "chr1"),
            (3000, 4000, "chr1"),
        ]
        assert find_largest_contiguous_region(bins, "chr1", {"chr1", "chr2"}) == (2000, 4000)

    def test_internal_nodes_count_as_compatible(self) -> None:
        bins = [
            (0, 1000, "chr1"),
            (1000, 2000, "autosome"),  # internal node, compatible
            (2000, 3000, "chr1"),
        ]
        assert find_largest_contiguous_region(bins, "chr1", {"chr1", "chr2"}) == (0, 3000)

    def test_no_main_returns_full_extent(self) -> None:
        bins = [(0, 1000, "novel"), (1000, 2000, "categorized")]
        assert find_largest_contiguous_region(bins, None, {"chr1"}) == (0, 2000)


class TestHalfRegionTotals:
    def test_basic_split(self) -> None:
        # window = [0, 1000); breakpoint = 500.
        bins = [(0, 500, "p_arm_short"), (500, 1000, "q_arm_short")]
        out = half_region_totals(bins, 0, 1000)
        assert out["p_arm"] == [500, 0]
        assert out["q_arm"] == [0, 500]

    def test_clamps_to_window(self) -> None:
        # A bin that extends past the window only counts the overlap.
        bins = [(0, 2000, "p_arm")]
        out = half_region_totals(bins, 500, 1500)
        # p_arm spans 500-1500, breakpoint = 1000 → first=500, second=500.
        assert out["p_arm"] == [500, 500]


class TestScaffoldRegionMajority:
    def test_majority(self) -> None:
        bins = [(0, 1000, "p_arm"), (1000, 1100, "q_arm")]
        assert scaffold_region_majority(bins, 0, 1100) == "p_arm"

    def test_empty(self) -> None:
        assert scaffold_region_majority([], 0, 0) == ""


class TestNeedToFlip:
    def _halfs(self, bins: list[Interval], length: int) -> dict[str, list[int]]:
        return half_region_totals(bins, 0, length)

    def test_correct_orientation_no_tels(self) -> None:
        # p_arm then centromere then q_arm: should NOT flip.
        bins = [
            (0, 3_000_000, "p_arm_short"),
            (3_000_000, 5_000_000, "centromere1"),
            (5_000_000, 10_000_000, "q_arm_long"),
        ]
        assert not need_to_flip(
            bins,
            self._halfs(bins, 10_000_000),
            region_start=0,
            region_end=10_000_000,
            scaffold_length=10_000_000,
            telo=TeloFlags(False, False),
            is_acrocentric=False,
        )

    def test_reversed_orientation_no_tels(self) -> None:
        # q_arm then centromere then p_arm: should flip.
        bins = [
            (0, 5_000_000, "q_arm_long"),
            (5_000_000, 7_000_000, "centromere1"),
            (7_000_000, 10_000_000, "p_arm_short"),
        ]
        assert need_to_flip(
            bins,
            self._halfs(bins, 10_000_000),
            region_start=0,
            region_end=10_000_000,
            scaffold_length=10_000_000,
            telo=TeloFlags(False, False),
            is_acrocentric=False,
        )

    def test_both_tels_centroid_correct(self) -> None:
        # Both ends telomered; p centroid < q centroid → don't flip.
        bins = [
            (0, 2_000_000, "p_arm"),
            (8_000_000, 10_000_000, "q_arm"),
        ]
        assert not need_to_flip(
            bins,
            self._halfs(bins, 10_000_000),
            region_start=0,
            region_end=10_000_000,
            scaffold_length=10_000_000,
            telo=TeloFlags(True, True),
            is_acrocentric=False,
        )

    def test_both_tels_centroid_reversed(self) -> None:
        # Both ends telomered; p centroid > q centroid → flip.
        bins = [
            (0, 2_000_000, "q_arm"),
            (8_000_000, 10_000_000, "p_arm"),
        ]
        assert need_to_flip(
            bins,
            self._halfs(bins, 10_000_000),
            region_start=0,
            region_end=10_000_000,
            scaffold_length=10_000_000,
            telo=TeloFlags(True, True),
            is_acrocentric=False,
        )

    def test_acrocentric_pure_telomere_with_stop_tel(self) -> None:
        # Acrocentric, only telomere content, only stop tel touches the
        # region end → flip (so the telomere ends up at the start).
        bins = [(9_500_000, 10_000_000, "telomere1")]
        assert need_to_flip(
            bins,
            self._halfs(bins, 10_000_000),
            region_start=0,
            region_end=10_000_000,
            scaffold_length=10_000_000,
            telo=TeloFlags(False, True),
            is_acrocentric=True,
        )


class TestFlipBins:
    def test_mirror_and_reverse(self) -> None:
        bins = [(0, 100, "A"), (200, 300, "B")]
        # length = 400; reversed → [(200,300,B), (0,100,A)]
        # mirrored → [(100, 200, B), (300, 400, A)]
        assert flip_bins(bins, 400) == [(100, 200, "B"), (300, 400, "A")]


class TestCategoryIndex:
    def test_start_tel(self) -> None:
        assert (
            category_index(
                p_total=0,
                cen_total=0,
                q_total=0,
                has_start_tel=True,
                has_stop_tel=False,
            )
            == 0
        )

    def test_stop_tel(self) -> None:
        assert (
            category_index(
                p_total=0,
                cen_total=0,
                q_total=0,
                has_start_tel=False,
                has_stop_tel=True,
            )
            == 7
        )

    def test_p_heavy_no_cen(self) -> None:
        assert (
            category_index(
                p_total=10,
                cen_total=0,
                q_total=5,
                has_start_tel=False,
                has_stop_tel=False,
            )
            == 1
        )

    def test_balanced(self) -> None:
        assert (
            category_index(
                p_total=10,
                cen_total=5,
                q_total=10,
                has_start_tel=False,
                has_stop_tel=False,
            )
            == 3
        )


# --- end-to-end on synthetic input ---------------------------------


def _ci(
    name: str,
    *,
    hap: str = "hap1",
    chrom_bins: list[Interval],
    region_bins: list[Interval],
    telo: TeloFlags = TeloFlags(False, False),
    input_file: str = "in.fa",
) -> ContigInput:
    length = max(
        (stop for _, stop, _ in chrom_bins + region_bins),
        default=0,
    )
    return ContigInput(
        input_name=hap,
        input_file=input_file,
        contig_name=name,
        length=length,
        chromosome_bins=chrom_bins,
        region_bins=region_bins,
        telo=telo,
    )


class TestClassifyAndOrient:
    def test_basic_keeps_long_contig(self) -> None:
        c = _ci(
            "ctgA",
            chrom_bins=[(0, 10_000_000, "chr1")],
            region_bins=[
                (0, 3_000_000, "p_arm"),
                (3_000_000, 5_000_000, "centromere1"),
                (5_000_000, 10_000_000, "q_arm"),
            ],
            telo=TeloFlags(True, True),
        )
        rows = classify_and_orient([c], chromosome_leaves={"chr1", "chr2"})
        assert len(rows) == 1
        r = rows[0]
        assert r.new_name == "chr1_hap1_ctgA"
        assert r.chromosome == "chr1"
        assert r.flipped is False
        assert r.stats == "TPCQT"

    def test_drops_short_contig_without_telo(self) -> None:
        c = _ci(
            "tinyA",
            chrom_bins=[(0, 100_000, "chr1")],
            region_bins=[(0, 100_000, "p_arm")],
        )
        rows = classify_and_orient([c], chromosome_leaves={"chr1"})
        assert rows == []

    def test_keeps_short_contig_with_telo(self) -> None:
        c = _ci(
            "tinyB",
            chrom_bins=[(0, 100_000, "chr1")],
            region_bins=[(0, 100_000, "p_arm")],
            telo=TeloFlags(True, False),
        )
        rows = classify_and_orient([c], chromosome_leaves={"chr1"})
        assert len(rows) == 1

    def test_drops_contig_with_no_leaf_chrom(self) -> None:
        c = _ci(
            "novel_only",
            chrom_bins=[(0, 10_000_000, "novel")],
            region_bins=[(0, 10_000_000, "novel")],
        )
        rows = classify_and_orient([c], chromosome_leaves={"chr1"})
        assert rows == []

    def test_flipped_emits_rc_suffix(self) -> None:
        c = _ci(
            "ctgRev",
            chrom_bins=[(0, 10_000_000, "chr1")],
            region_bins=[
                (0, 5_000_000, "q_arm"),
                (5_000_000, 7_000_000, "centromere1"),
                (7_000_000, 10_000_000, "p_arm"),
            ],
        )
        rows = classify_and_orient([c], chromosome_leaves={"chr1"})
        assert len(rows) == 1
        assert rows[0].new_name.endswith("_rc")
        assert rows[0].flipped is True

    def test_orders_by_chromosome_then_hap(self) -> None:
        # Two contigs on different chromosomes, both haps present.
        cs = [
            _ci(
                "c2h2",
                hap="hap2",
                chrom_bins=[(0, 10_000_000, "chr2")],
                region_bins=[(0, 10_000_000, "p_arm")],
                telo=TeloFlags(True, False),
            ),
            _ci(
                "c1h1",
                hap="hap1",
                chrom_bins=[(0, 10_000_000, "chr1")],
                region_bins=[(0, 10_000_000, "p_arm")],
                telo=TeloFlags(True, False),
            ),
            _ci(
                "c1h2",
                hap="hap2",
                chrom_bins=[(0, 10_000_000, "chr1")],
                region_bins=[(0, 10_000_000, "p_arm")],
                telo=TeloFlags(True, False),
            ),
        ]
        rows = classify_and_orient(cs, chromosome_leaves={"chr1", "chr2"})
        # chr1 before chr2; within chr1, hap1 before hap2.
        assert [r.new_name for r in rows] == [
            "chr1_hap1_c1h1",
            "chr1_hap2_c1h2",
            "chr2_hap2_c2h2",
        ]


# --- BED rewriter ----------------------------------------------------


class TestRewriteBed:
    def test_no_flip(self, tmp_path: Path) -> None:
        src = tmp_path / "in.bed"
        src.write_text("ctgA\t0\t100\tX\nctgA\t100\t200\tY\n")
        dst = tmp_path / "out.bed"
        rows = [MapRow("chr1_hap1_ctgA", "ctgA", "in.fa", "hap1", "chr1", False, 200, "PQ")]
        rewrite_bed(src, dst, map_rows=rows)
        assert dst.read_text() == "chr1_hap1_ctgA\t0\t100\tX\nchr1_hap1_ctgA\t100\t200\tY\n"

    def test_flip_mirrors_and_reverses(self, tmp_path: Path) -> None:
        src = tmp_path / "in.bed"
        src.write_text("ctgA\t0\t100\tP\nctgA\t100\t300\tC\nctgA\t300\t500\tQ\n")
        dst = tmp_path / "out.bed"
        rows = [MapRow("chr1_hap1_ctgA_rc", "ctgA", "in.fa", "hap1", "chr1", True, 500, "QCP")]
        rewrite_bed(src, dst, map_rows=rows)
        assert dst.read_text() == (
            "chr1_hap1_ctgA_rc\t0\t200\tQ\n"
            "chr1_hap1_ctgA_rc\t200\t400\tC\n"
            "chr1_hap1_ctgA_rc\t400\t500\tP\n"
        )

    def test_contigs_not_in_map_are_dropped(self, tmp_path: Path) -> None:
        src = tmp_path / "in.bed"
        src.write_text("ctgA\t0\t100\tX\nctgUnmapped\t0\t50\tZ\n")
        dst = tmp_path / "out.bed"
        rows = [MapRow("chr1_hap1_ctgA", "ctgA", "in.fa", "hap1", "chr1", False, 100, "P")]
        rewrite_bed(src, dst, map_rows=rows)
        assert "ctgUnmapped" not in dst.read_text()

    def test_contigs_in_map_but_absent_from_bed_are_silently_skipped(
        self,
        tmp_path: Path,
    ) -> None:
        src = tmp_path / "in.bed"
        src.write_text("ctgA\t0\t100\tX\n")
        dst = tmp_path / "out.bed"
        rows = [
            MapRow("chr1_hap1_ctgA", "ctgA", "in.fa", "hap1", "chr1", False, 100, "P"),
            MapRow("chr1_hap1_ctgGhost", "ctgGhost", "in.fa", "hap1", "chr1", False, 100, "P"),
        ]
        rewrite_bed(src, dst, map_rows=rows)
        assert dst.read_text() == "chr1_hap1_ctgA\t0\t100\tX\n"


# --- FASTA rewriter --------------------------------------------------


class TestRewriteFasta:
    def _write_src(self, p: Path, records: dict[str, str]) -> None:
        p.write_text("".join(f">{n}\n{s}\n" for n, s in records.items()))

    def test_renames_without_flip(self, tmp_path: Path) -> None:
        src = tmp_path / "in.fa"
        self._write_src(src, {"ctgA": "ACGT", "ctgB": "GGCC"})
        dst = tmp_path / "out.fa"
        rows = [
            MapRow("chr1_hap1_ctgA", "ctgA", "in.fa", "hap1", "chr1", False, 4, "P"),
            MapRow("chr1_hap1_ctgB", "ctgB", "in.fa", "hap1", "chr1", False, 4, "P"),
        ]
        rewrite_fasta(src, dst, map_rows=rows, keep_unscaffolded=False)
        out = read_fasta_records(dst)
        assert list(out.keys()) == ["chr1_hap1_ctgA", "chr1_hap1_ctgB"]
        assert out["chr1_hap1_ctgA"] == "ACGT"
        assert out["chr1_hap1_ctgB"] == "GGCC"

    def test_reverse_complement_when_flipped(self, tmp_path: Path) -> None:
        src = tmp_path / "in.fa"
        self._write_src(src, {"ctgA": "AATTCG"})  # RC = CGAATT
        dst = tmp_path / "out.fa"
        rows = [
            MapRow("chr1_hap1_ctgA_rc", "ctgA", "in.fa", "hap1", "chr1", True, 6, "Q"),
        ]
        rewrite_fasta(src, dst, map_rows=rows, keep_unscaffolded=False)
        assert read_fasta_records(dst)["chr1_hap1_ctgA_rc"] == "CGAATT"

    def test_keep_unscaffolded_appends_originals(self, tmp_path: Path) -> None:
        # ctgA is scaffolded; tiny is left over.
        src = tmp_path / "in.fa"
        self._write_src(src, {"ctgA": "AAAA", "tiny": "TT"})
        dst = tmp_path / "out.fa"
        rows = [
            MapRow("chr1_hap1_ctgA", "ctgA", "in.fa", "hap1", "chr1", False, 4, "P"),
        ]
        rewrite_fasta(src, dst, map_rows=rows, keep_unscaffolded=True)
        out = read_fasta_records(dst)
        # Scaffolded contigs come first (in map order), then unscaffolded
        # under their ORIGINAL names.
        assert list(out.keys()) == ["chr1_hap1_ctgA", "tiny"]
        assert out["tiny"] == "TT"

    def test_drop_unscaffolded(self, tmp_path: Path) -> None:
        src = tmp_path / "in.fa"
        self._write_src(src, {"ctgA": "AAAA", "tiny": "TT"})
        dst = tmp_path / "out.fa"
        rows = [
            MapRow("chr1_hap1_ctgA", "ctgA", "in.fa", "hap1", "chr1", False, 4, "P"),
        ]
        rewrite_fasta(src, dst, map_rows=rows, keep_unscaffolded=False)
        out = read_fasta_records(dst)
        assert list(out.keys()) == ["chr1_hap1_ctgA"]
        assert "tiny" not in out

    def test_emits_in_map_order_not_source_order(self, tmp_path: Path) -> None:
        # Source: A, B; map emits B first (e.g. because of chrom ordering).
        src = tmp_path / "in.fa"
        self._write_src(src, {"ctgA": "AAAA", "ctgB": "CCCC"})
        dst = tmp_path / "out.fa"
        rows = [
            MapRow("chr1_h1_ctgB", "ctgB", "in.fa", "h1", "chr1", False, 4, "P"),
            MapRow("chr2_h1_ctgA", "ctgA", "in.fa", "h1", "chr2", False, 4, "P"),
        ]
        rewrite_fasta(src, dst, map_rows=rows, keep_unscaffolded=False)
        assert list(read_fasta_records(dst).keys()) == ["chr1_h1_ctgB", "chr2_h1_ctgA"]

    def test_map_contig_absent_from_fasta_is_skipped(self, tmp_path: Path) -> None:
        # ctgB is in the map but not in the source FASTA. Skip silently.
        src = tmp_path / "in.fa"
        self._write_src(src, {"ctgA": "AAAA"})
        dst = tmp_path / "out.fa"
        rows = [
            MapRow("chr1_h1_ctgA", "ctgA", "in.fa", "h1", "chr1", False, 4, "P"),
            MapRow("chr1_h1_ctgGhost", "ctgGhost", "in.fa", "h1", "chr1", False, 4, "P"),
        ]
        rewrite_fasta(src, dst, map_rows=rows, keep_unscaffolded=False)
        out = read_fasta_records(dst)
        assert list(out.keys()) == ["chr1_h1_ctgA"]

    def test_gzip_input_and_output(self, tmp_path: Path) -> None:
        import gzip as _gz

        src = tmp_path / "in.fa.gz"
        with _gz.open(src, "wt") as h:
            h.write(">ctgA\nACGT\n")
        dst = tmp_path / "out.fa.gz"
        rows = [MapRow("chr1_h1_ctgA", "ctgA", "in.fa.gz", "h1", "chr1", False, 4, "P")]
        rewrite_fasta(src, dst, map_rows=rows, keep_unscaffolded=False)
        assert read_fasta_records(dst)["chr1_h1_ctgA"] == "ACGT"


# --- combined-chromosome layout ------------------------------------


def _row(
    new_name: str,
    original: str,
    chrom: str,
    hap: str,
    *,
    flipped: bool = False,
    length: int = 100,
) -> MapRow:
    return MapRow(new_name, original, "in.fa", hap, chrom, flipped, length, "P")


class TestPlanCombinedLayout:
    def test_two_contigs_one_object_with_offsets(self) -> None:
        # E_A=7 (L=10), E_B=5 (L=8), gap=5 -> B starts at 10+5=15.
        rows = [
            _row("chr1_hap1_A", "A", "chr1", "hap1", length=7),
            _row("chr1_hap1_B", "B", "chr1", "hap1", length=5),
        ]
        true_lengths = {"A": 10, "B": 8}
        objs = plan_combined_layout(
            rows, true_lengths, gap_size=5, acrocentrics=set(), combine_acrocentrics=False
        )
        assert len(objs) == 1
        obj = objs[0]
        assert obj.name == "chr1_hap1"
        assert obj.combined is True
        assert [c.object_start for c in obj.components] == [0, 15]
        assert [c.true_length for c in obj.components] == [10, 8]
        assert [c.bed_extent for c in obj.components] == [7, 5]

    def test_separate_haps_are_separate_objects(self) -> None:
        rows = [
            _row("chr1_hap1_A", "A", "chr1", "hap1"),
            _row("chr1_hap2_B", "B", "chr1", "hap2"),
        ]
        objs = plan_combined_layout(rows, {"A": 100, "B": 100}, gap_size=10, acrocentrics=set())
        assert [o.name for o in objs] == ["chr1_hap1", "chr1_hap2"]

    def test_acrocentric_not_combined_by_default(self) -> None:
        rows = [
            _row("chr13_hap1_A", "A", "chr13", "hap1"),
            _row("chr13_hap1_B", "B", "chr13", "hap1"),
        ]
        objs = plan_combined_layout(
            rows,
            {"A": 100, "B": 100},
            gap_size=10,
            acrocentrics={"chr13"},
            combine_acrocentrics=False,
        )
        # Two singleton objects, each under its encoded name, no combining.
        assert [o.name for o in objs] == ["chr13_hap1_A", "chr13_hap1_B"]
        assert all(o.combined is False for o in objs)
        assert all(len(o.components) == 1 for o in objs)

    def test_acrocentric_combined_when_forced(self) -> None:
        rows = [
            _row("chr13_hap1_A", "A", "chr13", "hap1"),
            _row("chr13_hap1_B", "B", "chr13", "hap1"),
        ]
        objs = plan_combined_layout(
            rows,
            {"A": 100, "B": 100},
            gap_size=10,
            acrocentrics={"chr13"},
            combine_acrocentrics=True,
        )
        assert [o.name for o in objs] == ["chr13_hap1"]
        assert objs[0].combined is True

    def test_single_contig_group_renamed(self) -> None:
        rows = [_row("chr2_hap1_A", "A", "chr2", "hap1")]
        objs = plan_combined_layout(rows, {"A": 100}, gap_size=10, acrocentrics=set())
        assert objs[0].name == "chr2_hap1"
        assert objs[0].combined is True
        assert len(objs[0].components) == 1

    def test_contig_absent_from_fasta_skipped(self) -> None:
        rows = [
            _row("chr1_hap1_A", "A", "chr1", "hap1", length=7),
            _row("chr1_hap1_Ghost", "Ghost", "chr1", "hap1", length=5),
        ]
        objs = plan_combined_layout(rows, {"A": 10}, gap_size=5, acrocentrics=set())
        assert len(objs) == 1
        assert [c.row.original_name for c in objs[0].components] == ["A"]


# --- combined FASTA -------------------------------------------------


class TestWriteCombinedFasta:
    def _write_src(self, p: Path, records: dict[str, str]) -> None:
        p.write_text("".join(f">{n}\n{s}\n" for n, s in records.items()))

    def test_concatenates_with_gap(self, tmp_path: Path) -> None:
        rows = [
            _row("chr1_hap1_A", "A", "chr1", "hap1", length=7),
            _row("chr1_hap1_B", "B", "chr1", "hap1", length=5),
        ]
        records = {"A": "ACGTACGTAC", "B": "GGGGCCCC"}  # len 10, 8
        objs = plan_combined_layout(rows, {"A": 10, "B": 8}, gap_size=5, acrocentrics=set())
        out = tmp_path / "out.fa"
        leftovers = write_combined_fasta(records, objs, out, keep_unscaffolded=False)
        res = read_fasta_records(out)
        assert list(res.keys()) == ["chr1_hap1"]
        assert res["chr1_hap1"] == "ACGTACGTAC" + "N" * 5 + "GGGGCCCC"
        assert leftovers == []

    def test_flipped_component_reverse_complemented(self, tmp_path: Path) -> None:
        rows = [_row("chr1_hap1_A_rc", "A", "chr1", "hap1", flipped=True, length=3)]
        records = {"A": "AATTCG"}  # RC = CGAATT
        objs = plan_combined_layout(rows, {"A": 6}, gap_size=5, acrocentrics=set())
        out = tmp_path / "out.fa"
        write_combined_fasta(records, objs, out, keep_unscaffolded=False)
        assert read_fasta_records(out)["chr1_hap1"] == "CGAATT"

    def test_keep_unscaffolded_returns_leftovers(self, tmp_path: Path) -> None:
        rows = [_row("chr1_hap1_A", "A", "chr1", "hap1", length=7)]
        records = {"A": "ACGTACGTAC", "tiny": "TT"}
        objs = plan_combined_layout(rows, {"A": 10}, gap_size=5, acrocentrics=set())
        out = tmp_path / "out.fa"
        leftovers = write_combined_fasta(records, objs, out, keep_unscaffolded=True)
        res = read_fasta_records(out)
        assert list(res.keys()) == ["chr1_hap1", "tiny"]
        assert leftovers == [("tiny", 2)]

    def test_acrocentric_singletons_kept_separate(self, tmp_path: Path) -> None:
        rows = [
            _row("chr13_hap1_A", "A", "chr13", "hap1", length=7),
            _row("chr13_hap1_B", "B", "chr13", "hap1", length=5),
        ]
        records = {"A": "ACGTACGTAC", "B": "GGGGCCCC"}
        objs = plan_combined_layout(rows, {"A": 10, "B": 8}, gap_size=5, acrocentrics={"chr13"})
        out = tmp_path / "out.fa"
        write_combined_fasta(records, objs, out, keep_unscaffolded=False)
        res = read_fasta_records(out)
        assert list(res.keys()) == ["chr13_hap1_A", "chr13_hap1_B"]
        assert res["chr13_hap1_A"] == "ACGTACGTAC"


# --- combined BED ---------------------------------------------------


class TestRewriteBedCombined:
    def _layout(self, gap_size: int = 5):
        rows = [
            _row("chr1_hap1_A", "A", "chr1", "hap1", length=7),
            _row("chr1_hap1_B", "B", "chr1", "hap1", length=5),
        ]
        return plan_combined_layout(rows, {"A": 10, "B": 8}, gap_size=gap_size, acrocentrics=set())

    def test_shift_and_novel_gap(self, tmp_path: Path) -> None:
        # A tiles [0,7), B tiles [0,5). gap=5 -> B offset 15.
        src = tmp_path / "in.bed"
        src.write_text("A\t0\t7\tchr1\nB\t0\t5\tchr1\n")
        out = tmp_path / "out.bed"
        rewrite_bed_combined(src, out, objects=self._layout(), gzip_out=False)
        lines = [tuple(x.split("\t")) for x in out.read_text().splitlines()]
        assert lines == [
            ("chr1_hap1", "0", "7", "chr1"),
            ("chr1_hap1", "7", "15", "novel"),  # gap = G + (L-E) = 5 + 3 = 8
            ("chr1_hap1", "15", "20", "chr1"),
        ]

    def test_gap_length_is_G_plus_k_minus_1(self, tmp_path: Path) -> None:
        src = tmp_path / "in.bed"
        src.write_text("A\t0\t7\tchr1\nB\t0\t5\tchr1\n")
        out = tmp_path / "out.bed"
        rewrite_bed_combined(src, out, objects=self._layout(gap_size=100), gzip_out=False)
        gap = next(
            x.split("\t") for x in out.read_text().splitlines() if x.split("\t")[3] == "novel"
        )
        # gap spans E_A(7)+offset0 .. offset_B; length = 100 + (10-7) = 103.
        assert int(gap[2]) - int(gap[1]) == 103

    def test_adjacent_novel_merges_across_junction(self, tmp_path: Path) -> None:
        # A ends novel, B starts novel: the trailing novel, the gap, and
        # the leading novel coalesce into a single run.
        src = tmp_path / "in.bed"
        src.write_text("A\t0\t4\tchr1\nA\t4\t7\tnovel\nB\t0\t2\tnovel\nB\t2\t5\tchr1\n")
        out = tmp_path / "out.bed"
        rewrite_bed_combined(src, out, objects=self._layout(), gzip_out=False)
        lines = [tuple(x.split("\t")) for x in out.read_text().splitlines()]
        assert lines == [
            ("chr1_hap1", "0", "4", "chr1"),
            ("chr1_hap1", "4", "17", "novel"),  # 4..7 + gap 7..15 + 15..17
            ("chr1_hap1", "17", "20", "chr1"),
        ]

    def test_flipped_component_mirrored(self, tmp_path: Path) -> None:
        rows = [_row("chr1_hap1_A_rc", "A", "chr1", "hap1", flipped=True, length=7)]
        objs = plan_combined_layout(rows, {"A": 10}, gap_size=5, acrocentrics=set())
        src = tmp_path / "in.bed"
        # A: [0,2) p, [2,7) q. Flipped within [0,7): q->[0,5), p->[5,7).
        src.write_text("A\t0\t2\tp\nA\t2\t7\tq\n")
        out = tmp_path / "out.bed"
        rewrite_bed_combined(src, out, objects=objs, gzip_out=False)
        lines = [tuple(x.split("\t")) for x in out.read_text().splitlines()]
        assert lines == [
            ("chr1_hap1", "0", "5", "q"),
            ("chr1_hap1", "5", "7", "p"),
        ]

    def test_missing_contig_filled_with_novel(self, tmp_path: Path) -> None:
        # B has no records in this feature set: fill its extent with novel.
        src = tmp_path / "in.bed"
        src.write_text("A\t0\t7\tchr1\n")
        out = tmp_path / "out.bed"
        rewrite_bed_combined(src, out, objects=self._layout(), gzip_out=False)
        lines = [tuple(x.split("\t")) for x in out.read_text().splitlines()]
        # A [0,7) chr1 ; novel gap 7..15 ; B extent filled novel 15..20 -> merges to 7..20.
        assert lines == [
            ("chr1_hap1", "0", "7", "chr1"),
            ("chr1_hap1", "7", "20", "novel"),
        ]


# --- AGP from layout ------------------------------------------------


class TestToAgpObjects:
    def test_combined_object_has_gap_rows(self) -> None:
        rows = [
            _row("chr1_hap1_A", "A", "chr1", "hap1", length=7),
            _row("chr1_hap1_B", "B", "chr1", "hap1", flipped=True, length=5),
        ]
        objs = plan_combined_layout(rows, {"A": 10, "B": 8}, gap_size=5, acrocentrics=set())
        agp = _to_agp_objects(objs, [])
        assert len(agp) == 1
        parts = agp[0].parts
        # W(A) +, N gap, W(B) -.
        assert parts[0].component_id == "A"
        assert parts[0].orientation == "+"
        assert parts[0].object_end == 10
        assert parts[1].length == 5  # literal N gap, not G+k-1
        assert parts[1].object_start == 10 and parts[1].object_end == 15
        assert parts[2].component_id == "B"
        assert parts[2].orientation == "-"
        assert parts[2].object_start == 15 and parts[2].object_end == 23

    def test_leftovers_become_singleton_objects(self) -> None:
        rows = [_row("chr1_hap1_A", "A", "chr1", "hap1", length=7)]
        objs = plan_combined_layout(rows, {"A": 10}, gap_size=5, acrocentrics=set())
        agp = _to_agp_objects(objs, [("tiny", 2)])
        assert [o.name for o in agp] == ["chr1_hap1", "tiny"]
        assert agp[1].parts[0].length == 2
