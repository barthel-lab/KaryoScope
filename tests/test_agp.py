"""Unit tests for :mod:`karyoscope.core.io.agp`."""

from __future__ import annotations

from pathlib import Path

from karyoscope.core.io.agp import (
    GAP_LINKAGE,
    GAP_LINKAGE_EVIDENCE,
    GAP_TYPE,
    AgpComponent,
    AgpGap,
    AgpObject,
    write_agp,
)


def _lines(path: Path) -> list[str]:
    return path.read_text().splitlines()


class TestWriteAgp:
    def test_header_first(self, tmp_path: Path) -> None:
        out = tmp_path / "x.agp"
        write_agp([], out)
        assert _lines(out) == ["##agp-version\t2.1"]

    def test_single_component_object(self, tmp_path: Path) -> None:
        out = tmp_path / "x.agp"
        obj = AgpObject(
            name="ctgA",
            parts=[AgpComponent("ctgA", 0, 100, 100, "+")],
        )
        write_agp([obj], out)
        # 1-based inclusive: obj_beg=1, obj_end=100.
        assert _lines(out)[1] == "ctgA\t1\t100\t1\tW\tctgA\t1\t100\t+"

    def test_combined_object_with_gap(self, tmp_path: Path) -> None:
        out = tmp_path / "x.agp"
        # ctgA [0,10), gap [10,15), ctgB [15,23) in 0-based half-open.
        obj = AgpObject(
            name="chr1_hap1",
            parts=[
                AgpComponent("A", 0, 10, 10, "+"),
                AgpGap(10, 15, 5),
                AgpComponent("B", 15, 23, 8, "-"),
            ],
        )
        write_agp([obj], out)
        rows = _lines(out)[1:]
        assert rows[0] == "chr1_hap1\t1\t10\t1\tW\tA\t1\t10\t+"
        assert (
            rows[1]
            == f"chr1_hap1\t11\t15\t2\tN\t5\t{GAP_TYPE}\t{GAP_LINKAGE}\t{GAP_LINKAGE_EVIDENCE}"
        )
        assert rows[2] == "chr1_hap1\t16\t23\t3\tW\tB\t1\t8\t-"

    def test_part_numbers_restart_per_object(self, tmp_path: Path) -> None:
        out = tmp_path / "x.agp"
        objs = [
            AgpObject("o1", [AgpComponent("a", 0, 5, 5, "+")]),
            AgpObject("o2", [AgpComponent("b", 0, 7, 7, "+")]),
        ]
        write_agp(objs, out)
        rows = _lines(out)[1:]
        assert rows[0].split("\t")[3] == "1"
        assert rows[1].split("\t")[3] == "1"

    def test_gap_constants(self) -> None:
        # Lock the KaryoScope gap metadata: alignment-free, but ordering
        # is asserted against a same-genus human reference.
        assert GAP_TYPE == "scaffold"
        assert GAP_LINKAGE == "yes"
        assert GAP_LINKAGE_EVIDENCE == "align_genus"
