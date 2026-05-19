"""Unit tests for :mod:`karyoscope.core.io.scaffold_map`."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.core.io.scaffold_map import (
    MapRow,
    read_map,
    write_legacy_stats,
    write_map,
)
from karyoscope.exceptions import ScaffoldError


def _row(**kw: object) -> MapRow:
    defaults: dict[str, object] = {
        "new_name": "chr1_hap1_ctg1",
        "original_name": "ctg1",
        "input_file": "hap1.fa.gz",
        "hap": "hap1",
        "chromosome": "chr1",
        "flipped": False,
        "length": 1_000_000,
        "stats": "TPCQT",
    }
    defaults.update(kw)
    return MapRow(**defaults)  # type: ignore[arg-type]


class TestRoundTrip:
    def test_basic_roundtrip(self, tmp_path: Path) -> None:
        rows = [
            _row(),
            _row(new_name="chr1_hap2_ctg2_rc", original_name="ctg2", hap="hap2", flipped=True),
            _row(new_name="chr2_hap1_ctg3", original_name="ctg3", chromosome="chr2"),
        ]
        p = tmp_path / "map.tsv"
        write_map(rows, p)
        assert read_map(p) == rows

    def test_header_present(self, tmp_path: Path) -> None:
        p = tmp_path / "map.tsv"
        write_map([_row()], p)
        text = p.read_text()
        first_line = text.splitlines()[0]
        assert first_line.split("\t") == [
            "new_name",
            "original_name",
            "input_file",
            "hap",
            "chromosome",
            "flipped",
            "length",
            "stats",
        ]

    def test_flipped_serialisation(self, tmp_path: Path) -> None:
        p = tmp_path / "map.tsv"
        write_map([_row(flipped=False), _row(flipped=True)], p)
        # Confirm yes/no rather than True/False on disk.
        lines = p.read_text().splitlines()
        assert "\tno\t" in lines[1]
        assert "\tyes\t" in lines[2]


class TestErrors:
    def test_missing_file(self, tmp_path: Path) -> None:
        with pytest.raises(ScaffoldError, match="not found"):
            read_map(tmp_path / "nope.tsv")

    def test_bad_header(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.tsv"
        p.write_text("wrong\tcolumns\n")
        with pytest.raises(ScaffoldError, match="header"):
            read_map(p)

    def test_wrong_column_count(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.tsv"
        p.write_text(
            "new_name\toriginal_name\tinput_file\thap\tchromosome\tflipped\tlength\tstats\n"
            "only_one_column\n"
        )
        with pytest.raises(ScaffoldError, match="expected 8 columns"):
            read_map(p)

    def test_bad_flipped_value(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.tsv"
        p.write_text(
            "new_name\toriginal_name\tinput_file\thap\tchromosome\tflipped\tlength\tstats\n"
            "n\to\tf\th\tc\tmaybe\t1\ts\n"
        )
        with pytest.raises(ScaffoldError, match="flipped"):
            read_map(p)

    def test_bad_length_value(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.tsv"
        p.write_text(
            "new_name\toriginal_name\tinput_file\thap\tchromosome\tflipped\tlength\tstats\n"
            "n\to\tf\th\tc\tno\tnot_a_number\ts\n"
        )
        with pytest.raises(ScaffoldError, match="length"):
            read_map(p)


class TestLegacyStats:
    def test_legacy_stats_two_columns(self, tmp_path: Path) -> None:
        rows = [
            _row(new_name="chr1_hap1_ctg1", stats="TPCQT"),
            _row(new_name="chr2_hap2_ctg2_rc", stats="QC"),
        ]
        p = tmp_path / "stats.tsv"
        write_legacy_stats(rows, p)
        assert p.read_text() == "chr1_hap1_ctg1\tTPCQT\nchr2_hap2_ctg2_rc\tQC\n"
