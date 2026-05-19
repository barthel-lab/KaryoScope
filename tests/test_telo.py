"""Unit tests for :mod:`karyoscope.core.io.telo`."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.core.io.telo import TeloFlags, parse_telo_file


class TestParseTeloFile:
    def test_basic_three_seqs(self, tmp_path: Path) -> None:
        p = tmp_path / "x.telo"
        p.write_text("seq_a\t0\t500\nseq_b\t100000\t100500\nseq_c\t0\t500\nseq_c\t99500\t100000\n")
        flags = parse_telo_file(p)
        assert flags == {
            "seq_a": TeloFlags(start=True, stop=False),
            "seq_b": TeloFlags(start=False, stop=True),
            "seq_c": TeloFlags(start=True, stop=True),
        }

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            parse_telo_file(tmp_path / "missing.telo")

    def test_empty_file(self, tmp_path: Path) -> None:
        p = tmp_path / "empty.telo"
        p.write_text("")
        assert parse_telo_file(p) == {}

    def test_skips_blank_and_comment_lines(self, tmp_path: Path) -> None:
        p = tmp_path / "x.telo"
        p.write_text("# header comment\n\nseq_a\t0\t100\n")
        assert parse_telo_file(p) == {"seq_a": TeloFlags(start=True, stop=False)}

    def test_non_integer_start_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "x.telo"
        p.write_text("seq_a\tnot_a_number\t100\n")
        with pytest.raises(ValueError, match="start"):
            parse_telo_file(p)

    def test_multiple_starts_collapse_to_one_flag(self, tmp_path: Path) -> None:
        # Two rows starting at 0 on the same seq should yield start=True,
        # not break.
        p = tmp_path / "x.telo"
        p.write_text("seq_a\t0\t100\nseq_a\t0\t200\n")
        assert parse_telo_file(p) == {"seq_a": TeloFlags(start=True, stop=False)}
