"""Unit tests for :mod:`karyoscope.core.io.telo`."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

import karyoscope.core.io.telo as telo_mod
from karyoscope.core.io.telo import TeloFlags, parse_telo_file, run_seqtk_telo


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


class TestRunSeqtkTelo:
    def _patch(self, monkeypatch: pytest.MonkeyPatch, captured: list) -> None:
        monkeypatch.setattr(telo_mod, "require_tool", lambda *a, **k: "seqtk")

        def fake_run_tool(cmd: list[str], **kwargs: object) -> SimpleNamespace:
            captured.append(cmd)
            return SimpleNamespace(stdout="seq_a\t0\t100\n")

        monkeypatch.setattr(telo_mod, "run_tool", fake_run_tool)

    def test_default_omits_motif_flag(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        captured: list = []
        self._patch(monkeypatch, captured)
        fasta = tmp_path / "g.fa"
        run_seqtk_telo(fasta, tmp_path / "out.telo")
        assert captured[0] == ["seqtk", "telo", str(fasta)]

    def test_motif_passes_dash_m(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        captured: list = []
        self._patch(monkeypatch, captured)
        fasta = tmp_path / "g.fa"
        out = tmp_path / "out.telo"
        run_seqtk_telo(fasta, out, motif="CCCTAAA")
        assert captured[0] == ["seqtk", "telo", "-m", "CCCTAAA", str(fasta)]
        assert out.read_text() == "seq_a\t0\t100\n"
