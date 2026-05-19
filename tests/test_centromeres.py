"""Unit + integration tests for :mod:`karyoscope.core.centromeres`."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from karyoscope.cli import main
from karyoscope.core.centromeres import find_centromere_ranges

# --- pure algorithm -------------------------------------------------


class TestFindCentromereRangesCoarseOnly:
    def test_one_contig_one_centromere_block(self) -> None:
        bins = {
            "chr1_h1_ctgA": [
                (0, 1_000_000, "p_arm"),
                (1_000_000, 2_000_000, "centromere1"),
                (2_000_000, 3_000_000, "centromere1"),
                (3_000_000, 4_000_000, "q_arm"),
            ],
        }
        assert dict(find_centromere_ranges(bins)) == {
            "chr1_h1_ctgA": (1_000_000, 3_000_000),
        }

    def test_contig_without_centromere_is_omitted(self) -> None:
        bins = {
            "chr1_h1_ctgA": [(0, 1_000_000, "p_arm"), (1_000_000, 2_000_000, "q_arm")],
        }
        assert dict(find_centromere_ranges(bins)) == {}

    def test_multiple_contigs_preserves_insertion_order(self) -> None:
        bins = {
            "chr3_h1_ctgC": [(0, 1, "centromere")],
            "chr1_h1_ctgA": [(0, 1, "centromere")],
            "chr2_h1_ctgB": [(0, 1, "centromere")],
        }
        ranges = find_centromere_ranges(bins)
        # Output preserves insertion order from the input.
        assert list(ranges.keys()) == ["chr3_h1_ctgC", "chr1_h1_ctgA", "chr2_h1_ctgB"]

    def test_non_contiguous_centromere_bins_span_outer_bounds(self) -> None:
        # Two centromeric bins separated by a stray q_arm bin -- the
        # output span is min(start) to max(stop), faithful to the
        # archive's coarse-pass behaviour.
        bins = {
            "ctg": [
                (1_000_000, 2_000_000, "centromere"),
                (2_000_000, 3_000_000, "q_arm"),
                (3_000_000, 4_000_000, "centromere"),
            ],
        }
        assert dict(find_centromere_ranges(bins)) == {"ctg": (1_000_000, 4_000_000)}


class TestFindCentromereRangesWithFine:
    def test_fine_narrows_coarse_range(self) -> None:
        coarse = {
            "ctgA": [
                (0, 1_000_000, "p_arm"),
                (1_000_000, 2_000_000, "centromere"),
                (2_000_000, 3_000_000, "centromere"),
                (3_000_000, 4_000_000, "q_arm"),
            ],
        }
        # Fine bins show centromere actually starts at 1.5M and ends at 2.8M.
        fine = {
            "ctgA": [
                (1_400_000, 1_500_000, "p_arm"),
                (1_500_000, 1_600_000, "centromere"),
                (2_700_000, 2_800_000, "centromere"),
                (2_800_000, 2_900_000, "q_arm"),
            ],
        }
        ranges = find_centromere_ranges(coarse, fine)
        assert ranges["ctgA"] == (1_500_000, 2_800_000)

    def test_fine_outside_search_window_ignored(self) -> None:
        # Coarse range is 1M-3M. With default 1Mb buffer the window is
        # 0-4M. A fine centromere bin at 5M-6M is outside; it should be
        # ignored.
        coarse = {"ctgA": [(1_000_000, 3_000_000, "centromere")]}
        fine = {
            "ctgA": [
                (1_500_000, 1_600_000, "centromere"),  # inside window
                (5_000_000, 6_000_000, "centromere"),  # outside, ignored
            ],
        }
        ranges = find_centromere_ranges(coarse, fine)
        # The fine pass picks the in-window bin; the 5-6M bin doesn't
        # widen the range.
        assert ranges["ctgA"] == (1_500_000, 1_600_000)

    def test_fine_with_no_signal_falls_back_to_coarse(self) -> None:
        # Coarse identifies a centromere; fine has bins in the window
        # but none classified as centromere. Output keeps the coarse
        # range so the contig still has a centromere call.
        coarse = {"ctgA": [(1_000_000, 3_000_000, "centromere")]}
        fine = {
            "ctgA": [
                (1_500_000, 1_600_000, "p_arm"),
                (2_400_000, 2_500_000, "q_arm"),
            ],
        }
        ranges = find_centromere_ranges(coarse, fine)
        assert ranges["ctgA"] == (1_000_000, 3_000_000)

    def test_fine_contig_missing_from_fine_bins(self) -> None:
        # Coarse has ctgA but fine doesn't. Output falls back to coarse.
        coarse = {"ctgA": [(1_000_000, 3_000_000, "centromere")]}
        fine: dict[str, list[tuple[int, int, str]]] = {}
        ranges = find_centromere_ranges(coarse, fine)
        assert ranges["ctgA"] == (1_000_000, 3_000_000)

    def test_zero_refinement_buffer_strict_inside(self) -> None:
        coarse = {"ctgA": [(1_000_000, 3_000_000, "centromere")]}
        # A bin at the edge of coarse range is OK; one outside fails.
        fine = {
            "ctgA": [
                (1_500_000, 1_600_000, "centromere"),
                (3_500_000, 3_600_000, "centromere"),  # outside coarse
            ],
        }
        ranges = find_centromere_ranges(coarse, fine, refinement_buffer=0)
        # Only the in-coarse bin counts.
        assert ranges["ctgA"] == (1_500_000, 1_600_000)


# --- CLI parsing (no external tools) --------------------------------


class TestCliSurface:
    def test_help_runs(self, cli_runner: CliRunner) -> None:
        result = cli_runner.invoke(main, ["centromeres", "--help"])
        assert result.exit_code == 0
        assert "centromere coordinates" in result.output.lower()

    def test_no_input_fails(self, cli_runner: CliRunner) -> None:
        result = cli_runner.invoke(main, ["centromeres"])
        assert result.exit_code != 0

    def test_negative_coarse_bin_size_rejected(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
    ) -> None:
        fa = tmp_path / "x.fa"
        fa.write_text(">a\nACGT\n")
        result = cli_runner.invoke(
            main,
            ["centromeres", "-i", str(fa), "--coarse-bin-size", "0"],
        )
        assert result.exit_code != 0

    def test_negative_fine_bin_size_rejected(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
    ) -> None:
        fa = tmp_path / "x.fa"
        fa.write_text(">a\nACGT\n")
        result = cli_runner.invoke(
            main,
            ["centromeres", "-i", str(fa), "--fine-bin-size", "-1"],
        )
        assert result.exit_code != 0


# --- end-to-end against the dummy DB -------------------------------


def _has(name: str) -> bool:
    return shutil.which(name) is not None


_required = pytest.mark.skipif(
    not (_has("seqtk") and _has("bgzip")),
    reason="needs seqtk and bgzip on PATH",
)


@pytest.fixture
def dummy_assembly_fasta(tmp_path: Path) -> Path:
    """Same fixture shape as test_command_scaffold.py."""
    seed = "ACGTGCTAGCTAGGCTATCGTAC"
    fa = tmp_path / "assembly.fa"
    fa.write_text(f">seq_for_chr1\n{seed[:21]}\n>seq_for_chr2\n{seed[2:23]}\n")
    return fa


@_required
@pytest.mark.integration
def test_centromeres_auto_derive_end_to_end(
    cli_runner: CliRunner,
    populated_db_root: Path,
    dummy_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    """Full cascade from FASTA produces a centromeres BED.

    The dummy db's region features (rA, rB, rC) all classify as
    'centromere' under get_simple_region (no 'p_arm'/'q_arm'/'tel'
    substring), so every kept contig gets a centromere call.
    """
    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "centromeres",
            "-i",
            f"hap1={dummy_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--min-scaffold-length",
            "1",
            "--coarse-bin-size",
            "10",
            "--fine-bin-size",
            "0",  # disable fine pass; sequences are tiny
            "--no-bgzip",
        ],
    )
    assert result.exit_code == 0, result.output

    out = out_dir / "assembly.KS_dummy_test_v1.centromeres.bed"
    assert out.is_file()
    lines = [line for line in out.read_text().splitlines() if line]
    assert len(lines) >= 1
    for line in lines:
        parts = line.split("\t")
        assert len(parts) == 3, f"expected 3 columns, got {parts!r}"
        # cols 2 and 3 are integers
        int(parts[1])
        int(parts[2])


@_required
@pytest.mark.integration
def test_no_auto_errors_when_scaffold_missing(
    cli_runner: CliRunner,
    populated_db_root: Path,
    dummy_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    """--no-auto + no pre-existing scaffolded BEDs = clean error."""
    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "centromeres",
            "-i",
            f"hap1={dummy_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--no-auto",
            "--no-bgzip",
        ],
    )
    assert result.exit_code != 0
    assert "missing" in result.output.lower()
