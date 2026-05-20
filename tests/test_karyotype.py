"""Unit + integration tests for :mod:`karyoscope.core.karyotype`."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import ClassVar

import pytest
from click.testing import CliRunner

from karyoscope.cli import main
from karyoscope.core.io.scaffold_map import MapRow
from karyoscope.core.karyotype import (
    PREDEFINED_SEX_SYSTEMS,
    RenderInput,
    get_expected_haps,
    render_karyotype,
)
from karyoscope.exceptions import KaryotypeError

# --- get_expected_haps ----------------------------------------------


class TestGetExpectedHaps:
    haps: ClassVar[list[str]] = ["hap1", "hap2"]

    def test_autosome_returns_all_haps(self) -> None:
        assert get_expected_haps("chr1", "male", self.haps, "XY") == self.haps
        assert get_expected_haps("chr1", "female", self.haps, "XY") == self.haps

    def test_xy_male_chrY_first_hap_only(self) -> None:
        assert get_expected_haps("chrY", "male", self.haps, "XY") == ["hap1"]

    def test_xy_male_chrX_second_hap_only(self) -> None:
        # chrX is at index 1 in ['chrY', 'chrX']; with two haps that maps
        # to ['hap2'].
        assert get_expected_haps("chrX", "male", self.haps, "XY") == ["hap2"]

    def test_xy_female_chrY_empty(self) -> None:
        assert get_expected_haps("chrY", "female", self.haps, "XY") == []

    def test_xy_female_chrX_both_haps(self) -> None:
        assert get_expected_haps("chrX", "female", self.haps, "XY") == self.haps

    def test_unknown_sex_no_expectation_for_sex_chroms(self) -> None:
        assert get_expected_haps("chrY", None, self.haps, "XY") == []
        assert get_expected_haps("chr1", None, self.haps, "XY") == self.haps

    def test_reference_uses_first_hap(self) -> None:
        assert get_expected_haps("chr1", "reference", self.haps, "XY") == ["hap1"]

    def test_zw_female_chrW_first_hap(self) -> None:
        # ZW: heterogametic_sex='female', sex_chromosomes=['chrZ', 'chrW'].
        # Index of chrW is 1, so first hap.
        assert get_expected_haps("chrW", "female", self.haps, "ZW") == ["hap2"]
        assert get_expected_haps("chrZ", "female", self.haps, "ZW") == ["hap1"]
        assert get_expected_haps("chrZ", "male", self.haps, "ZW") == self.haps
        assert get_expected_haps("chrW", "male", self.haps, "ZW") == []

    def test_unknown_system_raises(self) -> None:
        with pytest.raises(KaryotypeError, match="unknown sex-determination"):
            get_expected_haps("chr1", "male", self.haps, "ABC")

    def test_dict_system_accepted(self) -> None:
        # Custom dict matching the predefined schema.
        custom = PREDEFINED_SEX_SYSTEMS["XY"]
        assert get_expected_haps("chrY", "male", self.haps, custom) == ["hap1"]


# --- render_karyotype: unit-level (no external tools) ---------------


def _row(
    name: str,
    *,
    chrom: str,
    hap: str,
    stats: str = "PCQ",
    flipped: bool = False,
    length: int = 1000,
) -> MapRow:
    return MapRow(
        new_name=name,
        original_name=name.split("_")[-1],
        input_file="x.fa",
        hap=hap,
        chromosome=chrom,
        flipped=flipped,
        length=length,
        stats=stats,
    )


def _read_svg_text(p: Path) -> str:
    return p.read_text()


class TestRenderKaryotypeUnit:
    def test_renders_genome_mode_minimal(self, tmp_path: Path) -> None:
        # One contig, chr1, hap1, genome mode.
        ri = RenderInput(
            map_rows=[_row("chr1_hap1_a", chrom="chr1", hap="hap1", stats="TPCQT", length=1000)],
            binned_bed={
                "chr1_hap1_a": [(0, 500, "rA"), (500, 1000, "rB")],
            },
        )
        out = tmp_path / "out.svg"
        render_karyotype(
            [ri],
            colors={"rA": "#2ca02c", "rB": "#d62728"},
            mode="genome",
            output_path=out,
        )
        assert out.is_file()
        text = _read_svg_text(out)
        assert text.startswith("<?xml") or text.startswith("<svg")
        assert "<svg" in text
        # Should contain at least one rectangle (for the rendered bins
        # plus the background plus the sequence outline).
        assert text.count("<rect") >= 2
        # Colours appear in the output.
        assert "#2ca02c" in text or "#d62728" in text

    def test_renders_unknown_mode_raises(self, tmp_path: Path) -> None:
        ri = RenderInput(map_rows=[], binned_bed={})
        with pytest.raises(KaryotypeError, match="unknown mode"):
            render_karyotype(
                [ri],
                colors={},
                mode="wat",  # type: ignore[arg-type]
                output_path=tmp_path / "x.svg",
            )

    def test_centromere_mode_requires_centromere_ranges(self, tmp_path: Path) -> None:
        ri = RenderInput(
            map_rows=[_row("chr1_hap1_a", chrom="chr1", hap="hap1", length=1000)],
            binned_bed={"chr1_hap1_a": [(0, 500, "rA")]},
            centromere_ranges=None,
        )
        with pytest.raises(KaryotypeError, match="centromere mode requires"):
            render_karyotype(
                [ri],
                colors={"rA": "#000"},
                mode="centromere",
                output_path=tmp_path / "x.svg",
            )

    def test_centromere_mode_with_ranges_renders(self, tmp_path: Path) -> None:
        ri = RenderInput(
            map_rows=[_row("chr1_hap1_a", chrom="chr1", hap="hap1", length=2000)],
            binned_bed={"chr1_hap1_a": [(800, 1200, "rA")]},
            centromere_ranges={"chr1_hap1_a": (800, 1200)},
        )
        out = tmp_path / "x.svg"
        render_karyotype(
            [ri],
            colors={"rA": "#abcdef"},
            mode="centromere",
            output_path=out,
        )
        assert out.is_file()
        assert "#abcdef" in out.read_text()

    def test_subtelomere_mode_skips_non_telomere_contigs(self, tmp_path: Path) -> None:
        ri = RenderInput(
            map_rows=[
                # Has telomere
                _row("chr1_hap1_a", chrom="chr1", hap="hap1", stats="TPCQ", length=1_000_000),
                # No telomere
                _row("chr1_hap1_b", chrom="chr1", hap="hap1", stats="PCQ", length=1_000_000),
            ],
            binned_bed={
                "chr1_hap1_a": [(0, 100_000, "rA")],
                "chr1_hap1_b": [(0, 100_000, "rB")],
            },
        )
        out = tmp_path / "x.svg"
        render_karyotype(
            [ri],
            colors={"rA": "#abcdef", "rB": "#fedcba"},
            mode="subtelomere",
            output_path=out,
        )
        text = out.read_text()
        # The contig with a telomere gets rendered; the one without doesn't.
        assert "#abcdef" in text
        assert "#fedcba" not in text

    def test_unknown_feature_renders_white(self, tmp_path: Path) -> None:
        # Missing colour -> warning + white. We don't have a hook to
        # observe the warning here, but the function must not raise.
        ri = RenderInput(
            map_rows=[_row("chr1_hap1_a", chrom="chr1", hap="hap1", length=1000)],
            binned_bed={"chr1_hap1_a": [(0, 1000, "mystery_feature")]},
        )
        out = tmp_path / "x.svg"
        render_karyotype([ri], colors={}, mode="genome", output_path=out)
        assert out.is_file()


# --- CLI parsing ---------------------------------------------------


class TestCliSurface:
    def test_help_runs(self, cli_runner: CliRunner) -> None:
        result = cli_runner.invoke(main, ["karyotype", "--help"])
        assert result.exit_code == 0
        assert "render" in result.output.lower()

    def test_no_input_fails(self, cli_runner: CliRunner) -> None:
        result = cli_runner.invoke(main, ["karyotype"])
        assert result.exit_code != 0

    def test_outdir_and_output_conflict(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
    ) -> None:
        fa = tmp_path / "x.fa"
        fa.write_text(">a\nACGT\n")
        result = cli_runner.invoke(
            main,
            [
                "karyotype",
                "-i",
                str(fa),
                "--outdir",
                str(tmp_path / "a"),
                "--output",
                str(tmp_path / "b.svg"),
            ],
        )
        assert result.exit_code != 0
        assert "--outdir and --output" in result.output


# --- integration tests against the dummy DB ------------------------


def _has(name: str) -> bool:
    return shutil.which(name) is not None


_required = pytest.mark.skipif(
    not (_has("seqtk") and _has("bgzip")),
    reason="needs seqtk and bgzip on PATH",
)


@pytest.fixture
def dummy_assembly_fasta(tmp_path: Path) -> Path:
    seed = "ACGTGCTAGCTAGGCTATCGTAC"
    fa = tmp_path / "assembly.fa"
    fa.write_text(f">seq_for_chr1\n{seed[:21]}\n>seq_for_chr2\n{seed[2:23]}\n")
    return fa


@_required
@pytest.mark.integration
def test_karyotype_genome_mode_against_dummy_db(
    cli_runner: CliRunner,
    populated_db_root: Path,
    dummy_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "karyotype",
            "-i",
            f"hap1={dummy_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--mode",
            "genome",
            "--bin-size",
            "10",  # tiny -- sequences are 21 bp
            "--min-scaffold-length",
            "1",
            "--feature-set",
            "region",
            "--no-human-chroms",  # dummy db has only chr1, chr2
        ],
    )
    assert result.exit_code == 0, result.output

    svgs = list(out_dir.glob("*.karyotype.svg"))
    assert svgs, f"no SVGs produced; got: {list(out_dir.iterdir())}"
    text = svgs[0].read_text()
    assert "<svg" in text
    assert "</svg>" in text
    # Colour from the dummy db's region set.
    assert "#2ca02c" in text or "#d62728" in text or "#9467bd" in text


@_required
@pytest.mark.integration
def test_karyotype_centromere_mode_against_dummy_db(
    cli_runner: CliRunner,
    populated_db_root: Path,
    dummy_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "karyotype",
            "-i",
            f"hap1={dummy_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--mode",
            "centromere",
            "--bin-size",
            "10",
            "--min-scaffold-length",
            "1",
            "--feature-set",
            "region",
            "--no-human-chroms",
        ],
    )
    # Centromere mode requires non-empty centromere ranges across all
    # inputs. The dummy db's rA/rB/rC region features all classify as
    # centromere via get_simple_region, so we should get a non-empty SVG.
    assert result.exit_code == 0, result.output
    svgs = list(out_dir.glob("*.centromere.region.karyotype.svg"))
    assert svgs
