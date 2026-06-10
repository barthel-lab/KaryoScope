"""Integration tests for ``karyoscope scaffold``.

Marked ``@pytest.mark.integration`` because the end-to-end runs need
the C++ ``get_featureIDs`` binary (via the auto-run ``annotate``
step), ``bgzip`` for the default output codepath, and ``seqtk`` for
the telomere detection step.

Unit-level behaviour of the algorithm lives in ``tests/test_scaffold.py``.
Unit-level CLI argument parsing lives below in the bottom section
(those tests don't need any external tools).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from karyoscope.cli import main
from karyoscope.commands.scaffold import _parse_named_path
from karyoscope.core.io.scaffold_map import read_map

# --- CLI parsing (no external tools) --------------------------------


class TestParseNamedPath:
    def test_name_and_path(self) -> None:
        name, path = _parse_named_path("hap1=foo.fa")
        assert name == "hap1"
        assert path == Path("foo.fa")

    def test_bare_path(self) -> None:
        name, path = _parse_named_path("foo.fa")
        assert name is None
        assert path == Path("foo.fa")

    def test_empty_name_rejected(self) -> None:
        import click as _click

        with pytest.raises(_click.BadParameter):
            _parse_named_path("=foo.fa")


class TestCliSurface:
    def test_help_runs(self, cli_runner: CliRunner) -> None:
        result = cli_runner.invoke(main, ["scaffold", "--help"])
        assert result.exit_code == 0
        assert "Order, orient" in result.output

    def test_no_input_fails(self, cli_runner: CliRunner) -> None:
        result = cli_runner.invoke(main, ["scaffold"])
        assert result.exit_code != 0

    def test_telo_without_name_rejected(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
    ) -> None:
        fa = tmp_path / "x.fa"
        fa.write_text(">a\nACGT\n")
        # --telo requires NAME=PATH form
        result = cli_runner.invoke(
            main,
            ["scaffold", "-i", f"hap1={fa}", "--telo", str(tmp_path / "x.telo")],
        )
        assert result.exit_code != 0
        assert "NAME=PATH" in result.output

    def test_unmatched_telo_name_rejected(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
    ) -> None:
        fa = tmp_path / "x.fa"
        fa.write_text(">a\nACGT\n")
        result = cli_runner.invoke(
            main,
            [
                "scaffold",
                "-i",
                f"hap1={fa}",
                "--telo",
                f"unmatched={tmp_path / 'x.telo'}",
            ],
        )
        assert result.exit_code != 0
        assert "no matching --input" in result.output


# --- end-to-end against the dummy DB -------------------------------


def _has(name: str) -> bool:
    return shutil.which(name) is not None


# get_featureIDs is found by ``karyoscope.core.io.kmc`` (which walks
# from the package source up to ``native/get_featureIDs/build/...``),
# so it doesn't need to be on $PATH. Only the tools we invoke via
# ``require_tool`` need to be there.
pytestmark_required = pytest.mark.skipif(
    not (_has("seqtk") and _has("bgzip")),
    reason="needs seqtk and bgzip on PATH",
)


@pytest.fixture
def dummy_assembly_fasta(tmp_path: Path) -> Path:
    """A tiny FASTA whose three sequences land on dummy-db features.

    Sequences chosen so that:
    * ``seq_for_chr1`` contains the dummy db's featureID-1 21-mer (chr1/rA);
    * ``seq_novel_chr2`` mostly novel with a featureID-3 21-mer
      embedded (chr2/rC), to give classify_and_orient something on chr2.

    The dummy db's seed is ``ACGTGCTAGCTAGGCTATCGTAC`` per
    ``tests/data/build_dummy_db.py``. The relevant 21-mers are:
    pos 0 (chr1/rA), pos 1 (chr1/rB), pos 2 (chr2/rC).
    """
    seed = "ACGTGCTAGCTAGGCTATCGTAC"
    fa = tmp_path / "assembly.fa"
    fa.write_text(
        f">seq_for_chr1\n{seed[:21]}\n"  # one 21-mer → featureID 1 (chr1/rA)
        f">seq_for_chr2\n{seed[2:23]}\n"  # one 21-mer → featureID 3 (chr2/rC)
    )
    return fa


@pytestmark_required
@pytest.mark.integration
def test_scaffold_auto_derives_everything_bed_mode(
    cli_runner: CliRunner,
    populated_db_root: Path,
    dummy_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    """``--mode bed`` should produce per-feature-set scaffolded BEDs.

    The dummy db's contigs are far below the default
    ``--min-scaffold-length`` (5 Mb) so we lower the threshold so the
    map file is non-empty and we can confirm the rewrite path.
    """
    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "scaffold",
            "-i",
            f"hap1={dummy_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--mode",
            "bed",
            "--min-scaffold-length",
            "1",
            "--bin-size",
            "10",  # tiny — sequences are ~21 bp
            "--no-bgzip",
        ],
    )
    assert result.exit_code == 0, result.output

    # Expect map + stats + per-feature-set rewritten BED.
    map_path = out_dir / "assembly.KS_dummy_test_v1.scaffold_map.tsv"
    stats_path = out_dir / "assembly.KS_dummy_test_v1.scaffold_stats.tsv"
    assert map_path.is_file()
    assert stats_path.is_file()

    rows = read_map(map_path)
    assert len(rows) >= 1  # at least one contig classified.
    names = {r.new_name for r in rows}
    # Each entry's encoded name starts with its chromosome.
    for r in rows:
        assert r.new_name.startswith(f"{r.chromosome}_{r.hap}_")
    # All names should be unique.
    assert len(names) == len(rows)

    # The per-input map + legacy stats agree on names.
    stats_lines = stats_path.read_text().splitlines()
    assert len(stats_lines) == len(rows)
    stat_names = {line.split("\t", 1)[0] for line in stats_lines}
    assert stat_names == names

    # A scaffolded BED exists for each feature set in the manifest.
    scaffolded = sorted(out_dir.glob("assembly.KS_dummy_test_v1.*.smoothed.scaffolded.bed"))
    assert scaffolded, f"no scaffolded BEDs produced; got: {list(out_dir.iterdir())}"

    # In bed mode, no scaffolded FASTA should be written.
    assert not list(out_dir.glob("assembly.KS_dummy_test_v1.scaffolded.fa*"))


@pytestmark_required
@pytest.mark.integration
def test_no_auto_errors_when_inputs_missing(
    cli_runner: CliRunner,
    populated_db_root: Path,
    dummy_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    """--no-auto + no pre-existing annotation BEDs = clean error."""
    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "scaffold",
            "-i",
            f"hap1={dummy_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--no-auto",
            "--no-bgzip",
        ],
    )
    assert result.exit_code != 0
    assert "annotation BEDs" in result.output.lower() or "missing" in result.output.lower()


@pytestmark_required
@pytest.mark.integration
def test_map_file_is_authoritative_format(
    cli_runner: CliRunner,
    populated_db_root: Path,
    dummy_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    """Map file round-trips cleanly via read_map."""
    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "scaffold",
            "-i",
            f"hap1={dummy_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--min-scaffold-length",
            "1",
            "--bin-size",
            "10",
            "--no-bgzip",
        ],
    )
    assert result.exit_code == 0, result.output
    map_path = out_dir / "assembly.KS_dummy_test_v1.scaffold_map.tsv"
    rows = read_map(map_path)
    for r in rows:
        # input_file is a basename (not the full path) per the design.
        assert "/" not in r.input_file
        assert r.input_file == dummy_assembly_fasta.name
        # The encoded name conforms to chrom_hap_contig[_rc].
        suffix = "_rc" if r.flipped else ""
        assert r.new_name == f"{r.chromosome}_{r.hap}_{r.original_name}{suffix}"


# --- FASTA mode (Stage 5d-1b) ---------------------------------------


def test_feature_set_with_mode_fasta_is_rejected(
    cli_runner: CliRunner,
    tmp_path: Path,
) -> None:
    """--feature-set has no effect in FASTA mode and should error cleanly."""
    fa = tmp_path / "x.fa"
    fa.write_text(">a\nACGT\n")
    result = cli_runner.invoke(
        main,
        [
            "scaffold",
            "-i",
            f"hap1={fa}",
            "--mode",
            "fasta",
            "--feature-set",
            "chromosome",
        ],
    )
    assert result.exit_code != 0
    assert "feature-set" in result.output.lower()
    assert "fasta" in result.output.lower()


@pytestmark_required
@pytest.mark.integration
def test_fasta_mode_writes_scaffolded_fasta(
    cli_runner: CliRunner,
    populated_db_root: Path,
    dummy_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    """Default (--mode fasta) writes one scaffolded FASTA per input, no BEDs."""
    from karyoscope.core.io.fasta import read_fasta_records

    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "scaffold",
            "-i",
            f"hap1={dummy_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--min-scaffold-length",
            "1",
            "--bin-size",
            "10",
            "--no-bgzip",
        ],
    )
    assert result.exit_code == 0, result.output

    # Map + stats always written.
    assert (out_dir / "assembly.KS_dummy_test_v1.scaffold_map.tsv").is_file()
    # Scaffolded FASTA exists.
    fa_out = out_dir / "assembly.KS_dummy_test_v1.scaffolded.fa"
    assert fa_out.is_file(), f"missing scaffolded FASTA; got: {list(out_dir.iterdir())}"
    # No scaffolded BEDs in fasta mode.
    assert not list(out_dir.glob("*.smoothed.scaffolded.bed*"))

    # Output FASTA names match the map's new_name column for placed contigs.
    rows = read_map(out_dir / "assembly.KS_dummy_test_v1.scaffold_map.tsv")
    fa_records = read_fasta_records(fa_out)
    for r in rows:
        assert r.new_name in fa_records


@pytestmark_required
@pytest.mark.integration
def test_both_mode_writes_fasta_and_beds(
    cli_runner: CliRunner,
    populated_db_root: Path,
    dummy_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    """--mode both produces both per-feature-set BEDs and the FASTA."""
    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "scaffold",
            "-i",
            f"hap1={dummy_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--mode",
            "both",
            "--min-scaffold-length",
            "1",
            "--bin-size",
            "10",
            "--no-bgzip",
        ],
    )
    assert result.exit_code == 0, result.output
    assert (out_dir / "assembly.KS_dummy_test_v1.scaffolded.fa").is_file()
    assert list(out_dir.glob("*.smoothed.scaffolded.bed"))


@pytestmark_required
@pytest.mark.integration
def test_drop_unscaffolded_omits_leftover_contigs(
    cli_runner: CliRunner,
    populated_db_root: Path,
    tmp_path: Path,
) -> None:
    """--drop-unscaffolded leaves out contigs that aren't in the map."""
    from karyoscope.core.io.fasta import read_fasta_records

    # Build a FASTA with one matching contig and one tiny novel decoy.
    seed = "ACGTGCTAGCTAGGCTATCGTAC"
    fa = tmp_path / "asm.fa"
    fa.write_text(
        f">scaffolded_one\n{seed[:21]}\n"  # featureID 1 -> chr1/rA
        f">decoy_junk\n{'AAAAA' * 4}\n"  # all-novel, gets dropped
    )
    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "scaffold",
            "-i",
            f"hap1={fa}",
            "--outdir",
            str(out_dir),
            "--min-scaffold-length",
            "1",
            "--bin-size",
            "10",
            "--no-bgzip",
            "--drop-unscaffolded",
        ],
    )
    assert result.exit_code == 0, result.output
    fa_out = read_fasta_records(out_dir / "asm.KS_dummy_test_v1.scaffolded.fa")
    # decoy_junk should not appear under its original name.
    assert "decoy_junk" not in fa_out


# --- combine-chromosomes (Stage 5d) ---------------------------------


def test_combine_chromosomes_mode_bed_rejected(
    cli_runner: CliRunner,
    tmp_path: Path,
) -> None:
    """--combine-chromosomes needs a FASTA output; --mode bed errors cleanly."""
    fa = tmp_path / "x.fa"
    fa.write_text(">a\nACGT\n")
    result = cli_runner.invoke(
        main,
        ["scaffold", "-i", f"hap1={fa}", "--mode", "bed", "--combine-chromosomes"],
    )
    assert result.exit_code != 0
    assert "combine-chromosomes" in result.output.lower()
    assert "bed" in result.output.lower()


def test_negative_scaffold_gap_size_rejected(
    cli_runner: CliRunner,
    tmp_path: Path,
) -> None:
    fa = tmp_path / "x.fa"
    fa.write_text(">a\nACGT\n")
    result = cli_runner.invoke(
        main,
        ["scaffold", "-i", f"hap1={fa}", "--combine-chromosomes", "--scaffold-gap-size", "-1"],
    )
    assert result.exit_code != 0
    assert "scaffold-gap-size" in result.output.lower()


def test_help_lists_combine_flags(cli_runner: CliRunner) -> None:
    result = cli_runner.invoke(main, ["scaffold", "--help"])
    assert result.exit_code == 0
    assert "--combine-chromosomes" in result.output
    assert "--scaffold-gap-size" in result.output
    assert "--combine-acrocentrics" in result.output


@pytest.fixture
def two_chr1_assembly_fasta(tmp_path: Path) -> Path:
    """A FASTA with two chr1 contigs and one chr2 contig.

    Lets the combine path concatenate the two chr1 contigs into a single
    ``chr1_hap1`` record with an N gap between them. 21-mer positions in
    the dummy seed: pos 0 -> featureID 1 (chr1/rA), pos 1 -> featureID 2
    (chr1/rB), pos 2 -> featureID 3 (chr2/rC).
    """
    seed = "ACGTGCTAGCTAGGCTATCGTAC"
    fa = tmp_path / "asm.fa"
    fa.write_text(
        f">chr1_part_a\n{seed[0:21]}\n"  # chr1/rA
        f">chr1_part_b\n{seed[1:22]}\n"  # chr1/rB
        f">chr2_part\n{seed[2:23]}\n"  # chr2/rC
    )
    return fa


@pytestmark_required
@pytest.mark.integration
def test_combine_chromosomes_fasta_and_agp(
    cli_runner: CliRunner,
    populated_db_root: Path,
    two_chr1_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    """--combine-chromosomes writes a combined FASTA + AGP, with N gaps."""
    from karyoscope.core.io.fasta import read_fasta_records

    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "scaffold",
            "-i",
            f"hap1={two_chr1_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--mode",
            "fasta",
            "--combine-chromosomes",
            "--scaffold-gap-size",
            "10",
            "--min-scaffold-length",
            "1",
            "--bin-size",
            "10",
            "--no-bgzip",
        ],
    )
    assert result.exit_code == 0, result.output

    base = "asm.KS_dummy_test_v1"
    fa_out = out_dir / f"{base}.scaffolded.combined_chromosomes.fa"
    agp_out = out_dir / f"{base}.scaffolded.combined_chromosomes.agp"
    assert fa_out.is_file(), list(out_dir.iterdir())
    assert agp_out.is_file()
    # The plain (non-combined) scaffolded FASTA must NOT be written.
    assert not (out_dir / f"{base}.scaffolded.fa").exists()

    recs = read_fasta_records(fa_out)
    # The two chr1 contigs combined into one chr1_hap1 record with an N run.
    assert "chr1_hap1" in recs
    assert "N" * 10 in recs["chr1_hap1"]
    # chr2 is its own (single-contig) object, renamed, no N gap.
    assert "chr2_hap1" in recs
    assert "N" not in recs["chr2_hap1"]

    agp = agp_out.read_text()
    assert agp.startswith("##agp-version\t2.1")
    # Exactly one N (gap) row, for the chr1 junction.
    n_rows = [ln for ln in agp.splitlines() if ln.split("\t")[4:5] == ["N"]]
    assert len(n_rows) == 1
    assert n_rows[0].split("\t")[5] == "10"  # gap_length
    assert n_rows[0].split("\t")[8] == "align_genus"


@pytestmark_required
@pytest.mark.integration
def test_combine_both_writes_only_combined_beds(
    cli_runner: CliRunner,
    populated_db_root: Path,
    two_chr1_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    """--mode both --combine-chromosomes writes combined BEDs, not plain ones."""
    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "scaffold",
            "-i",
            f"hap1={two_chr1_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--mode",
            "both",
            "--combine-chromosomes",
            "--scaffold-gap-size",
            "10",
            "--min-scaffold-length",
            "1",
            "--bin-size",
            "10",
            "--no-bgzip",
        ],
    )
    assert result.exit_code == 0, result.output
    assert list(out_dir.glob("*.smoothed.scaffolded.combined_chromosomes.bed"))
    # The plain scaffolded BEDs must NOT be produced in a combine run.
    plain = [
        p for p in out_dir.glob("*.smoothed.scaffolded.bed") if "combined_chromosomes" not in p.name
    ]
    assert not plain
