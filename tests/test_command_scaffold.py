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
def test_scaffold_auto_derives_everything(
    cli_runner: CliRunner,
    populated_db_root: Path,
    dummy_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    """One bare ``-i FA`` invocation should produce a complete output set.

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
