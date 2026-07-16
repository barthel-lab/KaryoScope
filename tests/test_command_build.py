"""Tests for the ``karyoscope build`` command wiring (no HKS binary needed).

The full construction path is covered by ``test_build_core.py`` (skipped when
``hks`` is absent); here we stub :func:`karyoscope.core.build.build_database` to
exercise spec assembly, flag/spec mutual exclusion, and error surfacing.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

import karyoscope.commands.build as build_cmd
from karyoscope.cli import main
from karyoscope.core.build import BuildResult
from karyoscope.core.buildspec import BuildSpec


@pytest.fixture
def cli_runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def inputs(tmp_path: Path) -> tuple[Path, Path]:
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGTACGTACGT\n")
    bed = tmp_path / "r.bed"
    bed.write_text("chr1\t0\t8\tLINE\n")
    return genome, bed


@pytest.fixture
def stub_build(monkeypatch: pytest.MonkeyPatch) -> dict:
    """Replace build_database with a recorder; return the captured spec/kwargs."""
    captured: dict = {}

    def _fake(spec: BuildSpec, db_root: Path, **kwargs: object) -> BuildResult:
        captured["spec"] = spec
        captured["db_root"] = db_root
        captured["kwargs"] = kwargs
        return BuildResult(
            db_id=spec.id, db_dir=db_root / spec.id, registered=kwargs.get("register") is not False
        )

    monkeypatch.setattr(build_cmd, "build_database", _fake)
    return captured


def test_simple_form_builds_spec(cli_runner: CliRunner, inputs, stub_build, tmp_path: Path) -> None:
    genome, bed = inputs
    result = cli_runner.invoke(
        main,
        [
            "build",
            "--id",
            "HKS_x",
            "--sequence",
            str(genome),
            "--feature-set",
            f"repeat={bed}",
            "--background",
            "repeat=nonrepeat",
            "-s",
            "15",
            "--db-root",
            str(tmp_path / "db"),
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    spec = stub_build["spec"]
    assert spec.id == "HKS_x"
    assert spec.s == 15
    (fs,) = spec.feature_sets
    assert fs.name == "repeat" and fs.background == "nonrepeat"


def test_spec_and_flags_mutually_exclusive(
    cli_runner: CliRunner, inputs, stub_build, tmp_path: Path
) -> None:
    genome, bed = inputs
    spec_file = tmp_path / "spec.yaml"
    spec_file.write_text(
        f'id: X\nversion: "1"\nsequence: {genome}\nfeature_sets:\n  - name: r\n    bed: {bed}\n'
    )
    result = cli_runner.invoke(
        main,
        ["build", "--spec", str(spec_file), "--id", "HKS_x", "--feature-set", f"repeat={bed}"],
    )
    assert result.exit_code != 0
    assert "cannot be combined" in result.output


def test_missing_id_errors(cli_runner: CliRunner, inputs, stub_build) -> None:
    genome, bed = inputs
    result = cli_runner.invoke(
        main, ["build", "--sequence", str(genome), "--feature-set", f"r={bed}"]
    )
    assert result.exit_code != 0
    assert "--id is required" in result.output


def test_no_register_flag_passed_through(
    cli_runner: CliRunner, inputs, stub_build, tmp_path: Path
) -> None:
    genome, bed = inputs
    result = cli_runner.invoke(
        main,
        [
            "build",
            "--id",
            "HKS_x",
            "--sequence",
            str(genome),
            "--feature-set",
            f"repeat={bed}",
            "--no-register",
            "--db-root",
            str(tmp_path / "db"),
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    assert stub_build["kwargs"]["register"] is False


def test_spec_form_reads_yaml(cli_runner: CliRunner, inputs, stub_build, tmp_path: Path) -> None:
    genome, bed = inputs
    spec_file = tmp_path / "spec.yaml"
    spec_file.write_text(
        f'id: HKS_spec\nversion: "3.1.0"\nsequence: {genome}\n'
        f"feature_sets:\n  - name: repeat\n    bed: {bed}\n"
    )
    result = cli_runner.invoke(
        main,
        ["build", "--spec", str(spec_file), "--db-root", str(tmp_path / "db")],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    assert stub_build["spec"].id == "HKS_spec"
    assert stub_build["spec"].version == "3.1.0"
