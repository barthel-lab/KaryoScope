"""Smoke tests for the top-level CLI and all subcommands."""

from __future__ import annotations

import pytest
from click.testing import CliRunner

from karyoscope import __version__
from karyoscope.cli import main


def test_top_level_help_lists_all_commands(cli_runner: CliRunner) -> None:
    """`karyoscope --help` should mention every registered subcommand."""
    result = cli_runner.invoke(main, ["--help"])
    assert result.exit_code == 0
    for cmd_name in (
        "download",
        "annotate",
        "scaffold",
        "bin",
        "centromeres",
        "karyotype",
        "info",
        "version",
    ):
        assert cmd_name in result.output, f"`{cmd_name}` missing from --help output"


def test_version_flag(cli_runner: CliRunner) -> None:
    """`karyoscope --version` should print the package version."""
    result = cli_runner.invoke(main, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


@pytest.mark.parametrize(
    "subcommand",
    [
        "download",
        "annotate",
        "scaffold",
        "bin",
        "centromeres",
        "karyotype",
        "info",
        "version",
    ],
)
def test_subcommand_help_succeeds(cli_runner: CliRunner, subcommand: str) -> None:
    """Each subcommand should respond to --help without error."""
    result = cli_runner.invoke(main, [subcommand, "--help"])
    assert result.exit_code == 0, f"`karyoscope {subcommand} --help` failed: {result.output}"


@pytest.mark.parametrize(
    "subcommand",
    [
        "centromeres",
        "karyotype",
    ],
)
def test_stub_subcommands_exit_cleanly(cli_runner: CliRunner, subcommand: str) -> None:
    """Stubbed subcommands should exit non-zero with a clear message,
    not raise an unhandled exception."""
    result = cli_runner.invoke(main, [subcommand])
    assert result.exit_code != 0
    assert "not yet implemented" in result.output.lower()
    # No traceback should leak to the user (click.ClickException is clean).
    assert "Traceback" not in result.output


def test_version_command_runs(cli_runner: CliRunner) -> None:
    """`karyoscope version` should run cleanly and report the version."""
    result = cli_runner.invoke(main, ["version"])
    assert result.exit_code == 0
    assert __version__ in result.output
    assert "Python" in result.output


def test_info_command_runs(cli_runner: CliRunner) -> None:
    """`karyoscope info` (no args) should run cleanly and mention the db root."""
    result = cli_runner.invoke(main, ["info"])
    assert result.exit_code == 0
    assert "database root" in result.output.lower()
