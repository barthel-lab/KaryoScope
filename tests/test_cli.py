"""Smoke tests for the top-level CLI and all subcommands."""

from __future__ import annotations

from pathlib import Path

import click
import pytest
from click.testing import CliRunner

from karyoscope import __version__
from karyoscope.cli import _LAZY_COMMANDS, main


def test_top_level_help_lists_all_commands(cli_runner: CliRunner) -> None:
    """`karyoscope --help` should mention every registered subcommand."""
    result = cli_runner.invoke(main, ["--help"])
    assert result.exit_code == 0
    for cmd_name in _LAZY_COMMANDS:
        assert cmd_name in result.output, f"`{cmd_name}` missing from --help output"


def test_version_flag(cli_runner: CliRunner) -> None:
    """`karyoscope --version` should print the package version."""
    result = cli_runner.invoke(main, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


@pytest.mark.parametrize("subcommand", sorted(_LAZY_COMMANDS))
def test_subcommand_help_succeeds(cli_runner: CliRunner, subcommand: str) -> None:
    """Each subcommand should respond to --help without error."""
    result = cli_runner.invoke(main, [subcommand, "--help"])
    assert result.exit_code == 0, f"`karyoscope {subcommand} --help` failed: {result.output}"


# (The Stage 5d series turned every subcommand into a real
# implementation; the old test_stub_subcommands_exit_cleanly check
# has been removed because there are no stubs left.)


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


# --- --no-resource-check / --no-space-check ---------------------------


def test_resource_check_flag_prefers_the_new_name() -> None:
    from karyoscope.commands._options import resolve_resource_check_flag

    assert resolve_resource_check_flag(True, False, command="annotate") is True
    assert resolve_resource_check_flag(False, False, command="annotate") is False


def test_the_deprecated_alias_still_works_and_warns(caplog: pytest.LogCaptureFixture) -> None:
    """--no-space-check shipped in 2.1.0 and is documented, so it cannot just vanish.

    The README promises deprecations ship with warnings and a back-compatible
    transition; this is that transition.
    """
    from karyoscope.commands._options import resolve_resource_check_flag

    with caplog.at_level("WARNING"):
        assert resolve_resource_check_flag(False, True, command="annotate") is True
    assert "deprecated" in caplog.text
    assert "--no-resource-check" in caplog.text


def test_supplying_both_spellings_is_a_usage_error() -> None:
    from karyoscope.commands._options import resolve_resource_check_flag

    with pytest.raises(click.UsageError, match="not both"):
        resolve_resource_check_flag(True, True, command="annotate")


def test_both_spellings_are_accepted_by_the_cli(tmp_path: Path) -> None:
    """The alias must remain wired, not merely defined.

    Deliberately does NOT use --help to prove it: --help short-circuits
    before the command body, so it happily passes while the body raises
    NameError on a missing import. That exact bug got through an earlier
    version of this test.
    """
    from click.testing import CliRunner

    from karyoscope.cli import main

    help_text = CliRunner().invoke(main, ["annotate", "--help"]).output
    assert "--no-resource-check" in help_text
    assert "--no-space-check" not in help_text  # hidden, but still parses

    fasta = tmp_path / "in.fa"
    fasta.write_text(">c\nACGT\n")
    for flag in ("--no-resource-check", "--no-space-check"):
        result = CliRunner().invoke(
            main,
            [
                "annotate",
                "-i",
                str(fasta),
                "--outdir",
                str(tmp_path / "out"),
                "--db-root",
                str(tmp_path / "nodb"),
                flag,
            ],
        )
        # It must get far enough to fail on the *database*, which proves the
        # flag parsed and the command body ran.
        assert "NameError" not in str(result.exception or ""), flag
        assert result.exit_code != 0, flag


def test_supplying_both_spellings_is_rejected_by_the_cli(tmp_path: Path) -> None:
    from click.testing import CliRunner

    from karyoscope.cli import main

    fasta = tmp_path / "in.fa"
    fasta.write_text(">c\nACGT\n")
    result = CliRunner().invoke(
        main,
        [
            "annotate",
            "-i",
            str(fasta),
            "--outdir",
            str(tmp_path / "out"),
            "--db-root",
            str(tmp_path / "nodb"),
            "--no-resource-check",
            "--no-space-check",
        ],
    )
    assert "not both" in result.output
