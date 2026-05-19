"""Tests for the top-level verbosity flags (``-v``/``-vv``/``-q``).

These verify the CLI parses the flags correctly and that the resulting
log level is what we expect. Actual log emission from subcommands is
tested incidentally by integration tests in :mod:`test_command_download`
etc.
"""

from __future__ import annotations

import logging

from click.testing import CliRunner

from karyoscope.cli import main


def _root_level() -> int:
    return logging.getLogger().level


def test_default_verbosity_is_warning(cli_runner: CliRunner) -> None:
    cli_runner.invoke(main, ["version"], catch_exceptions=False)
    assert _root_level() == logging.WARNING


def test_dash_v_sets_info(cli_runner: CliRunner) -> None:
    cli_runner.invoke(main, ["-v", "version"], catch_exceptions=False)
    assert _root_level() == logging.INFO


def test_dash_vv_sets_debug(cli_runner: CliRunner) -> None:
    cli_runner.invoke(main, ["-vv", "version"], catch_exceptions=False)
    assert _root_level() == logging.DEBUG


def test_more_v_clamps_to_debug(cli_runner: CliRunner) -> None:
    """Beyond -vv we just stay at DEBUG, not lower."""
    cli_runner.invoke(main, ["-vvvv", "version"], catch_exceptions=False)
    assert _root_level() == logging.DEBUG


def test_dash_q_sets_error(cli_runner: CliRunner) -> None:
    cli_runner.invoke(main, ["-q", "version"], catch_exceptions=False)
    assert _root_level() == logging.ERROR


def test_q_and_v_conflict(cli_runner: CliRunner) -> None:
    result = cli_runner.invoke(main, ["-q", "-v", "version"])
    assert result.exit_code != 0
    assert "cannot be combined" in result.output


def test_verbosity_does_not_break_help(cli_runner: CliRunner) -> None:
    result = cli_runner.invoke(main, ["-v", "--help"])
    assert result.exit_code == 0
    # CliRunner reports the command name as the function name ('main'), so
    # we look for stable content from the docstring rather than the prog
    # name itself.
    assert "rapid, alignment-free sequence annotation" in result.output
