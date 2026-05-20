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


def _format_string() -> str:
    """The active root-logger format string (assumes single handler)."""
    root = logging.getLogger()
    assert root.handlers, "no logging handler installed"
    fmt = root.handlers[0].formatter
    assert fmt is not None
    return fmt._fmt  # type: ignore[attr-defined]


def test_default_format_has_no_timestamp(cli_runner: CliRunner) -> None:
    # The default-verbosity user mostly sees one or two log lines;
    # a timestamp would be noise. The format should be the minimal
    # LEVEL: message form.
    cli_runner.invoke(main, ["version"], catch_exceptions=False)
    fmt = _format_string()
    assert "%(asctime)s" not in fmt


def test_dash_v_format_includes_timestamp(cli_runner: CliRunner) -> None:
    # At -v (INFO) the per-step log lines are useful for extracting
    # wall times via diffing consecutive timestamps. The format
    # should include %(asctime)s.
    cli_runner.invoke(main, ["-v", "version"], catch_exceptions=False)
    fmt = _format_string()
    assert "%(asctime)s" in fmt
    assert "%(name)s" not in fmt


def test_dash_vv_format_includes_timestamp_and_module(cli_runner: CliRunner) -> None:
    # At -vv (DEBUG) we also include %(name)s so subsystem-level
    # debugging is easier to trace.
    cli_runner.invoke(main, ["-vv", "version"], catch_exceptions=False)
    fmt = _format_string()
    assert "%(asctime)s" in fmt
    assert "%(name)s" in fmt
