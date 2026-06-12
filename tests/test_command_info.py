"""Tests for the upgraded ``karyoscope info`` command."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from karyoscope.cli import main


def _invoke(runner: CliRunner, *args: str) -> object:
    return runner.invoke(main, list(args), catch_exceptions=False)


# --- info (no arguments) --------------------------------------------


def test_info_empty_db_root(cli_runner: CliRunner, isolated_db_root: Path) -> None:
    """With nothing installed, `info` should report cleanly."""
    result = _invoke(cli_runner, "info", "--db-root", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert str(isolated_db_root) in result.output
    # Either "No installed databases" or the "root does not exist yet" message.
    assert "No installed" in result.output or "does not exist" in result.output


def test_info_lists_installed_databases(cli_runner: CliRunner, populated_db_root: Path) -> None:
    result = _invoke(cli_runner, "info", "--db-root", str(populated_db_root))
    assert result.exit_code == 0, result.output
    assert "KS_dummy_test_v1" in result.output
    assert "Version:" in result.output
    assert "Size:" in result.output


# --- info <database_id> ---------------------------------------------


def test_info_detailed_view_of_installed_database(
    cli_runner: CliRunner, populated_db_root: Path
) -> None:
    result = _invoke(cli_runner, "info", "KS_dummy_test_v1", "--db-root", str(populated_db_root))
    assert result.exit_code == 0, result.output
    assert "KS_dummy_test_v1" in result.output
    assert "k-mer:" in result.output
    assert "size=21" in result.output
    assert "Feature sets:" in result.output
    assert "chromosome:" in result.output
    assert "region:" in result.output


def test_info_unknown_database_id(cli_runner: CliRunner, isolated_db_root: Path) -> None:
    result = _invoke(cli_runner, "info", "not_a_real_database", "--db-root", str(isolated_db_root))
    assert result.exit_code != 0
    assert "not installed" in result.output


# --- info <path> ----------------------------------------------------


def test_info_existing_database_directory(
    cli_runner: CliRunner,
    unpacked_dummy_db: Path,
    isolated_db_root: Path,
) -> None:
    """Pointing at a database directory on disk works even if not installed."""
    result = _invoke(cli_runner, "info", str(unpacked_dummy_db), "--db-root", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert "KaryoScope database directory" in result.output
    assert "KS_dummy_test_v1" in result.output


def test_info_existing_regular_file(
    cli_runner: CliRunner,
    tmp_path: Path,
    isolated_db_root: Path,
) -> None:
    a_file = tmp_path / "something.bed"
    a_file.write_text("chrom\tstart\tend\tname\n")
    result = _invoke(cli_runner, "info", str(a_file), "--db-root", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert "Type: file" in result.output


def test_info_missing_path(
    cli_runner: CliRunner,
    tmp_path: Path,
    isolated_db_root: Path,
) -> None:
    nope = tmp_path / "does-not-exist"
    result = _invoke(cli_runner, "info", str(nope), "--db-root", str(isolated_db_root))
    assert result.exit_code == 0
    assert "Does not exist" in result.output


# --- --db / --db-root flag handling ---------------------------------


def test_info_db_is_deprecated_alias_for_db_root(
    cli_runner: CliRunner, isolated_db_root: Path
) -> None:
    """The legacy --db flag still resolves the database root."""
    result = _invoke(cli_runner, "info", "--db", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert str(isolated_db_root) in result.output


def test_info_db_and_db_root_conflict(cli_runner: CliRunner, isolated_db_root: Path) -> None:
    result = _invoke(
        cli_runner, "info", "--db", str(isolated_db_root), "--db-root", str(isolated_db_root)
    )
    assert result.exit_code != 0
    assert "not both" in result.output
