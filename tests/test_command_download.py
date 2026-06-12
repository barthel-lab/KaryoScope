"""Tests for the ``karyoscope download`` CLI command."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from karyoscope.cli import main


def _invoke(runner: CliRunner, *args: str) -> object:
    """Helper: invoke the top-level CLI with the given args."""
    return runner.invoke(main, list(args), catch_exceptions=False)


# --- --list -----------------------------------------------------------


def test_list_shows_dummy_db(
    cli_runner: CliRunner, dummy_registry_url: str, isolated_db_root: Path
) -> None:
    result = _invoke(
        cli_runner,
        "download",
        "--list",
        "--registry-url",
        dummy_registry_url,
        "--db-root",
        str(isolated_db_root),
    )
    assert result.exit_code == 0, result.output
    assert "KS_dummy_test_v1" in result.output
    assert "(default)" in result.output


def test_list_organism_filter(
    cli_runner: CliRunner, dummy_registry_url: str, isolated_db_root: Path
) -> None:
    # 'Synthetic' matches; 'human' does not.
    matching = _invoke(
        cli_runner,
        "download",
        "--list",
        "--organism",
        "Synthetic",
        "--registry-url",
        dummy_registry_url,
        "--db-root",
        str(isolated_db_root),
    )
    assert matching.exit_code == 0
    assert "KS_dummy_test_v1" in matching.output

    not_matching = _invoke(
        cli_runner,
        "download",
        "--list",
        "--organism",
        "platypus",
        "--registry-url",
        dummy_registry_url,
        "--db-root",
        str(isolated_db_root),
    )
    assert not_matching.exit_code == 0
    assert "KS_dummy_test_v1" not in not_matching.output
    assert "No databases match" in not_matching.output


def test_list_tag_filter(
    cli_runner: CliRunner, dummy_registry_url: str, isolated_db_root: Path
) -> None:
    result = _invoke(
        cli_runner,
        "download",
        "--list",
        "--tag",
        "test",
        "--registry-url",
        dummy_registry_url,
        "--db-root",
        str(isolated_db_root),
    )
    assert result.exit_code == 0
    assert "KS_dummy_test_v1" in result.output

    result2 = _invoke(
        cli_runner,
        "download",
        "--list",
        "--tag",
        "definitely_not_a_tag",
        "--registry-url",
        dummy_registry_url,
        "--db-root",
        str(isolated_db_root),
    )
    assert result2.exit_code == 0
    assert "KS_dummy_test_v1" not in result2.output


# --- --info -----------------------------------------------------------


def test_info_displays_details(
    cli_runner: CliRunner, dummy_registry_url: str, isolated_db_root: Path
) -> None:
    result = _invoke(
        cli_runner,
        "download",
        "--info",
        "KS_dummy_test_v1",
        "--registry-url",
        dummy_registry_url,
        "--db-root",
        str(isolated_db_root),
    )
    assert result.exit_code == 0, result.output
    assert "KS_dummy_test_v1" in result.output
    assert "Version: 1.0.0" in result.output
    assert "chromosome" in result.output
    assert "SHA-256:" in result.output


def test_info_unknown_id_fails_cleanly(
    cli_runner: CliRunner, dummy_registry_url: str, isolated_db_root: Path
) -> None:
    result = _invoke(
        cli_runner,
        "download",
        "--info",
        "does_not_exist",
        "--registry-url",
        dummy_registry_url,
        "--db-root",
        str(isolated_db_root),
    )
    assert result.exit_code != 0
    assert "not in the registry" in result.output


# --- --status ---------------------------------------------------------


def test_status_empty(cli_runner: CliRunner, isolated_db_root: Path) -> None:
    result = _invoke(
        cli_runner,
        "download",
        "--status",
        "--db-root",
        str(isolated_db_root),
    )
    assert result.exit_code == 0
    assert "No databases installed" in result.output


def test_status_after_install(cli_runner: CliRunner, populated_db_root: Path) -> None:
    result = _invoke(
        cli_runner,
        "download",
        "--status",
        "--db-root",
        str(populated_db_root),
    )
    assert result.exit_code == 0
    assert "KS_dummy_test_v1" in result.output


# --- install (default action) ----------------------------------------


def test_install_default_database(
    cli_runner: CliRunner, dummy_registry_url: str, isolated_db_root: Path
) -> None:
    """Running `karyoscope download` with no args installs the default db."""
    result = _invoke(
        cli_runner,
        "download",
        "--registry-url",
        dummy_registry_url,
        "--db-root",
        str(isolated_db_root),
        "--quiet",
    )
    assert result.exit_code == 0, result.output
    assert "Installing KS_dummy_test_v1" in result.output
    assert (isolated_db_root / "KS_dummy_test_v1" / "manifest.yaml").is_file()


def test_install_by_id(
    cli_runner: CliRunner, dummy_registry_url: str, isolated_db_root: Path
) -> None:
    result = _invoke(
        cli_runner,
        "download",
        "KS_dummy_test_v1",
        "--registry-url",
        dummy_registry_url,
        "--db-root",
        str(isolated_db_root),
        "--quiet",
    )
    assert result.exit_code == 0, result.output
    assert (isolated_db_root / "KS_dummy_test_v1").is_dir()


def test_install_unknown_id_fails(
    cli_runner: CliRunner, dummy_registry_url: str, isolated_db_root: Path
) -> None:
    result = _invoke(
        cli_runner,
        "download",
        "ghost_db",
        "--registry-url",
        dummy_registry_url,
        "--db-root",
        str(isolated_db_root),
        "--quiet",
    )
    assert result.exit_code != 0
    assert "not in the registry" in result.output


def test_install_skips_if_already_present(
    cli_runner: CliRunner, dummy_registry_url: str, isolated_db_root: Path
) -> None:
    # First install
    _invoke(
        cli_runner,
        "download",
        "KS_dummy_test_v1",
        "--registry-url",
        dummy_registry_url,
        "--db-root",
        str(isolated_db_root),
        "--quiet",
    )
    # Second install: should skip
    result = _invoke(
        cli_runner,
        "download",
        "KS_dummy_test_v1",
        "--registry-url",
        dummy_registry_url,
        "--db-root",
        str(isolated_db_root),
        "--quiet",
    )
    assert result.exit_code == 0
    assert "already installed" in result.output


# --- --remove ---------------------------------------------------------


def test_remove_with_yes_flag(cli_runner: CliRunner, populated_db_root: Path) -> None:
    result = _invoke(
        cli_runner,
        "download",
        "--remove",
        "KS_dummy_test_v1",
        "-y",
        "--db-root",
        str(populated_db_root),
    )
    assert result.exit_code == 0, result.output
    assert "Removed KS_dummy_test_v1" in result.output
    assert not (populated_db_root / "KS_dummy_test_v1").exists()


def test_remove_unknown_id(cli_runner: CliRunner, isolated_db_root: Path) -> None:
    result = _invoke(
        cli_runner,
        "download",
        "--remove",
        "ghost",
        "-y",
        "--db-root",
        str(isolated_db_root),
    )
    assert result.exit_code != 0
    assert "not installed" in result.output


# --- mutually exclusive flags -----------------------------------------


def test_action_flags_are_mutually_exclusive(
    cli_runner: CliRunner, dummy_registry_url: str, isolated_db_root: Path
) -> None:
    result = cli_runner.invoke(
        main,
        [
            "download",
            "--list",
            "--status",
            "--registry-url",
            dummy_registry_url,
            "--db-root",
            str(isolated_db_root),
        ],
    )
    assert result.exit_code != 0
    assert "Only one of" in result.output


# --- --db / --db-root flag handling -----------------------------------


def test_db_is_deprecated_alias_for_db_root(cli_runner: CliRunner, isolated_db_root: Path) -> None:
    """The legacy --db flag still resolves the database root (here via --status)."""
    result = _invoke(cli_runner, "download", "--status", "--db", str(isolated_db_root))
    assert result.exit_code == 0, result.output


def test_db_and_db_root_conflict(cli_runner: CliRunner, isolated_db_root: Path) -> None:
    result = _invoke(
        cli_runner,
        "download",
        "--status",
        "--db",
        str(isolated_db_root),
        "--db-root",
        str(isolated_db_root),
    )
    assert result.exit_code != 0
    assert "not both" in result.output
