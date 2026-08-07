"""Tests for the ``karyoscope register`` command."""

from __future__ import annotations

import re
import tarfile
from pathlib import Path

import pytest
from click.testing import CliRunner

from karyoscope import installed as _installed
from karyoscope.cli import main
from karyoscope.exceptions import DatabaseNotFoundError
from karyoscope.installed import resolve_database

from .conftest import DUMMY_DB_ID, _extractall_compat


def _invoke(runner: CliRunner, *args: str) -> object:
    return runner.invoke(main, list(args), catch_exceptions=False)


@pytest.fixture
def unregistered_db_root(tmp_path: Path, dummy_db_tarball: Path) -> Path:
    """A database root with the dummy db extracted but NOT recorded in installed.json."""
    db_root = tmp_path / "karyoscope_db"
    db_root.mkdir()
    with tarfile.open(dummy_db_tarball, "r:gz") as tar:
        _extractall_compat(tar, db_root)
    assert (db_root / DUMMY_DB_ID).is_dir()
    # Sanity: nothing recorded yet.
    assert not _installed.load(db_root).databases
    return db_root


def test_register_by_id(cli_runner: CliRunner, unregistered_db_root: Path) -> None:
    result = _invoke(cli_runner, "register", DUMMY_DB_ID, "--db-root", str(unregistered_db_root))
    assert result.exit_code == 0, result.output
    assert f"Registered {DUMMY_DB_ID}" in result.output

    state = _installed.load(unregistered_db_root)
    assert DUMMY_DB_ID in state.databases
    rec = state.databases[DUMMY_DB_ID]
    assert rec.version == "1.0.0"  # from the dummy manifest
    assert rec.directory == DUMMY_DB_ID
    assert rec.source_url == "local"
    assert rec.source_sha256 == ""


def test_register_by_path(cli_runner: CliRunner, unregistered_db_root: Path) -> None:
    db_dir = unregistered_db_root / DUMMY_DB_ID
    result = _invoke(cli_runner, "register", str(db_dir), "--db-root", str(unregistered_db_root))
    assert result.exit_code == 0, result.output
    assert DUMMY_DB_ID in _installed.load(unregistered_db_root).databases


def test_register_installed_at_is_registration_time(
    cli_runner: CliRunner, unregistered_db_root: Path
) -> None:
    """installed_at records when register ran, in the same format as now_iso."""
    result = _invoke(cli_runner, "register", DUMMY_DB_ID, "--db-root", str(unregistered_db_root))
    assert result.exit_code == 0, result.output
    installed_at = _installed.load(unregistered_db_root).databases[DUMMY_DB_ID].installed_at
    # Well-formed ISO-8601 UTC timestamp, and a real one (not epoch).
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z", installed_at)
    assert installed_at >= "2025-01-01T00:00:00Z"


def test_register_makes_database_resolvable(
    cli_runner: CliRunner, unregistered_db_root: Path
) -> None:
    """After registering, the data-command resolver must find the database."""
    # Before: resolution fails because it isn't recorded.
    with pytest.raises(DatabaseNotFoundError):
        resolve_database(unregistered_db_root, DUMMY_DB_ID)

    result = _invoke(cli_runner, "register", DUMMY_DB_ID, "--db-root", str(unregistered_db_root))
    assert result.exit_code == 0, result.output

    db_id, db_dir = resolve_database(unregistered_db_root, DUMMY_DB_ID)
    assert db_id == DUMMY_DB_ID
    assert db_dir == (unregistered_db_root / DUMMY_DB_ID).resolve()


def test_register_already_registered_errors_without_force(
    cli_runner: CliRunner, unregistered_db_root: Path
) -> None:
    first = _invoke(cli_runner, "register", DUMMY_DB_ID, "--db-root", str(unregistered_db_root))
    assert first.exit_code == 0, first.output

    second = _invoke(cli_runner, "register", DUMMY_DB_ID, "--db-root", str(unregistered_db_root))
    assert second.exit_code != 0
    assert "already registered" in second.output
    assert "--force" in second.output


def test_register_force_overwrites(cli_runner: CliRunner, unregistered_db_root: Path) -> None:
    _invoke(cli_runner, "register", DUMMY_DB_ID, "--db-root", str(unregistered_db_root))
    result = _invoke(
        cli_runner, "register", DUMMY_DB_ID, "--db-root", str(unregistered_db_root), "--force"
    )
    assert result.exit_code == 0, result.output
    assert "Re-registered" in result.output


def test_register_unknown_id_errors(cli_runner: CliRunner, isolated_db_root: Path) -> None:
    isolated_db_root.mkdir(parents=True, exist_ok=True)
    result = _invoke(cli_runner, "register", "not_a_real_db", "--db-root", str(isolated_db_root))
    assert result.exit_code != 0
    assert "no database directory found" in result.output


def test_register_directory_outside_db_root_errors(
    cli_runner: CliRunner, unpacked_dummy_db: Path, isolated_db_root: Path
) -> None:
    """A database directory outside the db root cannot be registered."""
    isolated_db_root.mkdir(parents=True, exist_ok=True)
    result = _invoke(
        cli_runner, "register", str(unpacked_dummy_db), "--db-root", str(isolated_db_root)
    )
    assert result.exit_code != 0
    assert "not inside the database root" in result.output
