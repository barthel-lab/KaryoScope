"""Tests for ``karyoscope.paths`` (database location resolution)."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.paths import (
    DB_ENV_VAR,
    default_db_root,
    ensure_db_root,
    installed_databases,
)


def test_explicit_argument_wins(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """An explicit path should override the env var and the default."""
    monkeypatch.setenv(DB_ENV_VAR, str(tmp_path / "env"))
    explicit = tmp_path / "explicit"
    assert default_db_root(explicit) == explicit.resolve()


def test_env_var_used_when_no_explicit(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """When no explicit arg is provided, KARYOSCOPE_DB should be honored."""
    target = tmp_path / "from-env"
    monkeypatch.setenv(DB_ENV_VAR, str(target))
    assert default_db_root() == target.resolve()


def test_default_when_neither_set(monkeypatch: pytest.MonkeyPatch) -> None:
    """Falls back to ~/.karyoscope/db when nothing else is set."""
    monkeypatch.delenv(DB_ENV_VAR, raising=False)
    expected = Path.home() / ".karyoscope" / "db"
    assert default_db_root() == expected


def test_ensure_db_root_creates_dir(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """ensure_db_root creates the directory if missing."""
    target = tmp_path / "needs-creation"
    monkeypatch.setenv(DB_ENV_VAR, str(target))
    assert not target.exists()
    created = ensure_db_root()
    assert created == target.resolve()
    assert created.is_dir()


def test_installed_databases_empty_when_root_missing(tmp_path: Path) -> None:
    """A nonexistent root yields an empty list (no exception)."""
    nonexistent = tmp_path / "nope"
    assert installed_databases(nonexistent) == []


def test_installed_databases_finds_dirs_with_manifest(tmp_path: Path) -> None:
    """Only subdirectories that contain a manifest.yaml are reported."""
    # One real-looking database
    db1 = tmp_path / "KS_human_CHM13_v2"
    db1.mkdir()
    (db1 / "manifest.yaml").write_text("id: KS_human_CHM13_v2\n")

    # A subdirectory without a manifest — should be ignored
    (tmp_path / "junk").mkdir()

    # A stray file — should be ignored
    (tmp_path / "README.txt").write_text("not a db")

    found = installed_databases(tmp_path)
    assert found == [db1]
