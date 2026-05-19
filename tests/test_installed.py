"""Tests for :mod:`karyoscope.installed`."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from karyoscope.exceptions import KaryoscopeError
from karyoscope.installed import (
    INSTALLED_FILENAME,
    InstalledRecord,
    load,
    now_iso,
    record_install,
    remove_install_record,
    save,
    uninstall,
)


def _make_record(directory: str = "KS_test_v1") -> InstalledRecord:
    return InstalledRecord(
        version="1.0.0",
        installed_at=now_iso(),
        source_url="https://example.com/x.tar.gz",
        source_sha256="0" * 64,
        directory=directory,
    )


def test_load_missing_file_returns_empty(tmp_path: Path) -> None:
    state = load(tmp_path)
    assert state.databases == {}


def test_round_trip_save_and_load(tmp_path: Path) -> None:
    record = _make_record()
    record_install(tmp_path, "KS_test_v1", record)

    state = load(tmp_path)
    assert "KS_test_v1" in state.databases
    loaded = state.databases["KS_test_v1"]
    assert loaded.version == "1.0.0"
    assert loaded.source_url == "https://example.com/x.tar.gz"


def test_save_is_atomic_format(tmp_path: Path) -> None:
    record_install(tmp_path, "KS_test_v1", _make_record())
    raw = json.loads((tmp_path / INSTALLED_FILENAME).read_text())
    assert raw["schema_version"] == 1
    assert "KS_test_v1" in raw["databases"]


def test_corrupted_file_is_renamed_and_treated_as_empty(tmp_path: Path) -> None:
    (tmp_path / INSTALLED_FILENAME).write_text("this is not valid json {{{")
    state = load(tmp_path)
    assert state.databases == {}
    # The corrupted file should have been preserved with a .corrupt suffix.
    assert (tmp_path / (INSTALLED_FILENAME + ".corrupt")).is_file()


def test_unknown_schema_treated_as_empty(tmp_path: Path) -> None:
    (tmp_path / INSTALLED_FILENAME).write_text(
        json.dumps({"schema_version": 999, "databases": {"x": {}}})
    )
    state = load(tmp_path)
    assert state.databases == {}


def test_replace_existing_record(tmp_path: Path) -> None:
    record_install(tmp_path, "KS_test_v1", _make_record())
    new_record = InstalledRecord(
        version="2.0.0",
        installed_at=now_iso(),
        source_url="https://example.com/y.tar.gz",
        source_sha256="1" * 64,
        directory="KS_test_v1",
    )
    record_install(tmp_path, "KS_test_v1", new_record)

    state = load(tmp_path)
    assert state.databases["KS_test_v1"].version == "2.0.0"


def test_remove_install_record(tmp_path: Path) -> None:
    record_install(tmp_path, "KS_test_v1", _make_record())
    assert remove_install_record(tmp_path, "KS_test_v1") is True
    assert "KS_test_v1" not in load(tmp_path).databases


def test_remove_install_record_returns_false_when_missing(tmp_path: Path) -> None:
    assert remove_install_record(tmp_path, "ghost") is False


def test_uninstall_removes_files_and_record(tmp_path: Path) -> None:
    db_dir = tmp_path / "KS_test_v1"
    db_dir.mkdir()
    (db_dir / "manifest.yaml").write_text("dummy")

    record_install(tmp_path, "KS_test_v1", _make_record(directory="KS_test_v1"))

    assert uninstall(tmp_path, "KS_test_v1") is True
    assert not db_dir.exists()
    assert "KS_test_v1" not in load(tmp_path).databases


def test_uninstall_returns_false_when_not_installed(tmp_path: Path) -> None:
    assert uninstall(tmp_path, "not_there") is False


def test_uninstall_refuses_path_traversal(tmp_path: Path) -> None:
    # Plant a malicious installed.json record that points outside db_root.
    record = InstalledRecord(
        version="1.0.0",
        installed_at=now_iso(),
        source_url="x",
        source_sha256="0" * 64,
        directory="../escape_target",
    )
    record_install(tmp_path, "evil", record)

    with pytest.raises(KaryoscopeError, match="outside"):
        uninstall(tmp_path, "evil")


def test_save_is_atomic_on_crash(tmp_path: Path) -> None:
    """If save() fails mid-write, the existing file should remain intact."""
    record_install(tmp_path, "KS_test_v1", _make_record())
    original = (tmp_path / INSTALLED_FILENAME).read_text()

    # Make save fail by passing an invalid state (we'll patch by hand).
    from karyoscope.installed import InstalledState

    bad_state = InstalledState()
    # Sabotage: make the records non-serializable. asdict() will reject a
    # non-dataclass value with TypeError.
    bad_state.databases["x"] = "not a record"  # type: ignore[assignment]

    with pytest.raises((TypeError, AttributeError)):
        save(tmp_path, bad_state)

    # Existing file untouched.
    assert (tmp_path / INSTALLED_FILENAME).read_text() == original
