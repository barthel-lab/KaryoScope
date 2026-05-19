"""Tests for :mod:`karyoscope._fetch`."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from karyoscope._fetch import fetch
from karyoscope.exceptions import ChecksumError, FetchError, UnsupportedSchemeError


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def test_fetch_file_url_copies_content(tmp_path: Path) -> None:
    source = tmp_path / "source.bin"
    source.write_bytes(b"hello world")
    dest = tmp_path / "dest.bin"

    fetch(source.absolute().as_uri(), dest, show_progress=False)

    assert dest.is_file()
    assert dest.read_bytes() == b"hello world"


def test_fetch_file_url_verifies_checksum(tmp_path: Path) -> None:
    source = tmp_path / "source.bin"
    payload = b"some test content"
    source.write_bytes(payload)
    dest = tmp_path / "dest.bin"

    fetch(
        source.absolute().as_uri(),
        dest,
        expected_sha256=_sha256(payload),
        show_progress=False,
    )

    assert dest.read_bytes() == payload


def test_fetch_file_url_raises_on_bad_checksum(tmp_path: Path) -> None:
    source = tmp_path / "source.bin"
    source.write_bytes(b"hello")
    dest = tmp_path / "dest.bin"
    wrong_sha = "0" * 64

    with pytest.raises(ChecksumError):
        fetch(
            source.absolute().as_uri(),
            dest,
            expected_sha256=wrong_sha,
            show_progress=False,
        )


def test_fetch_file_url_creates_parent_dir(tmp_path: Path) -> None:
    source = tmp_path / "source.bin"
    source.write_bytes(b"x")
    dest = tmp_path / "deep" / "nested" / "dest.bin"

    fetch(source.absolute().as_uri(), dest, show_progress=False)

    assert dest.is_file()


def test_fetch_unsupported_scheme_raises(tmp_path: Path) -> None:
    with pytest.raises(UnsupportedSchemeError):
        fetch(
            "ftp://example.com/file.tar.gz",
            tmp_path / "dest.bin",
            show_progress=False,
        )


def test_fetch_file_url_missing_source_raises(tmp_path: Path) -> None:
    missing = tmp_path / "does_not_exist.bin"
    with pytest.raises(FetchError):
        fetch(
            missing.absolute().as_uri(),
            tmp_path / "dest.bin",
            show_progress=False,
        )


def test_fetch_checksum_is_case_insensitive(tmp_path: Path) -> None:
    source = tmp_path / "source.bin"
    payload = b"case-insensitive checksums"
    source.write_bytes(payload)
    dest = tmp_path / "dest.bin"

    # Upper-case digest should still match.
    fetch(
        source.absolute().as_uri(),
        dest,
        expected_sha256=_sha256(payload).upper(),
        show_progress=False,
    )
    assert dest.read_bytes() == payload
