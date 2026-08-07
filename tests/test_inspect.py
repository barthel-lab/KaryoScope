"""Unit tests for :mod:`karyoscope.inspect`.

The archive-staging logic used to live inside the ``info`` command and
was only testable through click; these test it directly. The command-
level behavior stays covered by ``test_command_info.py``.
"""

from __future__ import annotations

import io
import tarfile
from pathlib import Path, PurePosixPath

from karyoscope.inspect import dir_size, safe_relpath, stage_archive, tar_member_issue

# --- tar_member_issue -------------------------------------------------


def _member(name: str, *, kind: int = tarfile.REGTYPE) -> tarfile.TarInfo:
    m = tarfile.TarInfo(name)
    m.type = kind
    return m


def test_regular_files_and_dirs_pass() -> None:
    assert tar_member_issue(_member("db/manifest.yaml")) is None
    assert tar_member_issue(_member("db/index", kind=tarfile.DIRTYPE)) is None


def test_special_members_are_named() -> None:
    for kind in (tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.FIFOTYPE, tarfile.CHRTYPE):
        issue = tar_member_issue(_member("db/evil", kind=kind))
        assert issue is not None and "special files" in issue


def test_traversal_is_named() -> None:
    issue = tar_member_issue(_member("db/../../etc/passwd"))
    assert issue is not None and "'..'" in issue
    # An absolute path that climbs after the root is still caught.
    assert tar_member_issue(_member("/db/../x")) is not None


# --- safe_relpath ------------------------------------------------------


def test_safe_relpath_normalises_and_rejects() -> None:
    assert safe_relpath("db/./manifest.yaml") == PurePosixPath("db/manifest.yaml")
    assert safe_relpath("/db/manifest.yaml") is None
    assert safe_relpath("db/../x") is None
    assert safe_relpath(".") is None


# --- dir_size ----------------------------------------------------------


def test_dir_size_sums_files_recursively(tmp_path: Path) -> None:
    (tmp_path / "a").write_bytes(b"x" * 10)
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "b").write_bytes(b"x" * 32)
    assert dir_size(tmp_path) == 42


# --- stage_archive -----------------------------------------------------


def _make_archive(path: Path, members: list[tuple[str, bytes | None]]) -> None:
    """Write a tar.gz; ``None`` data means a directory member."""
    with tarfile.open(path, "w:gz") as tar:
        for name, data in members:
            if data is None:
                info = tarfile.TarInfo(name)
                info.type = tarfile.DIRTYPE
                tar.addfile(info)
            else:
                info = tarfile.TarInfo(name)
                info.size = len(data)
                tar.addfile(info, io.BytesIO(data))


def test_stage_archive_inlines_text_and_placeholders_the_rest(tmp_path: Path) -> None:
    archive = tmp_path / "db.tar.gz"
    _make_archive(
        archive,
        [
            ("db", None),
            ("db/manifest.yaml", b"id: X\n"),
            ("db/index", None),
            ("db/index/big.hksb", b"\x00" * 1024),
        ],
    )
    staging = tmp_path / "staging"
    staging.mkdir()
    top_level, n_files, total = stage_archive(archive, staging)

    assert top_level == {"db"}
    assert n_files == 2
    assert total == len(b"id: X\n") + 1024
    # The text member arrives in full; the index is an empty placeholder
    # standing in for gigabytes.
    assert (staging / "db/manifest.yaml").read_text() == "id: X\n"
    assert (staging / "db/index/big.hksb").stat().st_size == 0


def test_stage_archive_skips_unsafe_members(tmp_path: Path) -> None:
    archive = tmp_path / "evil.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        info = tarfile.TarInfo("db/ok.txt")
        info.size = 2
        tar.addfile(info, io.BytesIO(b"ok"))
        info = tarfile.TarInfo("../escape.txt")
        info.size = 2
        tar.addfile(info, io.BytesIO(b"no"))
        link = tarfile.TarInfo("db/link")
        link.type = tarfile.SYMTYPE
        link.linkname = "/etc/passwd"
        tar.addfile(link)
    staging = tmp_path / "staging"
    staging.mkdir()
    top_level, n_files, _ = stage_archive(archive, staging)

    assert top_level == {"db"}
    assert n_files == 1
    assert not (tmp_path / "escape.txt").exists()
    assert set(staging.rglob("*")) == {staging / "db", staging / "db" / "ok.txt"}
