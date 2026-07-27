"""Tests for :mod:`karyoscope.diskspace`.

The behaviour under test is arithmetic and message formatting, so free
space is stubbed rather than actually exhausted -- filling a real
filesystem in a test is neither portable nor kind to the machine running
it.
"""

from __future__ import annotations

import errno
from pathlib import Path

import pytest

from karyoscope import diskspace
from karyoscope.exceptions import InsufficientDiskSpaceError


@pytest.fixture
def stub_free(monkeypatch: pytest.MonkeyPatch):
    """Return a setter that pins :func:`diskspace.free_bytes` to a value."""

    def _set(n: int) -> None:
        monkeypatch.setattr(diskspace, "free_bytes", lambda _path: n)

    return _set


# --- format_bytes -----------------------------------------------------


@pytest.mark.parametrize(
    ("n", "expected"),
    [
        (0, "0 B"),
        (512, "512 B"),
        (1_500, "1.5 kB"),
        (2_500_000, "2.5 MB"),
        (13_337_554_753, "13.3 GB"),  # the real HKS archive
        (22_695_262_463, "22.7 GB"),  # the real extracted HKS database
        (2_500 * diskspace.GB, "2.5 TB"),
    ],
)
def test_format_bytes(n: int, expected: str) -> None:
    assert diskspace.format_bytes(n) == expected


# --- free_bytes / _nearest_existing -----------------------------------


def test_free_bytes_works_for_a_path_that_does_not_exist_yet(tmp_path: Path) -> None:
    """An --outdir is routinely checked before it's created."""
    missing = tmp_path / "not" / "created" / "yet"
    assert diskspace.free_bytes(missing) > 0


def test_nearest_existing_walks_up_to_a_real_directory(tmp_path: Path) -> None:
    assert diskspace._nearest_existing(tmp_path / "a" / "b" / "c") == tmp_path.resolve()


# --- directory_size ---------------------------------------------------


def test_directory_size_sums_files_recursively(tmp_path: Path) -> None:
    (tmp_path / "sub").mkdir()
    (tmp_path / "a.bin").write_bytes(b"x" * 100)
    (tmp_path / "sub" / "b.bin").write_bytes(b"y" * 250)
    assert diskspace.directory_size(tmp_path) == 350


def test_directory_size_of_a_missing_directory_is_zero(tmp_path: Path) -> None:
    assert diskspace.directory_size(tmp_path / "nope") == 0


# --- require_free_space -----------------------------------------------


def test_require_free_space_passes_when_there_is_room(tmp_path: Path, stub_free) -> None:
    stub_free(100 * diskspace.GB)
    diskspace.require_free_space(tmp_path, 10 * diskspace.GB, what="testing")


def test_require_free_space_raises_when_short(tmp_path: Path, stub_free) -> None:
    stub_free(12 * diskspace.GB)
    with pytest.raises(InsufficientDiskSpaceError) as excinfo:
        diskspace.require_free_space(tmp_path, 36 * diskspace.GB, what="installing a database")
    message = str(excinfo.value)
    assert "installing a database" in message
    assert "Short by:" in message
    assert str(tmp_path.resolve()) in message


def test_margin_is_applied_on_top_of_the_requirement(tmp_path: Path, stub_free) -> None:
    """A request that exactly matches free space must still fail.

    Driving a filesystem to precisely zero free bytes is how you get an
    ENOSPC on the very last write, which is the failure this check exists
    to avoid -- so "just barely fits" has to count as not fitting.
    """
    stub_free(10 * diskspace.GB)
    with pytest.raises(InsufficientDiskSpaceError):
        diskspace.require_free_space(
            tmp_path, 10 * diskspace.GB, what="testing", margin=diskspace.DEFAULT_MARGIN
        )
    # ...but with no margin requested, it is allowed through.
    diskspace.require_free_space(tmp_path, 10 * diskspace.GB, what="testing", margin=0.0)


def test_credit_bytes_counts_toward_available(tmp_path: Path, stub_free) -> None:
    """Space the operation frees first is space the operation can use.

    ``download --force`` removes the old install before fetching the new
    one, so a same-size reinstall must not be rejected on a full disk.
    """
    stub_free(1 * diskspace.GB)
    with pytest.raises(InsufficientDiskSpaceError):
        diskspace.require_free_space(tmp_path, 20 * diskspace.GB, what="reinstalling")
    diskspace.require_free_space(
        tmp_path, 20 * diskspace.GB, what="reinstalling", credit_bytes=25 * diskspace.GB
    )


def test_skip_bypasses_the_check(tmp_path: Path, stub_free) -> None:
    """--no-space-check must let a run through even when clearly short."""
    stub_free(0)
    diskspace.require_free_space(tmp_path, 100 * diskspace.GB, what="testing", skip=True)


def test_estimated_wording_appears_only_when_requested(tmp_path: Path, stub_free) -> None:
    stub_free(0)
    with pytest.raises(InsufficientDiskSpaceError) as measured:
        diskspace.require_free_space(tmp_path, diskspace.GB, what="testing")
    with pytest.raises(InsufficientDiskSpaceError) as guessed:
        diskspace.require_free_space(tmp_path, diskspace.GB, what="testing", estimated=True)
    assert "an estimated" not in str(measured.value)
    assert "an estimated" in str(guessed.value)


def test_hint_is_included_in_the_message(tmp_path: Path, stub_free) -> None:
    stub_free(0)
    with pytest.raises(InsufficientDiskSpaceError) as excinfo:
        diskspace.require_free_space(
            tmp_path, diskspace.GB, what="testing", hint="Try --outdir elsewhere."
        )
    assert "Try --outdir elsewhere." in str(excinfo.value)


# --- ENOSPC translation -----------------------------------------------


def test_is_enospc_recognises_disk_full_and_quota() -> None:
    assert diskspace.is_enospc(OSError(errno.ENOSPC, "No space left on device"))
    assert diskspace.is_enospc(OSError(errno.EDQUOT, "Disk quota exceeded"))
    assert not diskspace.is_enospc(OSError(errno.ENOENT, "No such file"))
    assert not diskspace.is_enospc(ValueError("unrelated"))


def test_enospc_error_names_the_file_that_failed(tmp_path: Path) -> None:
    target = tmp_path / "output.bed"
    exc = OSError(errno.ENOSPC, "No space left on device", str(target))
    message = str(diskspace.enospc_error(exc, what="annotating"))
    assert "annotating" in message
    assert str(target) in message


def test_reframe_enospc_passes_other_errors_through_untouched() -> None:
    """Only disk-full errors get reworded; everything else must not be masked."""
    original = OSError(errno.EACCES, "Permission denied")
    assert diskspace.reframe_enospc(original, what="testing") is original


def test_reframe_enospc_converts_disk_full() -> None:
    original = OSError(errno.ENOSPC, "No space left on device")
    reframed = diskspace.reframe_enospc(original, what="testing")
    assert isinstance(reframed, InsufficientDiskSpaceError)
