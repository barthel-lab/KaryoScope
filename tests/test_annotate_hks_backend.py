"""Tests for :func:`karyoscope.core.annotate._run_hks_backend`.

The backend runs the presmoothed convert and ``hks smooth`` concurrently
against the same raw lookup TSV, then deletes that TSV. These tests pin the
properties that concurrency puts at risk: that both halves actually run, that
neither is cut short by the other's failure, and that the shared input
survives until both are done with it.
"""

from __future__ import annotations

import threading
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from karyoscope.core import annotate as annotate_mod
from karyoscope.exceptions import KaryoscopeError


class _Manifest:
    def __init__(self) -> None:
        self.index = SimpleNamespace(basename="features", type="hks")


def _backend_kwargs(tmp_path: Path, **overrides: object) -> dict[str, object]:
    db_dir = tmp_path / "db"
    db_dir.mkdir(exist_ok=True)
    out_dir = tmp_path / "out"
    out_dir.mkdir(exist_ok=True)
    kwargs: dict[str, object] = {
        "manifest": _Manifest(),
        "db_dir": db_dir,
        "input_path": tmp_path / "sample.fa",
        "prefix": "sample.DB",
        "output_dir": out_dir,
        "requested": ["chromosome"],
        "smooth": True,
        "keep_presmoothed": True,
        "presmoothed_paths": {"chromosome": out_dir / "sample.chromosome.presmoothed.bed"},
        "smoothed_paths": {"chromosome": out_dir / "sample.chromosome.smoothed.bed"},
        "threads": 4,
        "k": 31,
    }
    kwargs.update(overrides)
    return kwargs


@pytest.fixture
def fake_lookup(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make ``run_hks_lookup`` write a plausible raw TSV instead of running hks."""

    def _lookup(*, output_path: Path, **_: object) -> Path:
        output_path.write_text("query_name\tfrom_kmer\tto_kmer\tlabel_name\nchr1\t0\t10\tchr1\n")
        return output_path

    monkeypatch.setattr(annotate_mod, "run_hks_lookup", _lookup)


def test_convert_and_smooth_both_run_and_tsv_is_cleaned(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fake_lookup: None
) -> None:
    seen: dict[str, Path] = {}

    def _convert(tsv: Path, bed: Path, **_: object) -> None:
        seen["convert"] = tsv
        bed.write_text("chr1\t0\t10\tchr1\n")

    def _smooth(*, input_path: Path, output_path: Path, **_: object) -> Path:
        seen["smooth"] = input_path
        output_path.write_text("chr1\t0\t10\tchr1\n")
        return output_path

    monkeypatch.setattr(annotate_mod, "convert_hks_tsv_to_bed", _convert)
    monkeypatch.setattr(annotate_mod, "run_hks_smooth", _smooth)

    kwargs = _backend_kwargs(tmp_path)
    annotate_mod._run_hks_backend(**kwargs)  # type: ignore[arg-type]

    assert seen["convert"] == seen["smooth"], "both halves read the same raw TSV"
    assert kwargs["presmoothed_paths"]["chromosome"].is_file()  # type: ignore[index]
    assert kwargs["smoothed_paths"]["chromosome"].is_file()  # type: ignore[index]
    assert not seen["convert"].exists(), "raw TSV removed once both halves finished"


def test_halves_actually_overlap(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fake_lookup: None
) -> None:
    """Each half must observe the other in flight, or they are still serial."""
    both_started = threading.Barrier(2, timeout=10)
    overlapped = {"convert": False, "smooth": False}

    def _convert(tsv: Path, bed: Path, **_: object) -> None:
        both_started.wait()
        overlapped["convert"] = True
        bed.write_text("x\n")

    def _smooth(*, input_path: Path, output_path: Path, **_: object) -> Path:
        both_started.wait()
        overlapped["smooth"] = True
        output_path.write_text("x\n")
        return output_path

    monkeypatch.setattr(annotate_mod, "convert_hks_tsv_to_bed", _convert)
    monkeypatch.setattr(annotate_mod, "run_hks_smooth", _smooth)

    # A serial implementation deadlocks on the barrier and fails via timeout.
    annotate_mod._run_hks_backend(**_backend_kwargs(tmp_path))  # type: ignore[arg-type]
    assert overlapped == {"convert": True, "smooth": True}


def test_smooth_failure_lets_convert_finish_and_propagates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fake_lookup: None
) -> None:
    """A crash in one half must not yank the TSV from the other mid-read."""
    convert_done = threading.Event()

    def _convert(tsv: Path, bed: Path, **_: object) -> None:
        time.sleep(0.2)
        # Would raise if the failing sibling had already unlinked the input.
        tsv.read_text()
        bed.write_text("chr1\t0\t10\tchr1\n")
        convert_done.set()

    def _smooth(**_: object) -> Path:
        raise KaryoscopeError("hks smooth exploded")

    monkeypatch.setattr(annotate_mod, "convert_hks_tsv_to_bed", _convert)
    monkeypatch.setattr(annotate_mod, "run_hks_smooth", _smooth)

    with pytest.raises(KaryoscopeError, match="hks smooth exploded"):
        annotate_mod._run_hks_backend(**_backend_kwargs(tmp_path))  # type: ignore[arg-type]

    assert convert_done.is_set(), "convert was joined, not abandoned"


def test_smooth_only_still_works(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fake_lookup: None
) -> None:
    """--no-keep-presmoothed leaves a single half; the pool must handle that."""

    def _smooth(*, input_path: Path, output_path: Path, **_: object) -> Path:
        output_path.write_text("chr1\t0\t10\tchr1\n")
        return output_path

    monkeypatch.setattr(annotate_mod, "run_hks_smooth", _smooth)
    monkeypatch.setattr(
        annotate_mod,
        "convert_hks_tsv_to_bed",
        lambda *a, **k: pytest.fail("convert must not run when keep_presmoothed is False"),
    )

    kwargs = _backend_kwargs(tmp_path, keep_presmoothed=False, presmoothed_paths={})
    annotate_mod._run_hks_backend(**kwargs)  # type: ignore[arg-type]
    assert kwargs["smoothed_paths"]["chromosome"].is_file()  # type: ignore[index]
