"""Unit tests for :mod:`karyoscope.core.annotate` helpers.

The CLI-level end-to-end tests live in ``test_command_annotate.py`` and
need the ``get_featureIDs`` binary; everything here is pure-Python and
runs on the default unit-test pass.
"""

from __future__ import annotations

import json
import os
import threading
from pathlib import Path
from types import SimpleNamespace

from karyoscope.core.annotate import _detect_worker_death, _quiet_worker_pipe_errors
from karyoscope.core.io.kmc import (
    clear_combined_marker,
    combined_bed_is_complete,
    combined_marker_path,
    write_combined_marker,
)

# --- _detect_worker_death (pure function, no threads, no os._exit) -----

# Background: ``multiprocessing.Pool`` does not reassign tasks whose
# worker dies from an external signal -- ``imap_unordered`` just blocks
# forever. The watchdog in :mod:`karyoscope.core.annotate` polls this
# helper to detect the death within a few seconds and fail loudly.
# These tests cover the detection logic without involving real Pool
# objects or threads (the ``os._exit`` path can only be verified by
# manual / subprocess testing -- documented in CHANGELOG).


def _fake_pool(*workers: tuple[int | None, int | None]) -> SimpleNamespace:
    """Build a Pool-shaped stand-in. Each worker is ``(pid, exitcode)``.

    The real :func:`_detect_worker_death` only reads ``pool._pool`` and
    each worker's ``.pid`` / ``.exitcode``, so ``SimpleNamespace`` is
    sufficient.
    """
    return SimpleNamespace(_pool=[SimpleNamespace(pid=pid, exitcode=ec) for pid, ec in workers])


def test_detect_worker_death_returns_none_when_all_alive() -> None:
    pool = _fake_pool((100, None), (101, None), (102, None))
    assert _detect_worker_death(pool, {100, 101, 102}) is None


def test_detect_worker_death_catches_fresh_dead_worker() -> None:
    """A worker with non-None exitcode is reported even before the pool
    replaces it."""
    pool = _fake_pool((100, None), (101, -9), (102, None))
    result = _detect_worker_death(pool, {100, 101, 102})
    assert result == ([101], [-9])


def test_detect_worker_death_catches_pool_replacement() -> None:
    """When the pool's _handle_workers thread has already replaced a
    dead worker, the new PID is evidence of the original's death."""
    # 101 died; pool spawned 200 to replace it.
    pool = _fake_pool((100, None), (200, None), (102, None))
    result = _detect_worker_death(pool, {100, 101, 102})
    # No fresh-dead worker (exitcodes empty); died inferred from
    # missing PID.
    assert result == ([101], [])


def test_detect_worker_death_reports_multiple_deaths() -> None:
    pool = _fake_pool((100, None), (101, -9), (102, -15))
    result = _detect_worker_death(pool, {100, 101, 102})
    assert result == ([101, 102], [-15, -9])


def test_detect_worker_death_handles_simultaneous_death_and_replacement() -> None:
    """Mixed: one dead-not-yet-replaced AND one already-replaced.

    Prefers the precise (fresh-dead) info in that case, since we have
    its exitcode -- the inference fallback only kicks in when there are
    zero fresh-dead workers to report on.
    """
    # 101 died and was replaced by 200; 102 died and not yet replaced.
    pool = _fake_pool((100, None), (200, None), (102, -9))
    result = _detect_worker_death(pool, {100, 101, 102})
    # Reports the worker we actually still have evidence for.
    assert result == ([102], [-9])


def test_detect_worker_death_ignores_zero_exitcode() -> None:
    """A worker with exitcode 0 ought to be unusual during normal
    operation, but is still a failure mode worth surfacing if it
    happens -- the pool hangs the same way regardless of the signal."""
    pool = _fake_pool((100, None), (101, 0))
    result = _detect_worker_death(pool, {100, 101})
    # Reports the pid; exitcode 0 is included in the exitcodes list.
    assert result == ([101], [0])


# --- combined-BED completion marker (pure-Python; no binary) -----------
#
# The marker is what makes "reuse an existing combined BED" safe: it is
# written only after get_featureIDs exits 0, and records the BED's size
# and mtime so a file that was truncated by a killed run (or otherwise
# modified after the fact) is never silently trusted. These tests pin
# the fail-closed behaviour.


def _make_combined(tmp_path: Path) -> Path:
    bed = tmp_path / "x.KS_test.combined.presmoothed.featureIDs.bed"
    bed.write_text("seqA\t0\t10\t1\n")
    return bed


def _write_marker(bed: Path) -> None:
    write_combined_marker(
        bed, prefix="x.KS_test", db_path=Path("/db/features"), input_path=Path("in.fa")
    )


def test_marker_roundtrip_marks_complete(tmp_path: Path) -> None:
    bed = _make_combined(tmp_path)
    assert not combined_bed_is_complete(bed)  # no marker yet
    _write_marker(bed)
    assert combined_bed_is_complete(bed)


def test_missing_bed_is_incomplete(tmp_path: Path) -> None:
    bed = tmp_path / "nope.combined.presmoothed.featureIDs.bed"
    assert not combined_bed_is_complete(bed)


def test_size_change_invalidates_marker(tmp_path: Path) -> None:
    """A combined BED appended to / truncated after the marker was
    written (e.g. a killed run that left a partial file) must not be
    trusted."""
    bed = _make_combined(tmp_path)
    _write_marker(bed)
    assert combined_bed_is_complete(bed)
    with bed.open("a") as fh:
        fh.write("seqB\t0\t5\t2\n")
    assert not combined_bed_is_complete(bed)


def test_mtime_change_invalidates_marker(tmp_path: Path) -> None:
    bed = _make_combined(tmp_path)
    _write_marker(bed)
    st = bed.stat()
    os.utime(bed, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))
    assert not combined_bed_is_complete(bed)


def test_corrupt_marker_is_incomplete(tmp_path: Path) -> None:
    bed = _make_combined(tmp_path)
    _write_marker(bed)
    combined_marker_path(bed).write_text("{ not valid json")
    assert not combined_bed_is_complete(bed)


def test_wrong_schema_is_incomplete(tmp_path: Path) -> None:
    bed = _make_combined(tmp_path)
    _write_marker(bed)
    marker = combined_marker_path(bed)
    data = json.loads(marker.read_text())
    data["schema"] = 999
    marker.write_text(json.dumps(data))
    assert not combined_bed_is_complete(bed)


def test_clear_marker_is_idempotent(tmp_path: Path) -> None:
    bed = _make_combined(tmp_path)
    _write_marker(bed)
    assert combined_marker_path(bed).is_file()
    clear_combined_marker(bed)
    assert not combined_marker_path(bed).is_file()
    clear_combined_marker(bed)  # second call must not raise


# --- _quiet_worker_pipe_errors -----------------------------------------


def _raise_in_thread(exc: BaseException) -> None:
    t = threading.Thread(target=lambda: (_ for _ in ()).throw(exc))
    t.start()
    t.join()


def test_quiet_worker_pipe_errors_suppresses_and_restores() -> None:
    """BrokenPipe/EOF from threads are swallowed inside the context;
    everything else is delegated to the previous hook, which is restored
    on exit."""
    original = threading.excepthook
    seen: list[type] = []

    def record(args: threading.ExceptHookArgs) -> None:
        if args.exc_type is not None:
            seen.append(args.exc_type)

    threading.excepthook = record
    try:
        with _quiet_worker_pipe_errors():
            _raise_in_thread(BrokenPipeError())
            _raise_in_thread(EOFError())
            _raise_in_thread(ValueError("real failure"))
        # Hook restored to whatever it was on entry (our `record`).
        assert threading.excepthook is record
    finally:
        threading.excepthook = original

    # Benign pipe errors suppressed; the real error still reached the
    # delegate.
    assert BrokenPipeError not in seen
    assert EOFError not in seen
    assert ValueError in seen


# --- query-k resolution (variable-k support) -------------------------

import pytest  # noqa: E402

from karyoscope.core.annotate import _resolve_query_k  # noqa: E402
from karyoscope.exceptions import KaryoscopeError  # noqa: E402


def _manifest(size: int, type_: str, max_size: int) -> SimpleNamespace:
    return SimpleNamespace(kmer=SimpleNamespace(size=size, type=type_, max_size=max_size))


def test_resolve_query_k_defaults_to_manifest_size() -> None:
    assert _resolve_query_k(_manifest(31, "fixed", 31), None, "db") == 31


def test_resolve_query_k_same_as_size_always_ok() -> None:
    assert _resolve_query_k(_manifest(31, "fixed", 31), 31, "db") == 31


def test_resolve_query_k_variable_allows_smaller() -> None:
    assert _resolve_query_k(_manifest(31, "variable", 31), 21, "db") == 21


def test_resolve_query_k_fixed_rejects_override() -> None:
    with pytest.raises(KaryoscopeError, match="fixed-k index"):
        _resolve_query_k(_manifest(31, "fixed", 31), 21, "db")


def test_resolve_query_k_variable_rejects_above_max() -> None:
    with pytest.raises(KaryoscopeError, match="exceeds"):
        _resolve_query_k(_manifest(31, "variable", 31), 40, "db")


def test_resolve_query_k_rejects_below_one() -> None:
    with pytest.raises(KaryoscopeError, match=">= 1"):
        _resolve_query_k(_manifest(31, "variable", 31), 0, "db")
