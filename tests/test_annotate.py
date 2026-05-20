"""Unit tests for :mod:`karyoscope.core.annotate` helpers.

The CLI-level end-to-end tests live in ``test_command_annotate.py`` and
need the ``get_featureIDs`` binary; everything here is pure-Python and
runs on the default unit-test pass.
"""

from __future__ import annotations

from types import SimpleNamespace

from karyoscope.core.annotate import _detect_worker_death

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
