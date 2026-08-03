"""Tests for :mod:`karyoscope.cpus`.

The bug these pin: ``--threads 0`` sized its worker pool from
``os.cpu_count()``, which reports the *machine's* cores and ignores every
mechanism that restricts a process to a subset of them. On a shared 36-core
SLURM node with one allocated CPU, the default spawned 36 workers.
"""

from __future__ import annotations

import logging
import os

import pytest

from karyoscope import cpus


@pytest.fixture(autouse=True)
def _clean_allocation_env(monkeypatch: pytest.MonkeyPatch):
    """Run with no scheduler variables unless a test sets them.

    Otherwise the suite's results depend on whether it was launched inside
    a SLURM job, which is exactly the confusion this module exists to fix.
    """
    for var in cpus._ALLOCATION_VARS:
        monkeypatch.delenv(var, raising=False)


# --- allocation takes precedence --------------------------------------


def test_slurm_cpus_per_task_wins_over_the_machine(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "4")
    monkeypatch.setattr(os, "cpu_count", lambda: 128)
    assert cpus.usable_cpus() == 4


def test_cpus_per_task_is_preferred_over_job_cpus_per_node(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The per-task grant is what this process may use."""
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "4")
    monkeypatch.setenv("SLURM_JOB_CPUS_PER_NODE", "36")
    assert cpus.usable_cpus() == 4


def test_job_cpus_per_node_range_syntax_is_tolerated(monkeypatch: pytest.MonkeyPatch) -> None:
    """SLURM writes things like '4(x2)'; take the leading integer."""
    monkeypatch.setenv("SLURM_JOB_CPUS_PER_NODE", "4(x2)")
    assert cpus.usable_cpus() == 4


@pytest.mark.parametrize("junk", ["", "none", "0", "x4"])
def test_unusable_allocation_values_fall_through(
    junk: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A malformed variable must not yield 0 workers or crash."""
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", junk)
    assert cpus.usable_cpus() >= 1


# --- affinity ---------------------------------------------------------


@pytest.mark.skipif(not hasattr(os, "sched_getaffinity"), reason="sched_getaffinity is Linux-only")
def test_affinity_is_used_when_no_allocation_is_set(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(os, "sched_getaffinity", lambda _pid: {0, 1, 2})
    monkeypatch.setattr(os, "cpu_count", lambda: 128)
    assert cpus.usable_cpus() == 3


def test_falls_back_to_cpu_count_without_affinity(monkeypatch: pytest.MonkeyPatch) -> None:
    """macOS exposes no affinity API, so the machine count is all there is."""
    monkeypatch.delattr(os, "sched_getaffinity", raising=False)
    monkeypatch.delattr(os, "process_cpu_count", raising=False)
    monkeypatch.setattr(os, "cpu_count", lambda: 10)
    assert cpus.usable_cpus() == 10


def test_never_returns_zero(monkeypatch: pytest.MonkeyPatch) -> None:
    """Callers size a worker pool with this directly."""
    monkeypatch.delattr(os, "sched_getaffinity", raising=False)
    monkeypatch.delattr(os, "process_cpu_count", raising=False)
    monkeypatch.setattr(os, "cpu_count", lambda: None)
    assert cpus.usable_cpus() == 1


# --- resolve_threads --------------------------------------------------


def test_auto_resolves_to_usable_not_machine_cores(monkeypatch: pytest.MonkeyPatch) -> None:
    """The reported bug, in one assertion."""
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
    monkeypatch.setattr(os, "cpu_count", lambda: 36)
    assert cpus.resolve_threads(0) == 1


def test_explicit_threads_are_honoured_exactly(monkeypatch: pytest.MonkeyPatch) -> None:
    """Not clamped: oversubscription is sometimes deliberate."""
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
    assert cpus.resolve_threads(16) == 16


@pytest.mark.parametrize("value", [-1, 0])
def test_non_positive_threads_mean_auto(value: int, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "7")
    assert cpus.resolve_threads(value) == 7


# --- oversubscription warning -----------------------------------------


def test_warns_when_threads_exceed_the_allocation(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    with caplog.at_level(logging.WARNING):
        cpus.warn_if_oversubscribed(16, what="annotate")
    assert "exceeds the 2 CPU(s)" in caplog.text
    assert "SLURM_CPUS_PER_TASK" in caplog.text
    assert "--threads 2" in caplog.text  # actionable suggestion


def test_no_warning_when_threads_fit(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
    with caplog.at_level(logging.WARNING):
        cpus.warn_if_oversubscribed(16, what="annotate")
    assert caplog.text == ""


def test_no_warning_for_auto(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """--threads 0 resolves to the usable count, so it can't oversubscribe."""
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
    with caplog.at_level(logging.WARNING):
        cpus.warn_if_oversubscribed(0, what="annotate")
    assert caplog.text == ""


def test_warning_does_not_raise(monkeypatch: pytest.MonkeyPatch) -> None:
    """It is advice, not a limit -- an M1 Max user may want -t 16 on 10 cores."""
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
    cpus.warn_if_oversubscribed(999, what="annotate")  # must not raise
