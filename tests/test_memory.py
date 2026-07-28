"""Tests for :mod:`karyoscope.memory`.

The check exists because an under-allocated HKS run is killed by the kernel
mid-query with no explanation. Its whole value is being right about when to
block — so most of what is pinned here is when it must *not*.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope import memory
from karyoscope.exceptions import InsufficientMemoryError

GB = memory.GB


def _no_resolvers(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)
    monkeypatch.setattr(memory, "_cgroup_limit_bytes", lambda: None)
    monkeypatch.setattr(memory, "_meminfo_available_bytes", lambda: None)


# --- resolution order -------------------------------------------------


def test_slurm_allocation_wins(monkeypatch: pytest.MonkeyPatch) -> None:
    """The scheduler's grant outranks the cgroup and the machine.

    A site can schedule by memory without cgroup enforcement, in which case
    the allocation is real but invisible to everything else — and exceeding
    it still gets the job killed.
    """
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "16384")  # MB
    monkeypatch.setattr(memory, "_cgroup_limit_bytes", lambda: 999 * GB)
    monkeypatch.setattr(memory, "_meminfo_available_bytes", lambda: 999 * GB)
    have, source = memory.available_bytes()
    assert have == 16384 * 1024 * 1024
    assert "SLURM" in source


def test_slurm_per_cpu_is_multiplied_by_the_cpu_count(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
    monkeypatch.setenv("SLURM_MEM_PER_CPU", "2048")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setattr(memory, "_cgroup_limit_bytes", lambda: None)
    have, _ = memory.available_bytes()
    assert have == 2048 * 8 * 1024 * 1024


def test_cgroup_is_used_when_slurm_is_silent(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)
    monkeypatch.setattr(memory, "_cgroup_limit_bytes", lambda: 12 * GB)
    monkeypatch.setattr(memory, "_meminfo_available_bytes", lambda: 999 * GB)
    have, source = memory.available_bytes()
    assert have == 12 * GB
    assert "cgroup" in source


def test_unlimited_cgroup_values_are_ignored(tmp_path: Path) -> None:
    """cgroup v2 writes 'max'; v1 writes a number near 2**63.

    Treating either as a real limit would report absurd availability, and
    (worse) would stop the resolver falling through to MemAvailable.
    """
    v2 = tmp_path / "memory.max"
    v2.write_text("max\n")
    assert memory._read_int(v2) is None
    v1 = tmp_path / "limit_in_bytes"
    v1.write_text(str(2**63 - 1))
    assert memory._read_int(v1) is None
    real = tmp_path / "real"
    real.write_text(str(8 * GB))
    assert memory._read_int(real) == 8 * GB


# --- fail-open --------------------------------------------------------


def test_undeterminable_memory_never_blocks(monkeypatch: pytest.MonkeyPatch) -> None:
    """macOS has no /proc, a bare machine has no cgroup, users run outside SLURM.

    Refusing to run because we could not measure something would break
    working setups; letting it through costs at worst the tool's own OOM
    message, which is now legible.
    """
    _no_resolvers(monkeypatch)
    assert memory.available_bytes() is None
    memory.require_memory(500 * GB, what="testing")  # must not raise


def test_skip_bypasses_the_check(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(memory, "available_bytes", lambda: (1 * GB, "test"))
    memory.require_memory(500 * GB, what="testing", skip=True)


def test_enough_memory_passes(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(memory, "available_bytes", lambda: (64 * GB, "test"))
    memory.require_memory(10 * GB, what="testing")


# --- the failure it exists to produce ---------------------------------


def test_too_little_memory_raises_with_the_numbers(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(memory, "available_bytes", lambda: (4 * GB, "$SLURM_MEM_PER_NODE"))
    with pytest.raises(InsufficientMemoryError) as excinfo:
        memory.require_memory(
            10 * GB, what="annotating HG002 against HKS_human_CHM13_v2", hint="  - use --mem=16G"
        )
    msg = str(excinfo.value)
    assert "annotating HG002 against HKS_human_CHM13_v2" in msg
    assert "4.0 GB" in msg  # what we have
    assert "$SLURM_MEM_PER_NODE" in msg  # and where that came from
    assert "short by" in msg
    assert "--mem=16G" in msg  # actionable


def test_the_margin_is_applied(monkeypatch: pytest.MonkeyPatch) -> None:
    """A run driven to exactly its limit is killed, not slowed.

    The index is not the whole footprint -- query buffers, allocator slack
    and the Python parent sit on top.
    """
    monkeypatch.setattr(memory, "available_bytes", lambda: (10 * GB, "test"))
    memory.require_memory(10 * GB, what="testing", margin=0.0)
    with pytest.raises(InsufficientMemoryError):
        memory.require_memory(10 * GB, what="testing", margin=0.15)


# --- the HKS requirement is measured, not extrapolated ----------------


def test_hks_requirement_is_base_plus_largest_labeling(tmp_path: Path) -> None:
    """hks loads the base plus ONE labeling at a time, so the peak is that sum.

    This is what makes the figure independent of input size and thread count.
    """
    from karyoscope.core.io.hks import estimate_hks_memory_bytes

    (tmp_path / "features.hksb").write_bytes(b"x" * 1000)
    (tmp_path / "features.chromosome.hksf").write_bytes(b"x" * 500)
    (tmp_path / "features.repeat.hksf").write_bytes(b"x" * 300)
    got = estimate_hks_memory_bytes(
        db_dir=tmp_path, basename="features", feature_sets=["chromosome", "repeat"]
    )
    assert got == 1000 + 500  # base + LARGEST labeling, not the sum of both

    # Requesting only the smaller set lowers it.
    assert (
        estimate_hks_memory_bytes(
            db_dir=tmp_path, basename="features", feature_sets=["repeat"]
        )
        == 1000 + 300
    )


def test_hks_requirement_is_none_when_the_index_cannot_be_sized(tmp_path: Path) -> None:
    """Fail open: a database we cannot inspect must not block the run."""
    from karyoscope.core.io.hks import estimate_hks_memory_bytes

    assert (
        estimate_hks_memory_bytes(
            db_dir=tmp_path, basename="features", feature_sets=["chromosome"]
        )
        is None
    )
    (tmp_path / "features.hksb").write_bytes(b"x" * 10)
    assert (
        estimate_hks_memory_bytes(
            db_dir=tmp_path, basename="features", feature_sets=["missing"]
        )
        is None
    )
