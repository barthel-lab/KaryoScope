"""Available-memory accounting for the steps that load a large index.

KaryoScope already refuses to start when the *disk* is too small
(:mod:`karyoscope.diskspace`). Memory had no equivalent, and it is the
easier of the two to get wrong: the HKS backend holds the index resident,
so a human-scale database needs ~10 GB whatever is being annotated, and a
run that under-requests is killed by the kernel — mid-query, with no
warning, after the index load has already been paid for.

The requirement is unusually knowable here. It is not a heuristic scaled
from the input like the output-size estimate; it is the size of the index
files themselves, readable from disk before anything runs. See
:func:`karyoscope.core.io.hks.estimate_hks_memory_bytes`.

Fail-open by design
===================

Every resolver below can legitimately return nothing — macOS has no
``/proc``, a bare machine has no cgroup, a user may run outside SLURM. When
the limit cannot be determined, :func:`require_memory` logs and returns
rather than blocking: refusing to run because we could not measure
something would break working setups, while letting it through costs at
worst the OOM message the tool would have produced anyway (which
:mod:`karyoscope.core.external` now makes legible).
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

from karyoscope.exceptions import InsufficientMemoryError

logger = logging.getLogger(__name__)

#: Headroom over the measured requirement. The index is not the whole
#: footprint — there are query buffers, the allocator's slack, and the
#: Python parent — and a run driven to exactly its limit is killed rather
#: than slowed. Measured overhead above the index size was ~0.3-1.0 GB on
#: human data, comfortably inside this.
DEFAULT_MARGIN = 0.15

#: Bytes per decimal GB, matching :mod:`karyoscope.diskspace`.
GB = 1_000_000_000

_CGROUP_V2_MAX = Path("/sys/fs/cgroup/memory.max")
_CGROUP_V1_MAX = Path("/sys/fs/cgroup/memory/memory.limit_in_bytes")
_MEMINFO = Path("/proc/meminfo")

#: A cgroup limit at or above this is the kernel's "unlimited" sentinel in
#: disguise rather than a real cap, so it is ignored. cgroup v1 reports
#: something near 2^63; v2 writes the literal string ``max``.
_EFFECTIVELY_UNLIMITED = 1 << 62


def _read_int(path: Path) -> int | None:
    try:
        text = path.read_text().strip()
    except OSError:
        return None
    if text == "max":
        return None
    try:
        value = int(text)
    except ValueError:
        return None
    return value if 0 < value < _EFFECTIVELY_UNLIMITED else None


def _slurm_allocation_bytes() -> int | None:
    """Memory SLURM granted this job, if it says so.

    Checked before the cgroup because a site can schedule by memory without
    enforcing it via cgroups, in which case the allocation is real but
    invisible to every other resolver — and exceeding it still gets the job
    killed.
    """
    per_node = os.environ.get("SLURM_MEM_PER_NODE")
    if per_node and per_node.isdigit():
        return int(per_node) * 1024 * 1024  # SLURM reports MB
    per_cpu = os.environ.get("SLURM_MEM_PER_CPU")
    cpus = os.environ.get("SLURM_CPUS_PER_TASK") or os.environ.get("SLURM_JOB_CPUS_PER_NODE")
    if per_cpu and per_cpu.isdigit() and cpus and cpus.isdigit():
        return int(per_cpu) * int(cpus) * 1024 * 1024
    return None


def _cgroup_limit_bytes() -> int | None:
    """The cgroup memory ceiling, which is what actually kills the process.

    Covers containers and SLURM sites that bind with cgroups. v2 first,
    since v1 lingers only on older kernels.
    """
    return _read_int(_CGROUP_V2_MAX) or _read_int(_CGROUP_V1_MAX)


def _meminfo_available_bytes() -> int | None:
    """``MemAvailable`` from ``/proc/meminfo``.

    Deliberately ``MemAvailable`` rather than ``MemFree``: the kernel's own
    estimate of what a new allocation could obtain, counting reclaimable
    page cache. ``MemFree`` on a busy node reads near zero because the cache
    is doing its job, and would reject every run.
    """
    try:
        for line in _MEMINFO.read_text().splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) * 1024  # kB
    except (OSError, ValueError, IndexError):
        return None
    return None


def available_bytes() -> tuple[int, str] | None:
    """Return ``(bytes, source)`` for the memory this process may use.

    Resolved in order of authority, mirroring :mod:`karyoscope.cpus`:

    1. ``$SLURM_MEM_PER_NODE`` / ``$SLURM_MEM_PER_CPU`` — what the scheduler
       granted, correct even where cgroups do not enforce it.
    2. The cgroup limit — what will actually kill the process, in a
       container or a cgroup-bound job.
    3. ``MemAvailable`` — the machine, for an unconstrained run.

    ``None`` means "could not determine", never "no memory".
    """
    slurm = _slurm_allocation_bytes()
    if slurm is not None:
        return slurm, "$SLURM_MEM_PER_NODE/$SLURM_MEM_PER_CPU"
    cgroup = _cgroup_limit_bytes()
    if cgroup is not None:
        return cgroup, "the cgroup memory limit"
    meminfo = _meminfo_available_bytes()
    if meminfo is not None:
        return meminfo, "MemAvailable in /proc/meminfo"
    return None


def require_memory(
    required_bytes: int,
    *,
    what: str,
    margin: float = DEFAULT_MARGIN,
    hint: str | None = None,
    skip: bool = False,
) -> None:
    """Raise unless this process may use ``required_bytes`` (plus margin).

    Parameters
    ----------
    required_bytes
        The measured requirement — for HKS, the size of the index files
        that must be resident.
    what
        What the memory is for, used in the message.
    margin
        Fractional headroom added on top of ``required_bytes``.
    hint
        Extra advice appended to the failure message.
    skip
        If True, log and return without raising. Wired to the same
        ``--no-space-check`` escape as the disk check.

    Raises
    ------
    InsufficientMemoryError
        Only when the available memory could actually be determined and is
        below the requirement. An undeterminable limit is never an error.
    """
    needed = int(required_bytes * (1 + margin))
    resolved = available_bytes()

    if skip:
        logger.info("memory check skipped for %s: need ~%.1f GB", what, needed / GB)
        return
    if resolved is None:
        logger.debug(
            "could not determine available memory; not enforcing the ~%.1f GB %s needs",
            needed / GB,
            what,
        )
        return

    have, source = resolved
    if have >= needed:
        logger.info(
            "memory check passed for %s: need ~%.1f GB, %.1f GB available (per %s)",
            what,
            needed / GB,
            have / GB,
            source,
        )
        return

    lines = [
        f"not enough memory for {what}: "
        f"needs ~{needed / GB:.1f} GB but only {have / GB:.1f} GB is available "
        f"(per {source}); short by {(needed - have) / GB:.1f} GB."
    ]
    if hint:
        lines.append(hint)
    lines.append(
        "This is checked up front because the alternative is the kernel "
        "killing the run mid-query with no explanation."
    )
    raise InsufficientMemoryError("\n".join(lines))
