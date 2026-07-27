"""How many CPUs this process may actually use.

``os.cpu_count()`` answers a different question than the one KaryoScope
needs. It reports the machine's logical CPUs and ignores every mechanism
that restricts a process to a subset of them — cgroups, CPU affinity, and
in particular a SLURM allocation. On a shared 36-core node with
``--cpus-per-task=1``, ``os.cpu_count()`` returns 36 while the job may use
exactly one:

    os.cpu_count()            : 36
    len(os.sched_getaffinity(0)) : 1
    $SLURM_JOB_CPUS_PER_NODE  : 1

Since ``--threads 0`` (the default) sized its worker pool from
``os.cpu_count()``, that default spawned 36 workers to contend for one
CPU. This module answers the question that was actually meant.

Ordering rationale
==================

Each source is consulted only when the more authoritative ones are
unavailable:

1. ``$SLURM_CPUS_PER_TASK`` — what the scheduler granted this task.
   Authoritative on our cluster, and correct even when the site runs
   without cgroup binding (in which case affinity would over-report).
2. ``os.sched_getaffinity`` — Linux only; reflects ``taskset``, cgroup
   CPU sets, and SLURM's affinity binding.
3. ``os.process_cpu_count`` — does the right thing, but is Python 3.13+
   and we support 3.10.
4. ``os.cpu_count`` — the machine. The last resort, and what macOS falls
   back to since it exposes no affinity API.

What this deliberately does not do
==================================

It does not cap ``--threads``. Oversubscription is sometimes the right
call — I/O-bound stages, and heterogeneous CPUs (an Apple M1 Max is 8
performance + 2 efficiency cores, where extra threads give the scheduler
slack to keep the fast cores fed). A hard cap would remove a legitimate
tuning knob. :func:`warn_if_oversubscribed` says something and gets out of
the way.

It also cannot see core *asymmetry*: macOS reports 10 CPUs on an M1 Max
with no way to learn that only 8 are performance cores. That is another
reason "number of CPUs" is advice rather than a limit.
"""

from __future__ import annotations

import logging
import os

logger = logging.getLogger(__name__)

#: Environment variables that describe a scheduler allocation, in order of
#: preference. ``SLURM_CPUS_PER_TASK`` is per-task and the one that matches
#: how KaryoScope jobs are submitted; ``SLURM_JOB_CPUS_PER_NODE`` is a
#: fallback that can carry a range syntax we don't try to parse beyond its
#: first integer.
_ALLOCATION_VARS: tuple[str, ...] = (
    "SLURM_CPUS_PER_TASK",
    "SLURM_JOB_CPUS_PER_NODE",
)


def _from_allocation() -> int | None:
    """CPU count from a scheduler's environment, or None if not under one."""
    for var in _ALLOCATION_VARS:
        raw = os.environ.get(var)
        if not raw:
            continue
        # SLURM_JOB_CPUS_PER_NODE can look like "4(x2)" or "2,4". Take the
        # leading integer: for a single-node task that is this node's count,
        # and a wrong-but-plausible number here is still far closer than the
        # whole machine's core count.
        digits = ""
        for ch in raw.strip():
            if ch.isdigit():
                digits += ch
            else:
                break
        if digits:
            n = int(digits)
            if n > 0:
                logger.debug("usable CPUs from $%s: %d", var, n)
                return n
    return None


def usable_cpus() -> int:
    """Return the number of CPUs this process may actually use (>= 1).

    Prefers a scheduler allocation, then OS affinity, then the machine's
    core count. Never returns 0, so callers can size a pool with it
    directly.
    """
    allocated = _from_allocation()
    if allocated is not None:
        return allocated

    getaffinity = getattr(os, "sched_getaffinity", None)
    if getaffinity is not None:
        try:
            n = len(getaffinity(0))
            if n > 0:
                logger.debug("usable CPUs from sched_getaffinity: %d", n)
                return n
        except OSError:  # pragma: no cover — platform-specific
            pass

    process_cpu_count = getattr(os, "process_cpu_count", None)
    if process_cpu_count is not None:  # Python 3.13+
        n = process_cpu_count()
        if n:
            logger.debug("usable CPUs from os.process_cpu_count: %d", n)
            return n

    return os.cpu_count() or 1


def resolve_threads(threads: int) -> int:
    """Turn a ``--threads`` value into a concrete worker count.

    ``threads <= 0`` means "auto" and resolves to :func:`usable_cpus`.
    An explicit value is honoured as given — see the module docstring for
    why this does not clamp.
    """
    if threads > 0:
        return threads
    n = usable_cpus()
    logger.debug("--threads 0 resolved to %d usable CPU(s)", n)
    return n


def warn_if_oversubscribed(threads: int, *, what: str) -> None:
    """Log a warning when an explicit ``--threads`` exceeds what's usable.

    A warning rather than an error: oversubscription is occasionally
    faster, and the figure we compare against can itself be wrong on an
    unusual setup. But the common case — a SLURM job asking for more
    threads than it was allocated — is a real performance bug the user
    almost certainly didn't intend, and it is invisible without this.
    """
    if threads <= 0:
        return
    available = usable_cpus()
    if threads <= available:
        return
    logger.warning(
        "--threads %d exceeds the %d CPU(s) this process can use%s; %s will "
        "oversubscribe. This is allowed (it can help on heterogeneous CPUs, "
        "e.g. Apple silicon's performance/efficiency mix) but is usually "
        "slower. Pass --threads %d to match, or request more CPUs.",
        threads,
        available,
        _allocation_note(),
        what,
        available,
    )


def _allocation_note() -> str:
    """Name the source of the limit, so the warning is actionable."""
    for var in _ALLOCATION_VARS:
        if os.environ.get(var):
            return f" (limited by ${var})"
    if getattr(os, "sched_getaffinity", None) is not None:
        try:
            if len(os.sched_getaffinity(0)) < (os.cpu_count() or 1):
                return " (limited by CPU affinity, not the machine's core count)"
        except OSError:  # pragma: no cover — platform-specific
            pass
    return ""
