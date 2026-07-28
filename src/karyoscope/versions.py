"""Comparing dotted version strings.

Two unrelated checks need the same ordering: whether this KaryoScope is new
enough for a database (``karyoscope_min_version``), and whether an external
binary is new enough for the flags KaryoScope passes it. Both want "is A at
least B" over strings like ``2.1.0``, and neither wants a dependency on
``packaging`` for it.
"""

from __future__ import annotations


def parse_version(s: str) -> tuple[int, ...]:
    """Parse a dotted version string into a tuple of ints, for comparison.

    Non-numeric components are treated as 0 (so ``1.0.0.dev0`` < ``1.0.0``).
    This is a small intentional subset of PEP 440 — enough to order the
    release versions we actually compare, not a general-purpose comparator.
    """
    parts: list[int] = []
    for part in s.split("."):
        # Strip any non-digit suffix (e.g., "0dev0" -> 0).
        n = 0
        for ch in part:
            if not ch.isdigit():
                break
            n = n * 10 + int(ch)
        parts.append(n)
    return tuple(parts)


def at_least(have: str, required: str) -> bool:
    """Is version ``have`` greater than or equal to ``required``?

    Shorter versions compare as if zero-padded, because tuple comparison
    already does that: ``(0, 3) < (0, 3, 1)``.
    """
    return parse_version(have) >= parse_version(required)
