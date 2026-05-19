"""Parser for KaryoScope ``hierarchy.tsv`` files.

A hierarchy file is a 3-column tab-separated file with a header row:

    feature_set	feature	parent

* ``feature_set`` is the named partition the feature belongs to
  (e.g., ``chromosome``, ``region``, ``centromere``).
* ``feature`` is a unique name within its feature set.
* ``parent`` is either another feature (the parent of this one for
  smoothing purposes) or ``.`` to indicate "no parent" — i.e., this
  feature is a root of the hierarchy.

Example::

    feature_set	feature	parent
    chromosome	chr1	.
    chromosome	chr2	.
    region	1p	chr1
    region	1q	chr1
    region	2p	chr2
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from karyoscope.exceptions import KaryoscopeError


class HierarchyError(KaryoscopeError):
    """Problems parsing a ``hierarchy.tsv`` file."""


#: Sentinel used in the parent column to mean "no parent" (this row is a root).
NO_PARENT = "."


@dataclass(frozen=True)
class HierarchyRow:
    """One row in a hierarchy file."""

    feature_set: str
    feature: str
    parent: str | None  # ``None`` when the parent column was ``.``


@dataclass
class Hierarchy:
    """The parsed contents of a ``hierarchy.tsv`` file."""

    rows: list[HierarchyRow]

    def feature_sets(self) -> list[str]:
        """Return all feature-set names, in the order they first appear."""
        seen: list[str] = []
        for row in self.rows:
            if row.feature_set not in seen:
                seen.append(row.feature_set)
        return seen

    def features_in(self, feature_set: str) -> list[str]:
        """Return all features that belong to ``feature_set``, in input order."""
        return [r.feature for r in self.rows if r.feature_set == feature_set]

    def count_by_feature_set(self) -> dict[str, int]:
        """Return ``{feature_set: count}`` for every feature set."""
        counts: dict[str, int] = {}
        for row in self.rows:
            counts[row.feature_set] = counts.get(row.feature_set, 0) + 1
        return counts

    def roots(self, feature_set: str | None = None) -> list[HierarchyRow]:
        """Return all rows with no parent.

        If ``feature_set`` is given, restrict to that feature set.
        """
        return [
            r
            for r in self.rows
            if r.parent is None and (feature_set is None or r.feature_set == feature_set)
        ]


_REQUIRED_HEADER = ("feature_set", "feature", "parent")


def parse_hierarchy(path: Path) -> Hierarchy:
    """Parse ``hierarchy.tsv`` at ``path`` and return a :class:`Hierarchy`.

    Raises :class:`HierarchyError` on:

    * missing file
    * empty file or missing header
    * wrong column count on any line
    * unexpected header columns
    """
    if not path.is_file():
        raise HierarchyError(f"hierarchy file not found: {path}")

    try:
        text = path.read_text()
    except OSError as e:
        raise HierarchyError(f"could not read hierarchy file {path}: {e}") from e

    lines = text.splitlines()
    if not lines:
        raise HierarchyError(f"hierarchy file is empty: {path}")

    header = lines[0].split("\t")
    if tuple(h.strip() for h in header) != _REQUIRED_HEADER:
        raise HierarchyError(
            f"hierarchy file {path} has unexpected header: "
            f"got {header!r}, expected {list(_REQUIRED_HEADER)!r}"
        )

    rows: list[HierarchyRow] = []
    for i, raw in enumerate(lines[1:], start=2):
        # Skip blank lines (tolerate trailing newlines and stray empties).
        if not raw.strip():
            continue
        parts = raw.split("\t")
        if len(parts) != 3:
            raise HierarchyError(
                f"{path}:{i}: expected 3 tab-separated columns, got {len(parts)}: {raw!r}"
            )
        feature_set, feature, parent_raw = parts
        parent: str | None = None if parent_raw == NO_PARENT else parent_raw
        rows.append(
            HierarchyRow(
                feature_set=feature_set,
                feature=feature,
                parent=parent,
            )
        )
    return Hierarchy(rows=rows)
