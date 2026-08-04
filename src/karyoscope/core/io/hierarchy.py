"""Parser for KaryoScope ``hierarchy.tsv`` files.

A hierarchy file declares the directed tree of feature relationships
within each feature set. Each feature set has its own independent
hierarchy; the file format is::

    feature_set	child	parent
    chromosome	autosome	categorized
    chromosome	chr1	autosome
    chromosome	chr2	autosome
    region	centromeric	categorized
    region	aSat	centromeric
    region	HSat	centromeric

Rows declare ``child → parent`` edges within the ``feature_set`` named in
the first column. Each feature set must form a single connected tree
whose root is named ``"categorized"`` (the v0.1 convention; we may
relax this in a later version). Node names can repeat across feature
sets — those are different nodes living in different trees — but
within a single feature set every node has at most one parent.

The same name can appear as both a child and a parent within one
feature set, on different rows (e.g., the acrocentric set has rows
declaring ``acrocentric → categorized`` and ``array → acrocentric``).
That's fine; the parser doesn't treat the visual repetition specially.

Two layers:

* :func:`parse_hierarchy` is a pure structural parser. It validates the
  file format (header, column counts) but does *not* validate tree
  shape — connectedness, cycles, root naming, etc.
* :func:`validate_hierarchy` is a separate validation pass. It takes a
  parsed :class:`Hierarchy` and (optionally) a feature-name index, and
  reports a list of issues. Callers decide what to do with the
  issues: ``info`` prints them as warnings, ``annotate`` raises.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from karyoscope.exceptions import KaryoscopeError


class HierarchyError(KaryoscopeError):
    """Problems parsing a ``hierarchy.tsv`` file."""


#: The required root node name (v0.1 convention).
REQUIRED_ROOT = "categorized"

_REQUIRED_HEADER = ("feature_set", "child", "parent")


@dataclass(frozen=True)
class HierarchyRow:
    """One edge in a hierarchy file: ``child → parent`` within ``feature_set``."""

    feature_set: str
    child: str
    parent: str


@dataclass
class Hierarchy:
    """The parsed contents of a ``hierarchy.tsv`` file.

    Provides per-feature-set lookups. The heavier algorithmic queries
    (ancestors, LCA, depth) live on
    :class:`karyoscope.core.smooth.HierarchyIndex`, which is built from
    a :class:`Hierarchy` on demand.
    """

    rows: list[HierarchyRow] = field(default_factory=list)

    def feature_sets(self) -> list[str]:
        """Return all feature-set names, in the order they first appear."""
        seen: list[str] = []
        for row in self.rows:
            if row.feature_set not in seen:
                seen.append(row.feature_set)
        return seen

    def rows_in(self, feature_set: str) -> list[HierarchyRow]:
        """Return all rows whose first column is ``feature_set``."""
        return [r for r in self.rows if r.feature_set == feature_set]

    def parent_map(self, feature_set: str) -> dict[str, str]:
        """Return ``{child: parent}`` for one feature set.

        The root node (which has no parent edge) does not appear as a
        key. Raises :class:`HierarchyError` if a child appears twice
        within this feature set (would silently lose data otherwise).
        """
        out: dict[str, str] = {}
        for row in self.rows_in(feature_set):
            if row.child in out:
                raise HierarchyError(
                    f"feature set {feature_set!r} has duplicate child {row.child!r}"
                )
            out[row.child] = row.parent
        return out

    def nodes(self, feature_set: str) -> set[str]:
        """Return the set of all node names referenced in ``feature_set``.

        Includes both child and parent names (so the root, which appears
        only as a parent, is included).
        """
        out: set[str] = set()
        for row in self.rows_in(feature_set):
            out.add(row.child)
            out.add(row.parent)
        return out

    def nodes_in_order(self, feature_set: str) -> list[str]:
        """Return every node in ``feature_set`` in hierarchy order: roots, then
        each child in edge order.

        The ordered counterpart to :meth:`nodes`. Anything written out per node
        must use this rather than iterating the set, or the row order varies
        between runs — Python randomises string hashing, so a rebuild from
        identical inputs produced a differently-ordered ``colors.tsv``. Row
        order is also load-bearing: legend groups are ordered by first
        appearance in ``colors.tsv``, so hierarchy order is what makes that
        ordering meaningful (for cytoband, the Giemsa intensity progression)
        instead of arbitrary.
        """
        rows = self.rows_in(feature_set)
        children = {row.child for row in rows}
        ordered: list[str] = []
        seen: set[str] = set()
        # Roots first: a parent that is never a child. Well-formed sets have
        # exactly one, but this does not depend on that.
        for row in rows:
            if row.parent not in children and row.parent not in seen:
                ordered.append(row.parent)
                seen.add(row.parent)
        for row in rows:
            if row.child not in seen:
                ordered.append(row.child)
                seen.add(row.child)
        return ordered

    def count_by_feature_set(self) -> dict[str, int]:
        """Return ``{feature_set: row_count}`` for every feature set."""
        counts: dict[str, int] = {}
        for row in self.rows:
            counts[row.feature_set] = counts.get(row.feature_set, 0) + 1
        return counts


def parse_hierarchy(path: Path) -> Hierarchy:
    """Parse ``hierarchy.tsv`` at ``path`` and return a :class:`Hierarchy`.

    Structural validation only — column count, header. Tree-shape
    validation (single root, connectedness, no cycles) is performed by
    :func:`validate_hierarchy`.

    Raises
    ------
    HierarchyError
        On missing file, empty file, malformed header, wrong column
        count on any row, or blank values in any column.
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
        if not raw.strip():
            continue
        parts = raw.split("\t")
        if len(parts) != 3:
            raise HierarchyError(
                f"{path}:{i}: expected 3 tab-separated columns, got {len(parts)}: {raw!r}"
            )
        feature_set, child, parent = (p.strip() for p in parts)
        if not feature_set or not child or not parent:
            raise HierarchyError(f"{path}:{i}: all three columns must be non-empty: {raw!r}")
        rows.append(HierarchyRow(feature_set=feature_set, child=child, parent=parent))
    return Hierarchy(rows=rows)


# --- validation ------------------------------------------------------


def validate_hierarchy(
    hierarchy: Hierarchy,
    *,
    feature_columns: dict[str, set[str]] | None = None,
) -> list[str]:
    """Validate a parsed hierarchy and return a list of issues.

    Checks performed (per feature set, independently):

    1. **Single parent**: no child appears as the left-hand side of two
       rows within one feature set.
    2. **Single root**: exactly one node has no parent within the set,
       and that node is named :data:`REQUIRED_ROOT`.
    3. **Connectedness / no cycles**: walking parent pointers from any
       node reaches the root in a bounded number of steps.
    4. **features.tsv consistency (optional)**: if ``feature_columns``
       is given, every feature name appearing in features.tsv's column
       for this set must be a node in the hierarchy. Pass
       ``{feature_set: {name, name, ...}}`` for the sets you want
       checked.

    Returns
    -------
    list[str]
        Human-readable issue messages. An empty list means the
        hierarchy is well-formed. Each issue is self-contained and
        suitable for printing as a warning or assembling into an error
        message.
    """
    issues: list[str] = []
    for feature_set in hierarchy.feature_sets():
        issues.extend(
            _validate_one_set(
                hierarchy,
                feature_set,
                feature_columns=feature_columns,
            )
        )
    return issues


def _validate_one_set(
    hierarchy: Hierarchy,
    feature_set: str,
    *,
    feature_columns: dict[str, set[str]] | None,
) -> list[str]:
    """Run validation checks for a single feature set's subtree."""
    issues: list[str] = []
    rows = hierarchy.rows_in(feature_set)
    if not rows:
        # Vacuously valid — caller probably shouldn't have asked.
        return issues

    # Check 1: single parent per child
    parent_of: dict[str, str] = {}
    seen_children: set[str] = set()
    for row in rows:
        if row.child in seen_children:
            issues.append(
                f"feature set {feature_set!r}: child {row.child!r} has multiple "
                "parent rows (each node must have at most one parent)"
            )
            continue
        seen_children.add(row.child)
        parent_of[row.child] = row.parent

    # All node names referenced (children + parents)
    all_nodes: set[str] = set(parent_of.keys()) | set(parent_of.values())

    # Check 2: single root, named REQUIRED_ROOT
    children = set(parent_of.keys())
    roots = all_nodes - children
    if len(roots) == 0:
        issues.append(
            f"feature set {feature_set!r}: no root found (every node has a "
            "parent — this indicates a cycle)"
        )
    elif len(roots) > 1:
        issues.append(
            f"feature set {feature_set!r}: multiple roots found "
            f"{sorted(roots)!r}; exactly one is required"
        )
    else:
        only_root = next(iter(roots))
        if only_root != REQUIRED_ROOT:
            issues.append(
                f"feature set {feature_set!r}: root is {only_root!r}, must be {REQUIRED_ROOT!r}"
            )

    # Check 3: no cycles / all nodes reachable from root.
    # Walk parent pointers from every child, with a bound to detect
    # cycles. If we walk more than (n_nodes) steps without reaching a
    # node with no parent, there's a cycle.
    bound = len(all_nodes) + 1
    for child in parent_of:
        node = child
        seen_on_walk: set[str] = set()
        cycle_found = False
        for _ in range(bound):
            if node not in parent_of:
                break
            if node in seen_on_walk:
                issues.append(
                    f"feature set {feature_set!r}: cycle detected involving node {node!r}"
                )
                cycle_found = True
                break
            seen_on_walk.add(node)
            node = parent_of[node]
        else:
            if not cycle_found:
                issues.append(
                    f"feature set {feature_set!r}: walk from {child!r} did not "
                    "terminate within bounded steps (likely cycle)"
                )

    # Check 4: features.tsv consistency
    if feature_columns is not None and feature_set in feature_columns:
        missing = sorted(feature_columns[feature_set] - all_nodes)
        if missing:
            shown = missing if len(missing) <= 5 else [*missing[:5], "..."]
            issues.append(
                f"feature set {feature_set!r}: {len(missing)} feature name(s) "
                f"in features.tsv have no row in hierarchy.tsv: {shown}"
            )

    return issues
