"""Parser for KaryoScope ``features.tsv`` files.

A features file describes the meaning of each integer feature id stored
in the KMC index. The schema is one row per feature id, with one column
per feature set::

    featureID	chromosome	region	gene
    1	chr1	rA	intergenic
    2	chr1	rB	gene_X
    3	chr2	rC	intergenic

Reading the second row: a k-mer with feature id ``2`` is simultaneously
``chr1`` in the ``chromosome`` set, ``rB`` in the ``region`` set, and
``gene_X`` in the ``gene`` set. The set of feature-set names is given by
all columns after ``featureID``.

Feature id ``0`` is implicit and never appears as a row. The C++
``get_featureIDs`` helper emits 0 for k-mers not present in the KMC
index (KMC's "miss" signal); :func:`render_feature` translates this
sentinel to ``"novel"``.

Any *other* feature id that appears in a BED record but has no row in
``features.tsv`` is treated as a hard error rather than papered over.
The previous archive's ``smooth_features.py`` used a fall-back of
``"Unknown"``, but real KaryoScope feature sets (e.g., the repeats set)
can include ``"Unknown"`` as a legitimate feature name, so silently
producing it for missing ids would be ambiguous — and a missing id
genuinely indicates a database / index mismatch that the user should
hear about.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from karyoscope.exceptions import KaryoscopeError


class FeaturesError(KaryoscopeError):
    """Problems parsing or looking up entries in a ``features.tsv`` file."""


#: The mandatory first column name.
_ID_COL = "featureID"

#: Sentinel name returned for ``feature_id == 0`` (k-mer not in index).
NOVEL_NAME = "novel"


@dataclass
class Features:
    """Parsed contents of a ``features.tsv`` file.

    Attributes
    ----------
    feature_sets
        The feature-set names, in the order they appear in the header.
    table
        ``{feature_id: {feature_set: feature_name}}`` lookup table. The
        outer keys exclude ``0`` (which is implicitly "novel").
    """

    feature_sets: list[str]
    table: dict[int, dict[str, str]] = field(default_factory=dict)

    def feature_in_set(self, feature_id: int, feature_set: str) -> str | None:
        """Return the feature name for ``feature_id`` in ``feature_set``.

        Returns ``None`` for ``feature_id == 0`` (the "novel" sentinel)
        or for any id not present in the table. Use
        :func:`render_feature` if you want a user-facing string —
        ``render_feature`` raises rather than returning ``None`` for the
        missing-id case, which is the right behaviour at the BED-writing
        stage but not always what callers of this method want.

        Raises
        ------
        FeaturesError
            If ``feature_set`` is not one of the file's columns.
        """
        if feature_set not in self.feature_sets:
            raise FeaturesError(
                f"unknown feature set {feature_set!r}; valid sets are {self.feature_sets!r}"
            )
        row = self.table.get(feature_id)
        if row is None:
            return None
        return row.get(feature_set)

    def names_for_set(self, feature_set: str) -> dict[int, str]:
        """Project the table onto one feature set as a flat ``{id: name}`` dict.

        For per-record hot loops: after this one-time projection, a
        record costs a single dict lookup instead of
        :meth:`feature_in_set`'s membership validation plus two nested
        lookups. ``0`` is not in the result (it is implicitly
        ``"novel"``, as everywhere).

        Raises
        ------
        FeaturesError
            If ``feature_set`` is not one of the file's columns.
        """
        if feature_set not in self.feature_sets:
            raise FeaturesError(
                f"unknown feature set {feature_set!r}; valid sets are {self.feature_sets!r}"
            )
        return {fid: row[feature_set] for fid, row in self.table.items()}


def render_feature(
    feature_id: int,
    feature_set: str,
    features: Features,
) -> str:
    """Convert a feature id into the user-facing name in ``feature_set``.

    Translation rules:

    * ``feature_id == 0`` → ``"novel"`` (k-mer not in KMC index)
    * Otherwise, if the id is in ``features.tsv`` → that row's value in
      ``feature_set``.
    * Otherwise → :class:`FeaturesError`. A non-zero id absent from
      ``features.tsv`` indicates the BED was produced by a different
      index than the features file describes, or that ``features.tsv``
      is incomplete — both situations the user should learn about
      rather than have silently mapped to a string that could collide
      with a legitimate feature name.

    Raises
    ------
    FeaturesError
        If ``feature_id`` is positive but absent from the table, or if
        ``feature_set`` is not one of the file's columns.
    """
    if feature_id == 0:
        return NOVEL_NAME
    name = features.feature_in_set(feature_id, feature_set)
    if name is None:
        raise FeaturesError(
            f"feature id {feature_id} is not in features.tsv. "
            "This usually means the BED was produced against a different "
            "KMC index than the one features.tsv was written for, or that "
            "features.tsv is missing rows. Check that the database's "
            "index and features.tsv are from the same build."
        )
    return name


def parse_features(path: Path) -> Features:
    """Parse ``features.tsv`` at ``path`` and return a :class:`Features`.

    Raises
    ------
    FeaturesError
        On any of: missing file, empty file, missing/malformed header,
        non-integer or duplicate feature id, wrong column count on any
        row.
    """
    if not path.is_file():
        raise FeaturesError(f"features file not found: {path}")

    try:
        text = path.read_text()
    except OSError as e:
        raise FeaturesError(f"could not read features file {path}: {e}") from e

    lines = text.splitlines()
    if not lines:
        raise FeaturesError(f"features file is empty: {path}")

    header = lines[0].split("\t")
    if not header or header[0] != _ID_COL:
        raise FeaturesError(
            f"features file {path}: first column must be {_ID_COL!r}, got {header[0]!r}"
        )
    feature_sets = [h.strip() for h in header[1:]]
    if not feature_sets:
        raise FeaturesError(f"features file {path}: no feature-set columns after {_ID_COL!r}")
    if len(set(feature_sets)) != len(feature_sets):
        raise FeaturesError(
            f"features file {path}: duplicate feature-set names in header: {feature_sets}"
        )

    table: dict[int, dict[str, str]] = {}
    expected_cols = 1 + len(feature_sets)
    for i, raw in enumerate(lines[1:], start=2):
        if not raw.strip():
            continue
        parts = raw.split("\t")
        if len(parts) != expected_cols:
            raise FeaturesError(
                f"{path}:{i}: expected {expected_cols} tab-separated columns, "
                f"got {len(parts)}: {raw!r}"
            )
        try:
            feature_id = int(parts[0])
        except ValueError as e:
            raise FeaturesError(
                f"{path}:{i}: featureID column must be an integer, got {parts[0]!r}"
            ) from e
        if feature_id == 0:
            raise FeaturesError(
                f"{path}:{i}: featureID 0 is reserved for the novel sentinel "
                "and cannot appear in features.tsv"
            )
        if feature_id in table:
            raise FeaturesError(f"{path}:{i}: duplicate featureID {feature_id}")
        table[feature_id] = {fs: parts[1 + idx].strip() for idx, fs in enumerate(feature_sets)}

    return Features(feature_sets=feature_sets, table=table)
