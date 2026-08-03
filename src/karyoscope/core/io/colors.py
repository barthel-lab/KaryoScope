"""Parser for KaryoScope ``colors.txt`` files.

The file lives inside each database directory and pins one hex colour
per ``(feature_set, feature)`` pair. Format (tab-separated, header row
required)::

    feature_set\tfeature\tcolor
    chromosome\tchr1\t#1f77b4
    chromosome\tchr2\t#ff7f0e
    region\trA\t#2ca02c
    ...

The same feature name can appear in multiple feature sets with
different colours (e.g. ``"acrocentric"`` may have one colour in the
``chromosome`` set and another in a ``region`` set), so the parser
returns a nested ``{feature_set: {feature: hex}}`` mapping. The
karyotype renderer flattens this to a single-feature-set view at
render time via :func:`colors_for_set`.

Two implicit colour conventions are kept consistent with the archive:

* ``novel`` always renders as white (``#ffffff``), whether or not the
  colours file declares it explicitly. The renderer applies this
  fallback inline so the colours file doesn't have to repeat it.
* Features absent from the file render as white (``#ffffff``) with a
  one-time warning per feature. This stops a missing-colour entry
  from being a hard error -- the figure still renders, just with a
  blank patch the user can investigate.
"""

from __future__ import annotations

from pathlib import Path

from karyoscope.core.io.features import NOVEL_NAME
from karyoscope.core.io.hierarchy import Hierarchy
from karyoscope.exceptions import KaryoscopeError


class ColorsError(KaryoscopeError):
    """Problems parsing a ``colors.txt`` file."""


_REQUIRED_HEADER: tuple[str, ...] = ("feature_set", "feature", "color")

#: Optional 4th column. When present, features sharing a value are collapsed
#: into one legend row labelled with it. Optional so every existing 3-column
#: file keeps parsing unchanged.
#:
#: The point is legibility of large feature sets: the CHM13 cytoband database
#: has 1105 features in 9 colours, and a per-feature legend both dwarfs the
#: figure and gets silently truncated to whatever fits the canvas. Grouping by
#: Giemsa stain turns 833 rendered features into 8 rows.
_OPTIONAL_GROUP_COLUMN = "legend_group"
_HEADER_WITH_GROUP: tuple[str, ...] = (*_REQUIRED_HEADER, _OPTIONAL_GROUP_COLUMN)


def parse_colors_and_groups(
    path: Path,
) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    """Parse ``colors.txt`` into ``(colours, legend_groups)``.

    Both are ``{feature_set: {feature: value}}``. ``legend_groups`` is empty
    for a 3-column file, which is every database shipped before the optional
    column existed.

    Validates the header and column count. Hex colour validation is
    intentionally loose -- any non-empty value is accepted, since the
    SVG layer renders it as-is and an invalid colour is harmless
    (just doesn't render).
    """
    if not path.is_file():
        raise ColorsError(f"colors file not found: {path}")

    try:
        text = path.read_text()
    except OSError as e:
        raise ColorsError(f"could not read colors file {path}: {e}") from e

    lines = text.splitlines()
    if not lines:
        raise ColorsError(f"colors file is empty: {path}")

    header = tuple(h.strip() for h in lines[0].split("\t"))
    if header not in (_REQUIRED_HEADER, _HEADER_WITH_GROUP):
        raise ColorsError(
            f"{path}: unexpected header. Got {header!r}, expected "
            f"{_REQUIRED_HEADER!r} or {_HEADER_WITH_GROUP!r}"
        )
    n_columns = len(header)

    out: dict[str, dict[str, str]] = {}
    groups: dict[str, dict[str, str]] = {}
    for i, raw in enumerate(lines[1:], start=2):
        if not raw.strip():
            continue
        parts = raw.split("\t")
        if len(parts) != n_columns:
            raise ColorsError(
                f"{path}:{i}: expected {n_columns} tab-separated columns, got {len(parts)}: {raw!r}"
            )
        fs, feature, color = (p.strip() for p in parts[:3])
        if not fs or not feature or not color:
            raise ColorsError(
                f"{path}:{i}: feature_set, feature and color must be non-empty: {raw!r}"
            )
        out.setdefault(fs, {})[feature] = color
        if n_columns == 4:
            group = parts[3].strip()
            if group:
                groups.setdefault(fs, {})[feature] = group
    return out, groups


def parse_colors(path: Path) -> dict[str, dict[str, str]]:
    """Parse ``colors.txt`` at ``path`` into ``{feature_set: {feature: hex_color}}``.

    The colours alone, for the many callers that do not care about legend
    grouping. See :func:`parse_colors_and_groups` for both halves.
    """
    colors, _ = parse_colors_and_groups(path)
    return colors


def colors_for_set(
    colors: dict[str, dict[str, str]],
    feature_set: str,
) -> dict[str, str]:
    """Return ``{feature: hex_color}`` for one feature set, with the
    implicit ``novel -> #ffffff`` sentinel pre-populated.

    Callers that look up features by name can layer their own fallback
    on top (typically "unknown feature -> #ffffff with a warning").
    """
    out: dict[str, str] = {"novel": "#ffffff"}
    out.update(colors.get(feature_set, {}))
    return out


def validate_colors(
    hierarchy: Hierarchy,
    colors: dict[str, dict[str, str]],
) -> list[str]:
    """Cross-check colors.tsv against hierarchy.tsv.

    Returns a list of human-readable issue strings (empty list when
    the colors file is complete). Each issue is one missing
    ``(feature_set, node)`` pair.

    The rule: every node that appears in ``hierarchy.tsv`` for a
    given feature set -- whether as a child or as a parent -- must
    have an entry in ``colors.tsv`` for that feature set. Smoothing
    can promote intervals to any node in the tree (including the
    ``categorized`` root, which only ever appears in the parent
    column), so all of them can show up in rendered output.

    The implicit ``novel`` sentinel is **not** required in
    ``colors.tsv`` -- it's always rendered white by the renderer.

    Used by:

    * :func:`karyoscope.commands.info.cmd` as a soft warning.
    * :func:`karyoscope.commands.download.cmd` as a hard error so
      community-built databases that are missing colors get rejected
      at install time, not karyotype time.
    * :func:`karyoscope.core.karyotype_run.karyotype_run` as a hard
      error so locally-installed databases that bypassed download
      validation still get caught before any SVG is rendered.
    """
    issues: list[str] = []
    for feature_set in hierarchy.feature_sets():
        nodes = hierarchy.nodes(feature_set)
        fs_colors = colors.get(feature_set, {})
        for node in sorted(nodes):
            if node == NOVEL_NAME:
                # Defensive: ``novel`` is the white sentinel and
                # should never appear in hierarchy.tsv. If it does
                # we let it through here (hierarchy validation is
                # the right place to flag that).
                continue
            if node not in fs_colors:
                issues.append(
                    f"feature set {feature_set!r}: hierarchy node {node!r} "
                    "has no entry in colors.tsv"
                )
    return issues


def validate_legend_groups(
    colors: dict[str, dict[str, str]],
    groups: dict[str, dict[str, str]],
) -> list[str]:
    """Check the invariant the legend grouping rests on: one colour per group.

    Returns human-readable issue strings, empty when the file is consistent
    (which includes every 3-column file, since those declare no groups).

    A legend row is one swatch and one label, so a group spanning two colours
    has no well-defined swatch -- the renderer would have to pick one and
    silently misrepresent the rest. That is exactly the failure the grouping
    exists to prevent, so it is caught when the database is validated rather
    than discovered in a figure.

    The reverse is deliberately allowed: two groups may share a colour. That
    is a legible figure (two labels, same swatch), just a redundant one.
    """
    issues: list[str] = []
    for feature_set, per_feature in sorted(groups.items()):
        colors_by_group: dict[str, dict[str, str]] = {}
        for feature, group in per_feature.items():
            color = colors.get(feature_set, {}).get(feature)
            if color is None:
                continue  # absent colour is already validate_colors' business
            colors_by_group.setdefault(group, {})[feature] = color
        for group, members in sorted(colors_by_group.items()):
            distinct = sorted(set(members.values()))
            if len(distinct) > 1:
                examples = ", ".join(f"{f}={c}" for f, c in sorted(members.items())[:4])
                issues.append(
                    f"feature set {feature_set!r}: legend group {group!r} spans "
                    f"{len(distinct)} colours ({', '.join(distinct)}) -- a legend "
                    f"row has one swatch, so every feature in a group must share "
                    f"a colour. Examples: {examples}"
                )
    return issues
