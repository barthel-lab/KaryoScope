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

from karyoscope.exceptions import KaryoscopeError


class ColorsError(KaryoscopeError):
    """Problems parsing a ``colors.txt`` file."""


_REQUIRED_HEADER: tuple[str, ...] = ("feature_set", "feature", "color")


def parse_colors(path: Path) -> dict[str, dict[str, str]]:
    """Parse ``colors.txt`` at ``path`` into ``{feature_set: {feature: hex_color}}``.

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
    if header != _REQUIRED_HEADER:
        raise ColorsError(
            f"{path}: unexpected header. Got {header!r}, expected {_REQUIRED_HEADER!r}"
        )

    out: dict[str, dict[str, str]] = {}
    for i, raw in enumerate(lines[1:], start=2):
        if not raw.strip():
            continue
        parts = raw.split("\t")
        if len(parts) != 3:
            raise ColorsError(
                f"{path}:{i}: expected 3 tab-separated columns, got {len(parts)}: {raw!r}"
            )
        fs, feature, color = (p.strip() for p in parts)
        if not fs or not feature or not color:
            raise ColorsError(f"{path}:{i}: all three columns must be non-empty: {raw!r}")
        out.setdefault(fs, {})[feature] = color
    return out


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
