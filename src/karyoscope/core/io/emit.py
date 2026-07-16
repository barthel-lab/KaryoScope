"""Writers for the artifact files a ``karyoscope build`` produces.

A built HKS database directory contains, alongside the ``index/`` files:

* ``hierarchy.tsv`` — the aggregate ``feature_set / child / parent`` edge list
  (:func:`write_hierarchy_tsv`), the source of truth consumed by ``bin`` /
  ``info`` / colors validation.
* ``index/features.<fs>.hierarchy.txt`` — the same edges for one feature set,
  header-less ``child <tab> parent`` (:func:`write_feature_set_hierarchy_txt`),
  passed to ``hks add-feature-set`` and read back by ``hks smooth``.
* ``colors.tsv`` — ``feature_set / feature / color`` (:func:`write_colors_tsv`),
  with colours auto-assigned by :func:`assign_colors` when the user supplies
  none.
* a per-feature-set ``<name> <priority>`` file for priority-mode builds
  (:func:`write_priorities_file`).
* ``manifest.yaml`` (:func:`write_manifest`).

Colour conventions match ``HKS_human_CHM13_v2``: leaves get distinct hues, the
background leaf is grey ``#808080``, and internal / root grouping nodes are the
structural blue-grey ``#B0C4DE``.
"""

from __future__ import annotations

import colorsys
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import yaml

#: Colour for internal / root grouping nodes (matches v2's convention).
STRUCTURAL_COLOR = "#B0C4DE"
#: Colour for a background / gap-fill leaf (matches v2's ``nonacrocentric``).
BACKGROUND_COLOR = "#808080"


def auto_palette(n: int) -> list[str]:
    """Return ``n`` visually distinct hex colours, evenly spaced around the wheel.

    Hues are spread uniformly; lightness/saturation are fixed at mid values so
    colours read clearly against a white background. Deterministic — the same
    ``n`` always yields the same list, keeping rebuilt databases byte-stable.
    """
    if n <= 0:
        return []
    colors: list[str] = []
    for i in range(n):
        hue = i / n
        # Alternate lightness slightly so neighbours with similar hue still differ.
        light = 0.55 if i % 2 == 0 else 0.45
        r, g, b = colorsys.hls_to_rgb(hue, light, 0.6)
        colors.append(f"#{round(r * 255):02X}{round(g * 255):02X}{round(b * 255):02X}")
    return colors


def assign_colors(
    *,
    nodes: Iterable[str],
    leaves: Sequence[str],
    background: str | None = None,
    provided: Mapping[str, str] | None = None,
) -> dict[str, str]:
    """Assign a colour to every node in one feature set.

    Rules, in precedence order:

    1. Any node present in ``provided`` keeps that colour.
    2. The ``background`` leaf (if any) is :data:`BACKGROUND_COLOR`.
    3. Remaining ``leaves`` receive distinct auto-palette colours (in the order
       given, so callers control assignment stability).
    4. Every other node (internal / root grouping node) is
       :data:`STRUCTURAL_COLOR`.

    Returns ``{node: hex}`` covering all of ``nodes``.
    """
    provided = dict(provided or {})
    node_set = set(nodes)

    auto_leaves = [leaf for leaf in leaves if leaf not in provided and leaf != background]
    palette = auto_palette(len(auto_leaves))
    leaf_colors = dict(zip(auto_leaves, palette, strict=True))

    out: dict[str, str] = {}
    for node in node_set:
        if node in provided:
            out[node] = provided[node]
        elif background is not None and node == background:
            out[node] = BACKGROUND_COLOR
        elif node in leaf_colors:
            out[node] = leaf_colors[node]
        else:
            out[node] = STRUCTURAL_COLOR
    return out


def write_hierarchy_tsv(path: Path, rows: Sequence[tuple[str, str, str]]) -> None:
    """Write the aggregate ``feature_set <tab> child <tab> parent`` file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write("feature_set\tchild\tparent\n")
        for feature_set, child, parent in rows:
            out.write(f"{feature_set}\t{child}\t{parent}\n")


def write_feature_set_hierarchy_txt(path: Path, edges: Sequence[tuple[str, str]]) -> None:
    """Write a header-less ``child <tab> parent`` edge list for one feature set."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        for child, parent in edges:
            out.write(f"{child}\t{parent}\n")


def write_colors_tsv(path: Path, rows: Sequence[tuple[str, str, str]]) -> None:
    """Write the ``feature_set <tab> feature <tab> color`` file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write("feature_set\tfeature\tcolor\n")
        for feature_set, feature, color in rows:
            out.write(f"{feature_set}\t{feature}\t{color}\n")


def write_priorities_file(path: Path, priorities: Mapping[str, int]) -> None:
    """Write a whitespace-separated ``<name> <priority>`` file for HKS.

    One line per node; lower value = higher priority. Nodes omitted by HKS
    default to ``0``, but we write every node explicitly for reproducibility.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        for name, prio in priorities.items():
            out.write(f"{name} {prio}\n")


def build_manifest_dict(
    *,
    db_id: str,
    version: str,
    karyoscope_min_version: str,
    basename: str,
    s: int,
    feature_sets: Sequence[str],
    kmer_type: str = "fixed",
    hierarchy: str = "hierarchy.tsv",
    colors: str = "colors.tsv",
    roles: Mapping[str, str] | None = None,
    smoothing: Mapping[str, object] | None = None,
) -> dict:
    """Assemble the manifest mapping for a built HKS database.

    ``features.tsv`` is intentionally omitted: the HKS backend reads label names
    from the index and never consults it (see :mod:`karyoscope.manifest`).

    ``kmer_type`` is ``"fixed"`` (queried only at k=s) or ``"variable"`` (the
    index supports any k<=s); ``size`` is the default query k and ``max_size``
    is s in both cases.
    """
    manifest: dict[str, object] = {
        "id": db_id,
        "version": version,
        "karyoscope_min_version": karyoscope_min_version,
        "index": {"type": "hks", "basename": basename},
        "hierarchy": hierarchy,
        "colors": colors,
        "kmer": {"size": s, "type": kmer_type, "max_size": s},
        "feature_sets": list(feature_sets),
    }
    if roles:
        manifest["roles"] = dict(roles)
    if smoothing:
        manifest["smoothing"] = dict(smoothing)
    return manifest


def write_manifest(path: Path, manifest: Mapping[str, object]) -> None:
    """Write ``manifest.yaml`` preserving key order (no alphabetical sort)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        yaml.safe_dump(dict(manifest), out, sort_keys=False, default_flow_style=False)
