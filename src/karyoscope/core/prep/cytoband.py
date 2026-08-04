"""A UCSC cytoband table -> a ``cytoband`` feature set.

Handles both shapes UCSC ships: the golden-path ``cytoBand.txt(.gz)``
(``chrom start end band stain``) used for hg19/hg38, and the 6-column
``cytoBandMapped`` BED used for CHM13. Both carry the Giemsa stain in column 5,
which is what the colours and the legend grouping are built from.

Band labels keep the chromosome (``chr1`` + ``p36.33`` -> ``1p36.33``) so a
label is never ambiguous across chromosomes, and the hierarchy nests them three
deep: chromosome -> band group -> band. Because a whole-genome cytoband set has
several hundred leaves, the emitted colours file groups the legend by stain —
the CHM13 set's 833 rendered features collapse to 9 legend rows.

Sequences with no cytogenetic banding (``_alt``/``_random``/``_fix``
scaffolds, ``chrUn_*`` contigs, and the mitochondrion) are not labelled. They
are reported for the build spec's top-level ``exclude:`` rather than given a
placeholder label, so no feature set claims them and they read as ``none``.
"""

from __future__ import annotations

import re
from collections import OrderedDict
from pathlib import Path

from karyoscope.core.prep.common import (
    ColorRow,
    Edge,
    PrepError,
    PrepResult,
    Record,
    SeqidRewriter,
    open_text,
    read_fai,
    write_bed,
    write_colors,
    write_hierarchy,
    write_priority,
)

CATEGORIZED = "categorized"

#: Sequences that carry cytogenetic banding. Deliberately strict: everything a
#: reference FASTA adds beyond the primary chromosomes (alt loci, patches,
#: unplaced contigs, organelles) has no banding to assign.
DEFAULT_PRIMARY_PATTERN = r"^chr([0-9]+|X|Y)$"

#: Colour for interior nodes (the root, chromosomes, band groups). They never
#: render — smoothed output carries leaves — but ``build`` colours every node.
INTERIOR_COLOR = "#B0C4DE"
INTERIOR_GROUP = "mixed"

#: Giemsa stain -> colour, as shipped in ``HKS_human_CHM13_cytoband``. Written
#: in intensity order so the legend, which orders groups by first appearance,
#: reads as the stain progression rather than alphabetically.
STAIN_COLORS: dict[str, str] = {
    "gneg": "#E0E0E0",
    "gpos25": "#A8A8A8",
    "gpos33": "#8C8C8C",  # interpolated; not present in the CHM13/hg38 tables
    "gpos50": "#707070",
    "gpos66": "#545454",  # interpolated; not present in the CHM13/hg38 tables
    "gpos75": "#383838",
    "gpos100": "#000000",
    "acen": "#E41A1C",
    "gvar": "#4292C6",
    "stalk": "#41AB5D",
}


def short_name(seqid: str) -> str:
    """``chr1`` -> ``1`` — the label prefix. The seqid itself is left alone."""
    return seqid[3:] if seqid.startswith("chr") else seqid


def full_band_name(seqid: str, band: str) -> str:
    """Qualify a band with its chromosome, unless it already is.

    UCSC ships the band column bare (``p36.33``) in both the 5-column
    golden-path table and the 6-column ``cytoBandMapped`` BED — the latter
    carries the qualified name in a *separate* column, which we ignore and
    rebuild, so one code path serves both. But a redistributed or hand-cut file
    may put the qualified name in column 4 directly, and blindly prefixing that
    would give ``11p36.33``.

    Bare cytogenetic bands always begin with the arm, ``p`` or ``q``, so a band
    that instead begins with its own chromosome name is already qualified.
    """
    prefix = short_name(seqid)
    if band.startswith(prefix) and band[len(prefix) :].startswith(("p", "q")):
        return band
    return f"{prefix}{band}"


def band_group(fullname: str) -> str | None:
    """Group node for a band, or ``None`` when it has no sub-band.

    ``1p36.33`` -> ``1p36``; ``1p33`` -> ``None`` (it hangs off its chromosome).
    """
    return fullname.split(".", 1)[0] if "." in fullname else None


def read_bands(path: Path) -> OrderedDict[str, list[tuple[int, int, str, str]]]:
    """Read a cytoband table into ``{seqid: [(start, end, band, stain), ...]}``.

    The 6-column BED's column 6 full name is ignored and rebuilt from chromosome
    plus band, so both input shapes take the same code path and produce
    identical labels.
    """
    by_seqid: OrderedDict[str, list[tuple[int, int, str, str]]] = OrderedDict()
    with open_text(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                raise PrepError(f"{path}:{lineno}: expected at least 4 columns, got {len(parts)}")
            band = parts[3].strip()
            if not band:  # e.g. chrM's single unnamed band
                continue
            stain = parts[4].strip() if len(parts) >= 5 else ""
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError as e:
                raise PrepError(f"{path}:{lineno}: non-integer coordinates") from e
            by_seqid.setdefault(parts[0], []).append((start, end, band, stain))
    if not by_seqid:
        raise PrepError(f"{path}: no cytoband rows found")
    return by_seqid


def hierarchy_edges(chrom_bands: OrderedDict[str, list[tuple[int, int, str, str]]]) -> list[Edge]:
    """Ordered edges: chromosome -> root, band group -> chromosome, band -> group."""
    edges: list[Edge] = []
    for seqid, bands in chrom_bands.items():
        node = short_name(seqid)
        edges.append((node, CATEGORIZED))
        seen_groups: set[str] = set()
        for _s, _e, fullname, _stain in bands:
            group = band_group(fullname)
            if group is None:
                edges.append((fullname, node))
            else:
                if group not in seen_groups:
                    edges.append((group, node))
                    seen_groups.add(group)
                edges.append((fullname, group))
    return edges


def verify(
    chrom_bands: OrderedDict[str, list[tuple[int, int, str, str]]], edges: list[Edge]
) -> None:
    """Check every band label is a hierarchy node with a clean path to the root.

    Cheap to do here and much cheaper than discovering it part-way through a
    multi-hour ``build``.
    """
    children = [c for c, _p in edges]
    duplicates = {c for c in children if children.count(c) > 1}
    if duplicates:
        raise PrepError(f"duplicate hierarchy node(s): {sorted(duplicates)}")
    parent = dict(edges)
    leaves = {full for bands in chrom_bands.values() for _s, _e, full, _st in bands}
    collisions = leaves & {g for g in parent if g not in leaves and g != CATEGORIZED}
    if collisions:
        raise PrepError(f"label is both a band and a group: {sorted(collisions)}")
    for leaf in leaves:
        node, depth = leaf, 0
        while node in parent and parent[node] != CATEGORIZED:
            node, depth = parent[node], depth + 1
            if depth > 5:
                raise PrepError(f"parent chain too deep for band {leaf!r}")
        if node not in parent:
            raise PrepError(f"band {leaf!r} has no path to {CATEGORIZED!r}")


def from_ucsc(
    *,
    input_path: Path,
    lengths: Path,
    output: Path,
    hierarchy: Path,
    colors: Path | None = None,
    priority: Path | None = None,
    name: str = "cytoband",
    primary_pattern: str = DEFAULT_PRIMARY_PATTERN,
    rename: SeqidRewriter | None = None,
) -> PrepResult:
    """Convert a UCSC cytoband table into a ``cytoband`` feature set."""
    rename = rename or SeqidRewriter()
    try:
        primary = re.compile(primary_pattern)
    except re.error as e:
        raise PrepError(f"--primary-pattern is not a valid regex: {e}") from e

    sizes = read_fai(lengths)
    raw_bands = read_bands(input_path)
    bands_by_seqid = {rename(seqid): bands for seqid, bands in raw_bands.items()}

    chrom_bands: OrderedDict[str, list[tuple[int, int, str, str]]] = OrderedDict()
    excluded: list[str] = []
    unbanded: list[str] = []
    for seqid in sizes:
        if not primary.match(seqid):
            excluded.append(seqid)
            continue
        if seqid not in bands_by_seqid:
            unbanded.append(seqid)
            continue
        chrom_bands[seqid] = [
            (start, end, full_band_name(seqid, band), stain)
            for start, end, band, stain in bands_by_seqid[seqid]
        ]
    if not chrom_bands:
        raise PrepError(
            f"no sequences kept — {primary_pattern!r} matched none of the {len(sizes)} "
            f"seqids in {lengths}; widen --primary-pattern or rename with --rename-prefix"
        )

    edges = hierarchy_edges(chrom_bands)
    verify(chrom_bands, edges)

    records: list[Record] = [
        (seqid, start, end, fullname)
        for seqid, bands in chrom_bands.items()
        for start, end, fullname, _stain in bands
    ]
    n_records = write_bed(output, records)
    n_edges = write_hierarchy(hierarchy, edges)
    if priority is not None:
        write_priority(priority, edges)

    notes: list[str] = []
    n_colors = 0
    if colors is not None:
        stain_of = {full: stain for bands in chrom_bands.values() for _s, _e, full, stain in bands}
        unknown_stains: set[str] = set()
        rows: list[ColorRow] = [(CATEGORIZED, INTERIOR_COLOR, INTERIOR_GROUP)]
        for child, _parent in edges:
            stain = stain_of.get(child)
            if stain is None:  # a chromosome or band-group node
                rows.append((child, INTERIOR_COLOR, INTERIOR_GROUP))
            elif stain in STAIN_COLORS:
                rows.append((child, STAIN_COLORS[stain], stain))
            else:
                unknown_stains.add(stain or "(blank)")
                rows.append((child, INTERIOR_COLOR, stain or INTERIOR_GROUP))
        n_colors = write_colors(colors, name, rows)
        if unknown_stains:
            notes.append(
                f"no colour known for stain(s) {', '.join(sorted(unknown_stains))} — "
                f"they were left at {INTERIOR_COLOR}; edit {colors} to give them one"
            )

    kept = len(chrom_bands)
    notes.insert(0, f"kept {kept} banded sequence(s), {n_records:,} bands")
    if unbanded:
        notes.append(
            f"{len(unbanded)} sequence(s) matched --primary-pattern but have no bands: "
            + ", ".join(unbanded[:5])
        )

    return PrepResult(
        name=name,
        bed=output,
        n_records=n_records,
        background=None,  # UCSC bands tile each banded chromosome end to end
        hierarchy=hierarchy,
        n_edges=n_edges,
        priority=priority,
        colors=colors,
        n_colors=n_colors,
        exclude=excluded + unbanded,
        notes=notes,
    )
