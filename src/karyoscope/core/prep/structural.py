"""Structural feature sets: whole-sequence ``chromosome``, and centromeric ``region``.

:func:`from_fai` is the trivial one — every sequence becomes a single record
labelled with its own name.

:func:`from_censat` (a CenSat annotation, as for human) and
:func:`from_satellite` (a bare satellite-monomer catalog, as for the Arabidopsis
CEN180 set) both build a ``region`` set, and are the two converters that
legitimately tile. The tiling is *semantic* rather than a gap-fill: the arms
either side of the centromere have to be told apart (``p_arm``/``q_arm``), only
the annotation knows where the centromere is, and ``build``'s gap-fill has
exactly one label where two are needed.

They differ in what the input gives them. CenSat already names its features, so
the work is normalising labels and locating the centromere from the ``ct``
transitions; a monomer catalog names nothing, so the work is merging monomers
into bands and finding the densest cluster.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from karyoscope.core.prep.common import (
    Edge,
    PrepError,
    PrepResult,
    Record,
    SeqidRewriter,
    coalesce,
    open_text,
    read_fai,
    write_bed,
    write_hierarchy,
    write_priority,
)

# -- chromosome -------------------------------------------------------


def from_fai(
    *,
    lengths: Path,
    output: Path,
    hierarchy: Path | None = None,
    name: str = "chromosome",
    rename: SeqidRewriter | None = None,
) -> PrepResult:
    """Emit one whole-length record per sequence, labelled with its own name.

    No hierarchy is written by default: how sequences group (autosome vs sex vs
    organelle, metacentric vs acrocentric, haplotype) is organism-specific
    curation that cannot be read off a ``.fai``. Pass ``hierarchy`` only if you
    intend to fill it in by hand afterwards.
    """
    rename = rename or SeqidRewriter()
    sizes = read_fai(lengths)
    records: list[Record] = [
        (rename(seqid), 0, length, rename(seqid)) for seqid, length in sizes.items()
    ]
    n_records = write_bed(output, records)

    n_edges = 0
    if hierarchy is not None:
        n_edges = write_hierarchy(hierarchy, [(rename(s), "categorized") for s in sizes])

    return PrepResult(
        name=name,
        bed=output,
        n_records=n_records,
        background=None,  # one record per sequence already tiles everything
        hierarchy=hierarchy,
        n_edges=n_edges,
        notes=[
            f"{n_records} sequence(s) from {lengths}",
            "no grouping hierarchy was derived — sequence grouping is organism-specific "
            "curation; write it by hand and pass it to build as 'hierarchy:'",
        ],
    )


# -- satellite / region -----------------------------------------------

INTERIOR = "cen_gap"
P_ARM = "p_arm"
Q_ARM = "q_arm"

# -- CenSat -----------------------------------------------------------

#: ``child parent`` edges of the CenSat v2.1 region tree, as used by the shipped
#: CHM13 v2 ``region`` set. ``p_arm``/``q_arm`` must stay leaves: centromere
#: detection reads anything that is not an arm as the centromere catch-all.
CENSAT_HIERARCHY: list[Edge] = [
    ("centromeric", "categorized"),
    ("aSat", "centromeric"),
    ("alpha_hor", "aSat"),
    ("active_hor", "alpha_hor"),
    ("dhor", "alpha_hor"),
    ("hor", "alpha_hor"),
    ("mixedAlpha", "alpha_hor"),
    ("mon", "aSat"),
    ("bSat", "centromeric"),
    ("cenSat", "centromeric"),
    ("ct", "centromeric"),
    ("gSat", "centromeric"),
    ("HSat", "centromeric"),
    ("HSat1", "HSat"),
    ("HSat1A", "HSat1"),
    ("HSat1B", "HSat1"),
    ("HSat2", "HSat"),
    ("HSat3", "HSat"),
    ("rDNA", "categorized"),
    ("arm", "categorized"),
    (P_ARM, "arm"),
    (Q_ARM, "arm"),
]

#: Priorities for the three top-level branches; everything else is 1. Siblings
#: must be all-equal or all-distinct, which this satisfies.
CENSAT_PRIORITIES = {"centromeric": 1, "rDNA": 2, "arm": 3}

#: Labels that are arms rather than centromeric features, when locating the
#: centromere by feature extent.
_NON_CENTROMERIC = frozenset({P_ARM, Q_ARM})


def censat_label(raw: str) -> str:
    """``gSat(TAR1)`` -> ``gSat``.

    CenSat qualifies each label with the specific arrays it contains, in
    parentheses — useful provenance, but hundreds of distinct values that all
    mean the same feature. The part before the parenthesis is the leaf.
    """
    return raw.split("(", 1)[0]


def _centromere_bounds(features: list[tuple[int, int, str]]) -> tuple[int, int] | None:
    """Centromere extent for one sequence: first ``ct`` start to last ``ct`` end.

    ``ct`` (centromeric transition) brackets the centromere, so it gives the
    boundary directly. Where a sequence has none, fall back to the extent of
    every centromeric feature — which is wider, but only used for deciding
    which arm a gap belongs to.
    """
    ct = [(s, e) for s, e, label in features if label == "ct"]
    if ct:
        return min(s for s, _e in ct), max(e for _s, e in ct)
    other = [(s, e) for s, e, label in features if label not in _NON_CENTROMERIC]
    if other:
        return min(s for s, _e in other), max(e for _s, e in other)
    return None


def from_censat(
    *,
    input_path: Path,
    lengths: Path,
    output: Path,
    hierarchy: Path,
    priority: Path | None = None,
    name: str = "region",
    rename: SeqidRewriter | None = None,
) -> PrepResult:
    """Turn a CenSat annotation into a fully-tiled ``region`` feature set.

    CenSat labels are kept (minus their parenthetical array detail) and every
    remaining base is labelled by which arm it falls on, split at the
    centromere. Like :func:`from_satellite` this tiles, and for the same
    reason: only the annotation knows where the centromere is, and a gap-fill
    has one label where two are needed.

    Scattered pericentromeric satellite remnants far out on an arm keep their
    own labels; only the gaps around them become arm.
    """
    rename = rename or SeqidRewriter()
    sizes = read_fai(lengths)

    by_seqid: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    unknown: set[str] = set()
    with open_text(input_path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            seqid = rename(parts[0])
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            if seqid not in sizes:
                unknown.add(seqid)
                continue
            by_seqid[seqid].append((start, end, censat_label(parts[3])))

    if not by_seqid:
        raise PrepError(
            f"{input_path}: no CenSat records matched a sequence in {lengths} — check "
            "whether the annotation seqids need --rename-prefix or --seqid-map"
        )

    records: list[Record] = []
    bare: list[str] = []
    for seqid, length in sizes.items():
        features = coalesce(
            [(seqid, s, e, label) for s, e, label in sorted(by_seqid.get(seqid, []))]
        )
        if not features:
            bare.append(seqid)
            continue
        bounds = _centromere_bounds([(s, e, label) for _c, s, e, label in features])
        last = 0
        for _c, start, end, label in features:
            if start > last:
                records.append((seqid, last, start, _arm_label(last, start, bounds)))
            records.append((seqid, start, end, label))
            last = end
        if last < length:
            records.append((seqid, last, length, _arm_label(last, length, bounds)))

    n_records = write_bed(output, records)
    n_edges = write_hierarchy(hierarchy, CENSAT_HIERARCHY)
    if priority is not None:
        with priority.open("w") as out:
            for child, parent in CENSAT_HIERARCHY:
                out.write(f"{child}\t{CENSAT_PRIORITIES.get(child, 1)}\t{parent}\n")

    notes = [
        f"{sum(len(v) for v in by_seqid.values()):,} CenSat records over {len(by_seqid)} sequence(s)"
    ]
    if unknown:
        notes.append(
            f"skipped {len(unknown)} annotation seqid(s) absent from {lengths}: "
            + ", ".join(sorted(unknown)[:5])
        )
    if bare:
        notes.append(
            f"{len(bare)} sequence(s) have no CenSat annotation and are left uncovered "
            "(list them under the spec's exclude:): " + ", ".join(bare[:5])
        )

    return PrepResult(
        name=name,
        bed=output,
        n_records=n_records,
        background=None,  # the arm split already tiles every annotated sequence
        hierarchy=hierarchy,
        n_edges=n_edges,
        priority=priority,
        exclude=bare,
        notes=notes,
    )


def _arm_label(start: int, end: int, bounds: tuple[int, int] | None) -> str:
    """Which arm a gap belongs to, given the centromere extent."""
    if bounds is None:
        return P_ARM
    cen_start, cen_end = bounds
    if end <= cen_start:
        return P_ARM
    if start >= cen_end:
        return Q_ARM
    # Inside the centromere: only reachable where CenSat coverage is incomplete.
    return P_ARM


def _hierarchy_for(satellite: str) -> list[Edge]:
    """Arms and satellite as leaves under ``arm``/``centromeric``.

    ``scaffold.get_simple_region`` reads ``p_arm``/``q_arm`` as arm and treats
    everything else as the centromere catch-all, and binning prefers leaf
    labels, so all four of these must be leaves rather than interior nodes.
    """
    return [
        (satellite, "centromeric"),
        (INTERIOR, "centromeric"),
        ("centromeric", "categorized"),
        (P_ARM, "arm"),
        (Q_ARM, "arm"),
        ("arm", "categorized"),
    ]


def merge_spans(spans: list[tuple[int, int]], gap: int) -> list[tuple[int, int]]:
    """Coalesce ``[start, end)`` spans separated by at most ``gap`` bases."""
    merged: list[list[int]] = []
    for start, end in sorted(spans):
        if merged and start - merged[-1][1] <= gap:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def core_cluster(bands: list[tuple[int, int]], cluster_gap: int) -> tuple[int, int]:
    """Return the span of the highest-coverage cluster of satellite bands.

    The centromere is the *densest* cluster, not the min-to-max extent: scattered
    pericentromeric remnants would blow a raw extent out across most of the arm.
    """
    best, best_coverage = (0, 0), -1
    for cluster_start, cluster_end in merge_spans(bands, cluster_gap):
        coverage = sum(e - s for s, e in bands if s >= cluster_start and e <= cluster_end)
        if coverage > best_coverage:
            best, best_coverage = (cluster_start, cluster_end), coverage
    return best


def _fill(
    records: list[Record], seqid: str, start: int, end: int, core: tuple[int, int] | None
) -> None:
    """Label a gap: ``p_arm`` before the core, ``cen_gap`` inside it, ``q_arm`` after."""
    if end <= start:
        return
    if core is None:
        records.append((seqid, start, end, P_ARM))
        return
    core_start, core_end = core
    if start < core_start:
        records.append((seqid, start, min(end, core_start), P_ARM))
    inner_start, inner_end = max(start, core_start), min(end, core_end)
    if inner_start < inner_end:
        records.append((seqid, inner_start, inner_end, INTERIOR))
    if end > core_end:
        records.append((seqid, max(start, core_end), end, Q_ARM))


def from_satellite(
    *,
    input_path: Path,
    lengths: Path,
    output: Path,
    hierarchy: Path,
    priority: Path | None = None,
    name: str = "region",
    satellite: str = "satellite",
    merge_gap: int = 10,
    cluster_gap: int = 500_000,
    rename: SeqidRewriter | None = None,
) -> PrepResult:
    """Turn a centromeric-satellite GFF/BED into a fully-tiled ``region`` set.

    Monomers within ``merge_gap`` bases coalesce into bands — the default bridges
    the 1-2 bp monomer-boundary artefacts tandem-repeat finders emit without
    merging real interior insertions, which run to kilobases. Bands then cluster
    within ``cluster_gap`` to locate the centromere core, and every remaining
    base is labelled by which side of that core it falls on.

    p and q are assigned by coordinate, so this assumes the assembly is oriented
    short-arm-first — the usual convention, but not guaranteed.
    """
    rename = rename or SeqidRewriter()
    sizes = read_fai(lengths)

    monomers: dict[str, list[tuple[int, int]]] = defaultdict(list)
    unknown: set[str] = set()
    is_gff = input_path.suffix.lower().lstrip(".").startswith("gff")
    with open_text(input_path) as fh:
        for raw in fh:
            if raw.startswith("#"):
                if raw.startswith("##FASTA"):
                    break
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            seqid = rename(parts[0])
            try:
                if is_gff:
                    if len(parts) < 5:
                        continue
                    start, end = int(parts[3]) - 1, int(parts[4])
                else:
                    start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue  # header rows and other non-coordinate lines
            if seqid not in sizes:
                unknown.add(seqid)
                continue
            monomers[seqid].append((start, end))

    if not monomers:
        raise PrepError(
            f"{input_path}: no satellite records matched a sequence in {lengths} — check "
            "whether the annotation seqids need --rename-prefix or --seqid-map"
        )

    records: list[Record] = []
    no_core: list[str] = []
    for seqid, length in sizes.items():
        bands = merge_spans(monomers.get(seqid, []), merge_gap)
        if not bands:
            # No satellite at all (an organelle, a small contig): all one arm.
            no_core.append(seqid)
            _fill(records, seqid, 0, length, None)
            continue
        core = core_cluster(bands, cluster_gap)
        last = 0
        for band_start, band_end in bands:
            _fill(records, seqid, last, band_start, core)
            records.append((seqid, band_start, band_end, satellite))
            last = band_end
        _fill(records, seqid, last, length, core)

    edges = _hierarchy_for(satellite)
    n_records = write_bed(output, coalesce(records))
    n_edges = write_hierarchy(hierarchy, edges)
    if priority is not None:
        write_priority(priority, edges)

    notes = [
        f"{sum(len(v) for v in monomers.values()):,} monomers merged into "
        f"{sum(len(merge_spans(v, merge_gap)) for v in monomers.values()):,} "
        f"{satellite} bands across {len(monomers)} sequence(s)"
    ]
    if no_core:
        notes.append(
            f"{len(no_core)} sequence(s) have no {satellite} and are labelled {P_ARM} "
            "end to end: " + ", ".join(no_core[:5])
        )
    if unknown:
        notes.append(
            f"skipped {len(unknown)} annotation seqid(s) absent from {lengths}: "
            + ", ".join(sorted(unknown)[:5])
        )

    return PrepResult(
        name=name,
        bed=output,
        n_records=n_records,
        background=None,  # the arm/core split already tiles every base
        hierarchy=hierarchy,
        n_edges=n_edges,
        priority=priority,
        notes=notes,
    )
