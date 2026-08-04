"""Structural feature sets: whole-sequence ``chromosome``, and satellite ``region``.

:func:`from_fai` is the trivial one — every sequence becomes a single record
labelled with its own name.

:func:`from_satellite` is the one converter that legitimately tiles, because the
tiling is *semantic* rather than a gap-fill: the arms either side of the
centromere have to be told apart (``p_arm``/``q_arm``), and only the satellite
annotation says where the centromere is. ``build``'s generic gap-fill cannot do
that — it has exactly one background label.
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
