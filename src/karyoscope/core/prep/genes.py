"""A gene annotation (GFF3 or GTF) -> an ``exon``/``intron``/``intergenic`` set.

Exons are read directly; introns are derived per transcript as the gaps between
its consecutive exons; everything else is intergenic. Where annotations of
different transcripts disagree about a base, the more specific label wins —
``exon`` over ``intron`` over ``intergenic`` — so overlapping transcripts and
alternative splicing never produce a base labelled twice.

The result tiles every sequence in the ``.fai``, so the set needs no gap-fill.
"""

from __future__ import annotations

import re
from collections import defaultdict
from itertools import pairwise
from pathlib import Path

from karyoscope.core.prep.common import (
    ColorRow,
    Edge,
    PrepError,
    PrepResult,
    Record,
    SeqidRewriter,
    iter_gff,
    read_fai,
    write_bed,
    write_colors,
    write_hierarchy,
)

EXON = "exon"
INTRON = "intron"
INTERGENIC = "intergenic"

#: All three labels are leaves directly under the root — there is no useful
#: intermediate grouping, and this matches the shipped CHM13 v2 ``gene`` set.
GENE_HIERARCHY: list[Edge] = [
    (EXON, "categorized"),
    (INTRON, "categorized"),
    (INTERGENIC, "categorized"),
]

#: The shipped CHM13 v2 ``gene`` palette.
GENE_COLORS: dict[str, str] = {
    EXON: "#090088",
    INTRON: "#198450",
    INTERGENIC: "#E89EB8",
}

_GTF_ATTR = re.compile(r'(\w+)\s+"([^"]*)"')
_GFF3_PARENT = re.compile(r"(?:^|;)\s*Parent=([^;]+)")


def transcript_ids(attributes: str) -> list[str]:
    """Extract the transcript(s) an exon belongs to, from GFF3 or GTF attributes.

    Tries three forms, so no ``--format`` flag is needed to tell dialects apart:
    ``transcript_id=...`` (GFF3 syntax), ``transcript_id "..."`` (GTF syntax),
    and finally ``Parent=`` (canonical GFF3, possibly comma-separated when one
    exon is shared by several transcripts).

    ``transcript_id`` wins over ``Parent`` when both are present, and that
    precedence is load-bearing. In canonical GFF3 an exon's ``Parent`` is its
    mRNA, but Liftoff-style output sets ``Parent`` to the **gene** while giving
    the transcript in ``transcript_id`` — the same file content is published
    both ways. Preferring ``Parent`` there would group every transcript of a
    gene together, and introns derived from a merged exon list span the gaps
    *between* transcripts, inflating intron at intergenic's expense. That is
    the defect this converter exists to avoid.
    """
    attrs = dict(_GTF_ATTR.findall(attributes))
    tid = attrs.get("transcript_id") or _gff3_attr(attributes, "transcript_id")
    if tid:
        return [tid]
    parent = _GFF3_PARENT.search(attributes)
    if parent:
        return [p.strip() for p in parent.group(1).split(",") if p.strip()]
    return []


def _gff3_attr(attributes: str, key: str) -> str | None:
    """Value of a GFF3 ``key=value`` attribute, or ``None``."""
    match = re.search(rf"(?:^|;)\s*{re.escape(key)}=([^;]+)", attributes)
    return match.group(1).strip() if match else None


def _partition(
    length: int,
    exons: list[tuple[int, int]],
    introns: list[tuple[int, int]],
) -> list[tuple[int, int, str]]:
    """Resolve one sequence into non-overlapping labelled spans covering ``[0, length)``.

    A boundary sweep carrying one depth counter per label: at each span the
    highest-precedence non-zero counter names it.
    """
    deltas: dict[int, list[int]] = defaultdict(lambda: [0, 0])
    for column, spans in ((0, exons), (1, introns)):
        for start, end in spans:
            start, end = max(0, start), min(length, end)
            if end <= start:
                continue
            deltas[start][column] += 1
            deltas[end][column] -= 1

    out: list[tuple[int, int, str]] = []
    depth_exon = depth_intron = 0
    previous = 0
    for position in sorted(deltas):
        if position > previous:
            if depth_exon > 0:
                label = EXON
            elif depth_intron > 0:
                label = INTRON
            else:
                label = INTERGENIC
            out.append((previous, position, label))
        d_exon, d_intron = deltas[position]
        depth_exon += d_exon
        depth_intron += d_intron
        previous = position
    if previous < length:
        out.append((previous, length, INTERGENIC))
    return out


def _merge_runs(spans: list[tuple[int, int, str]]) -> list[tuple[int, int, str]]:
    """Join abutting spans that ended up with the same label.

    The sweep emits a span per *boundary*, so a run of overlapping exons yields
    one span per exon edge even though they all label ``exon``. Merging them is
    not cosmetic: it is the difference between one record and dozens for every
    multi-transcript gene.
    """
    merged: list[list] = []
    for start, end, label in spans:
        if merged and merged[-1][2] == label and merged[-1][1] == start:
            merged[-1][1] = end
        else:
            merged.append([start, end, label])
    return [(s, e, lbl) for s, e, lbl in merged]


def from_gff(
    *,
    input_path: Path,
    lengths: Path,
    output: Path,
    hierarchy: Path,
    colors: Path | None = None,
    name: str = "gene",
    feature_type: str = "exon",
    rename: SeqidRewriter | None = None,
) -> PrepResult:
    """Convert a GFF3/GTF gene annotation into a fully-tiled gene feature set."""
    rename = rename or SeqidRewriter()
    sizes = read_fai(lengths)

    exons_by_transcript: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    unknown_seqids: set[str] = set()
    n_exons = 0

    for fields, lineno in iter_gff(input_path):
        if fields[2] != feature_type:
            continue
        seqid = rename(fields[0])
        if seqid not in sizes:
            unknown_seqids.add(seqid)
            continue
        try:
            start, end = int(fields[3]) - 1, int(fields[4])
        except ValueError as e:
            raise PrepError(f"{input_path}:{lineno}: non-integer coordinates") from e
        n_exons += 1
        ids = transcript_ids(fields[8])
        if not ids:
            # No parent to derive introns from, but the exon itself still counts.
            exons_by_transcript[(seqid, f"__orphan__:{lineno}")].append((start, end))
            continue
        for tid in ids:
            exons_by_transcript[(seqid, tid)].append((start, end))

    if n_exons == 0:
        raise PrepError(
            f"{input_path}: no '{feature_type}' features whose sequence is in {lengths} — "
            "check --feature-type and whether the annotation seqids need --rename-prefix "
            "or --seqid-map"
        )

    exons_by_seqid: dict[str, list[tuple[int, int]]] = defaultdict(list)
    introns_by_seqid: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for (seqid, _tid), spans in exons_by_transcript.items():
        exons_by_seqid[seqid].extend(spans)
        for (_s, prev_end), (next_start, _e) in pairwise(sorted(spans)):
            if next_start > prev_end:
                introns_by_seqid[seqid].append((prev_end, next_start))

    records: list[Record] = []
    for seqid in sizes:
        spans = _partition(
            sizes[seqid], exons_by_seqid.get(seqid, []), introns_by_seqid.get(seqid, [])
        )
        for start, end, label in _merge_runs(spans):
            records.append((seqid, start, end, label))

    n_records = write_bed(output, records)
    n_edges = write_hierarchy(hierarchy, GENE_HIERARCHY)
    n_colors = 0
    if colors is not None:
        rows: list[ColorRow] = [(label, GENE_COLORS[label], "") for label, _p in GENE_HIERARCHY]
        n_colors = write_colors(colors, name, rows)

    notes = [f"{n_exons:,} '{feature_type}' features over {len(exons_by_seqid)} sequence(s)"]
    if unknown_seqids:
        shown = ", ".join(sorted(unknown_seqids)[:5])
        more = f" (+{len(unknown_seqids) - 5} more)" if len(unknown_seqids) > 5 else ""
        notes.append(
            f"skipped {len(unknown_seqids)} annotation seqid(s) absent from {lengths}: {shown}{more}"
        )
    bare = sorted(s for s in sizes if s not in exons_by_seqid)
    if bare:
        shown = ", ".join(bare[:5])
        more = f" (+{len(bare) - 5} more)" if len(bare) > 5 else ""
        notes.append(
            f"{len(bare)} sequence(s) have no annotation and are all intergenic: {shown}{more}"
        )

    return PrepResult(
        name=name,
        bed=output,
        n_records=n_records,
        background=None,  # exon/intron/intergenic already tiles every base
        hierarchy=hierarchy,
        n_edges=n_edges,
        colors=colors,
        n_colors=n_colors,
        notes=notes,
    )
