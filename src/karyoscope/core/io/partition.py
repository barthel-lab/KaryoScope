"""BED partition helpers for the ``karyoscope build`` pipeline.

Turns a per-feature-set BED annotation (4th column = leaf label) into the
per-feature FASTA inputs that ``hks add-feature-set`` consumes, plus the
optional pieces that surround it:

* :func:`read_fai` — chromosome lengths from a samtools ``.fai`` index.
* :func:`compute_background_intervals` — the complement of the annotated
  regions, so gaps become a named *background* feature (e.g. ``nonrepeat``)
  rather than HKS's ``none`` miss sentinel. Union-merges the covered regions
  first, so it stays correct even when features overlap.
* :func:`flatten_by_priority` — optional pre-flattening of overlapping regions
  to a single highest-priority label per base (the old KaryoScope behaviour;
  usually unnecessary since HKS resolves overlaps per k-mer).
* :func:`slice_features_to_fastas` — extract each region's sequence into a
  per-label FASTA, streaming one contig at a time so peak memory is O(one
  contig) rather than O(genome). Each region is extended by ``k-1`` bp at the
  end so every k-mer that *starts* inside the region is captured.

Everything here is pure Python (no bedtools/pysam/samtools), matching the
repo's streaming FASTA conventions in :mod:`karyoscope.core.io.fasta`.
"""

from __future__ import annotations

import logging
from collections import Counter, defaultdict
from pathlib import Path

from karyoscope.core.io.fasta import _open_in
from karyoscope.exceptions import KaryoscopeError

logger = logging.getLogger(__name__)

#: A single BED region: ``(chrom, start, end, label)`` with 0-based, half-open
#: coordinates.
BedInterval = tuple[str, int, int, str]


class PartitionError(KaryoscopeError):
    """Problems reading or transforming a BED partition for ``build``."""


def parse_bed(path: Path) -> list[BedInterval]:
    """Read a BED file into a list of ``(chrom, start, end, label)`` tuples.

    Requires at least four columns; the 4th is the feature label. Comment
    (``#``) and blank lines are skipped. Supports plain or ``.gz`` input.

    Raises
    ------
    PartitionError
        On missing file, a row with fewer than four columns, or non-integer
        coordinates.
    """
    if not path.is_file():
        raise PartitionError(f"BED file not found: {path}")

    intervals: list[BedInterval] = []
    with _open_in(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise PartitionError(
                    f"{path}:{lineno}: expected at least 4 tab-separated columns "
                    f"(chrom, start, end, label), got {len(fields)}: {line!r}"
                )
            chrom, start_s, end_s, label = fields[0], fields[1], fields[2], fields[3]
            try:
                start, end = int(start_s), int(end_s)
            except ValueError as e:
                raise PartitionError(
                    f"{path}:{lineno}: start/end must be integers: {line!r}"
                ) from e
            if end < start:
                raise PartitionError(f"{path}:{lineno}: end ({end}) < start ({start}): {line!r}")
            if label == "":
                raise PartitionError(f"{path}:{lineno}: empty label column: {line!r}")
            intervals.append((chrom, start, end, label))
    return intervals


def read_fai(path: Path) -> dict[str, int]:
    """Read a samtools ``.fai`` index into ``{chrom: length}``.

    ``.fai`` format: ``name  length  offset  linebases  linewidth``.

    Raises
    ------
    PartitionError
        On missing file or a malformed length column.
    """
    if not path.is_file():
        raise PartitionError(
            f".fai index not found: {path}. Build it with `samtools faidx <fasta>`."
        )
    sizes: dict[str, int] = {}
    with _open_in(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise PartitionError(f"{path}:{lineno}: malformed .fai line: {line!r}")
            try:
                sizes[parts[0]] = int(parts[1])
            except ValueError as e:
                raise PartitionError(f"{path}:{lineno}: non-integer length: {line!r}") from e
    return sizes


def labels_in(intervals: list[BedInterval]) -> list[str]:
    """Return the distinct labels in ``intervals``, in first-seen order."""
    seen: list[str] = []
    seen_set: set[str] = set()
    for _, _, _, label in intervals:
        if label not in seen_set:
            seen_set.add(label)
            seen.append(label)
    return seen


def _group_by_chrom(intervals: list[BedInterval]) -> dict[str, list[tuple[int, int, str]]]:
    by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for chrom, start, end, label in intervals:
        by_chrom[chrom].append((start, end, label))
    return by_chrom


def compute_background_intervals(
    intervals: list[BedInterval],
    fai: dict[str, int],
    label: str,
) -> list[BedInterval]:
    """Return background regions covering every base no feature annotates.

    For each chromosome in ``fai``, the union of the (possibly overlapping)
    annotated intervals is merged, and the complement within ``[0, length)`` is
    emitted as ``label``. Chromosomes absent from ``fai`` are ignored (a warning
    is logged if they carried annotations, since their length is unknown).

    The result is genome-sorted by ``fai`` order. Feature intervals themselves
    are *not* returned — only the gap fill.
    """
    by_chrom = _group_by_chrom(intervals)

    unknown = sorted(set(by_chrom) - set(fai))
    if unknown:
        logger.warning(
            "%d annotated chromosome(s) absent from the .fai index and skipped "
            "for background fill: %s",
            len(unknown),
            unknown if len(unknown) <= 5 else [*unknown[:5], "..."],
        )

    background: list[BedInterval] = []
    for chrom, length in fai.items():
        covered = sorted((s, e) for s, e, _ in by_chrom.get(chrom, []))
        pos = 0
        for start, end in covered:
            start = max(0, min(start, length))
            end = max(0, min(end, length))
            if start > pos:
                background.append((chrom, pos, start, label))
            pos = max(pos, end)
        if pos < length:
            background.append((chrom, pos, length, label))
    return background


def flatten_by_priority(
    intervals: list[BedInterval],
    order: list[str],
) -> list[BedInterval]:
    """Flatten overlapping regions to one highest-priority label per base.

    ``order`` lists labels highest-priority first. Where regions overlap, the
    base is assigned the label appearing earliest in ``order``; labels absent
    from ``order`` rank last, ties broken alphabetically (matching the archive's
    ``priority_merge.py``). Adjacent segments with the same label are merged.

    This reproduces the old non-overlapping-partition behaviour. It is usually
    unnecessary under HKS, which resolves overlaps per k-mer via priority-aware
    LCA; prefer ``--feature-priorities`` at ``add-feature-set`` time instead.
    """
    rank = {label: i for i, label in enumerate(order)}

    def best(active: Counter[str]) -> str:
        return min(active, key=lambda lbl: (rank.get(lbl, len(order)), lbl))

    result: list[BedInterval] = []
    for chrom, chrom_intervals in _group_by_chrom(intervals).items():
        # Sweep the coordinate line; a Counter of active labels handles overlaps.
        # Each elementary segment [coords[i], coords[i+1]) is assigned the best
        # currently-active label, then adjacent equal-label segments are merged.
        events: dict[int, list[tuple[str, int]]] = defaultdict(list)  # pos -> [(label, +1/-1)]
        for start, end, label in chrom_intervals:
            if end <= start:
                continue
            events[start].append((label, 1))
            events[end].append((label, -1))
        coords = sorted(events)
        active: Counter[str] = Counter()
        seg_start: int | None = None
        seg_label: str | None = None
        for i in range(len(coords) - 1):
            pos = coords[i]
            for label, delta in events[pos]:
                active[label] += delta
                if active[label] <= 0:
                    del active[label]
            if not active:
                if seg_label is not None:
                    result.append((chrom, seg_start, pos, seg_label))
                    seg_start, seg_label = None, None
                continue
            label = best(active)
            if seg_label is None:
                seg_start, seg_label = pos, label
            elif label != seg_label:
                result.append((chrom, seg_start, pos, seg_label))
                seg_start, seg_label = pos, label
        if seg_label is not None:
            result.append((chrom, seg_start, coords[-1], seg_label))
    return result


def slice_features_to_fastas(
    fasta_path: Path,
    intervals: list[BedInterval],
    k: int,
    outdir: Path,
) -> dict[str, Path]:
    """Write one FASTA per label, extracting each region's sequence (+ ``k-1``).

    Streams ``fasta_path`` one contig at a time (plain or ``.gz``), so peak
    memory is O(largest contig). Each region ``[start, end)`` is written as the
    subsequence ``[start, min(end + k - 1, contig_len))`` — the ``k-1`` tail
    ensures every k-mer whose start lies in the region is present, matching the
    archive's ``subset_features_to_fastas.py``.

    Regions on contigs absent from the FASTA are skipped. Records are named
    ``>{chrom}:{start}-{end}_{label}``.

    Returns
    -------
    dict[str, Path]
        ``{label: fasta_path}`` for every label that produced at least one
        record, in first-seen order.
    """
    if k < 1:
        raise PartitionError(f"k must be >= 1, got {k}")
    outdir.mkdir(parents=True, exist_ok=True)
    extension = k - 1
    by_chrom = _group_by_chrom(intervals)

    handles: dict[str, object] = {}
    paths: dict[str, Path] = {}

    def handle_for(label: str):
        h = handles.get(label)
        if h is None:
            path = outdir / f"{label}.fasta"
            paths[label] = path
            h = path.open("w")
            handles[label] = h
        return h

    def flush_contig(name: str | None, seq_parts: list[str]) -> None:
        if name is None or name not in by_chrom:
            return
        seq = "".join(seq_parts)
        length = len(seq)
        for start, end, label in by_chrom[name]:
            s = max(0, min(start, length))
            e = max(s, min(end + extension, length))
            if e <= s:
                continue
            handle_for(label).write(f">{name}:{start}-{end}_{label}\n{seq[s:e]}\n")

    try:
        current_name: str | None = None
        current_seq: list[str] = []
        with _open_in(fasta_path) as fh:
            for raw in fh:
                if raw.startswith(">"):
                    flush_contig(current_name, current_seq)
                    head = raw[1:].strip()
                    current_name = head.split()[0] if head else ""
                    current_seq = []
                else:
                    current_seq.append(raw.strip())
            flush_contig(current_name, current_seq)
    finally:
        for h in handles.values():
            h.close()

    # Preserve first-seen label order for a stable feature-file-list.
    ordered = {label: paths[label] for label in labels_in(intervals) if label in paths}
    return ordered
