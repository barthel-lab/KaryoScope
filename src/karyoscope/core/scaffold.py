"""Classify, orient, and order assembly contigs into canonical scaffolds.

This module is a faithful port of the archive's ``scaffold_stats.py``
(the newer version, post-chrY-centroid fix), with two structural
differences that the Stage 5d design discussion locked in:

1. **Leaf detection is structural**, not string-suffix-based. The
   archive's ``_specific`` suffix hack is gone; chromosome leaves come
   from :func:`karyoscope.core.bin.leaves_for` against the chromosome
   feature set's hierarchy, the same source the binner uses.
2. **The output map file is the source of truth.** The archive
   round-tripped contig identity by string-parsing the renamed
   ``<chrom>_<hap>_<contig>[_rc]`` name; we keep that name format for
   human readability but emit an explicit
   :class:`karyoscope.core.io.scaffold_map.MapRow` per contig so
   downstream stages don't have to parse anything.

The :func:`classify_and_orient` entry point takes pre-built 1 Mb
binned BEDs (chromosome and region feature sets) plus telomere flags
and returns a list of :class:`MapRow`. The caller (typically the CLI
orchestrator in :mod:`karyoscope.commands.scaffold`) is responsible
for producing those binned BEDs and for writing the rewritten
per-feature-set scaffolded BEDs once it has the ordering.
"""

from __future__ import annotations

import gzip
import logging
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
from itertools import groupby
from pathlib import Path
from typing import IO

from karyoscope.core.io.scaffold_map import MapRow
from karyoscope.core.io.telo import TeloFlags
from karyoscope.exceptions import ScaffoldError

logger = logging.getLogger(__name__)


#: Default minimum scaffold length (bp). Shorter contigs without
#: telomeres are dropped from output. Matches the archive default.
DEFAULT_MIN_SCAFFOLD_LENGTH = 5_000_000

#: Threshold (bp) above which a simple-region category counts as
#: "present enough to influence the flip decision". Matches the
#: archive's hardcoded 1_000_000.
_LIKE_THRESHOLD = 1_000_000


#: The standard human acrocentric chromosomes. Used as the default
#: when the user doesn't pass ``--acrocentric``; the CLI warns when
#: falling back.
DEFAULT_HUMAN_ACROCENTRICS: frozenset[str] = frozenset(
    {"chr13", "chr14", "chr15", "chr21", "chr22"}
)


#: An interval as stored in the binned BEDs.
Interval = tuple[int, int, str]


# --- per-contig input -----------------------------------------------


@dataclass
class ContigInput:
    """Everything :func:`classify_and_orient` needs about one contig.

    Constructed by the CLI orchestrator after binning each input's
    annotation BEDs and running seqtk telo.
    """

    input_name: str  # hap label assigned by hap_inference
    input_file: str  # FASTA basename (for the map file)
    contig_name: str  # raw contig name in the source FASTA
    length: int  # max(stop) seen in either binned BED
    chromosome_bins: list[Interval]
    region_bins: list[Interval]
    telo: TeloFlags = field(default_factory=lambda: TeloFlags(False, False))


# --- pure helpers ---------------------------------------------------


def get_simple_region(feature_name: str) -> str:
    """Project a region feature name down to a coarse category.

    The categories are ``"p_arm"``, ``"q_arm"``, ``"telomere"``,
    ``"novel"``, ``"arm"`` (the ambiguous internal node), and
    ``"centromere"`` (the catch-all). Case-insensitive substring
    matching against the human convention; documented in the design
    as the human-genome assumption that we will eventually make
    manifest-driven.
    """
    name = feature_name.lower()
    if "p_arm" in name:
        return "p_arm"
    if "q_arm" in name:
        return "q_arm"
    if name == "arm":
        return "arm"
    if "tel" in name:
        return "telomere"
    if name == "novel":
        return "novel"
    return "centromere"


def assign_main_chromosome(
    chromosome_bins: Iterable[Interval],
    chromosome_leaves: set[str],
) -> str | None:
    """Pick the leaf chromosome with the largest weighted overlap.

    Bins whose feature is not a leaf are ignored. Ties broken by the
    natural :func:`chromosome_sort_key` ordering so the result is
    deterministic when two leaves have identical coverage.
    """
    counts: dict[str, int] = defaultdict(int)
    for start, stop, name in chromosome_bins:
        if name in chromosome_leaves:
            counts[name] += stop - start
    if not counts:
        return None
    best_overlap = max(counts.values())
    leaders = sorted(
        (name for name, v in counts.items() if v == best_overlap),
        key=chromosome_sort_key,
    )
    return leaders[0]


def chromosome_sort_key(chrom: str) -> tuple[int, int, str]:
    """Natural sort for chromosome names: numbers, then alphabetic.

    ``chr1, chr2, ..., chr22, chrX, chrY`` for human; works for any
    naming convention that mixes ``chr<int>`` and ``chr<letter>``.
    """
    suffix = chrom.removeprefix("chr")
    if suffix.isdigit():
        return (0, int(suffix), "")
    return (1, 0, suffix)


def find_largest_contiguous_region(
    chromosome_bins: list[Interval],
    main_chromosome: str | None,
    chromosome_leaves: set[str],
) -> tuple[int, int]:
    """Find the longest run of bins compatible with ``main_chromosome``.

    A bin is *compatible* if its feature is the main chromosome
    itself, an internal hierarchy node (anything not in
    ``chromosome_leaves``), ``"novel"``, or ``"categorized"`` (the
    hierarchy root). A run is *contiguous* if consecutive bins abut
    (``prev.end == next.start``).

    Returns ``(start, end)`` of the longest compatible-and-contiguous
    block, in source coordinates. Falls back to the full extent of
    the input when there is no main chromosome or no compatible bin.
    """
    if not main_chromosome or not chromosome_bins:
        total_end = max((stop for _, stop, _ in chromosome_bins), default=0)
        return 0, total_end

    def _compatible(name: str) -> bool:
        return (
            name == main_chromosome
            or name not in chromosome_leaves
            or name == "novel"
            or name == "categorized"
        )

    max_len = 0
    best_start = 0
    best_end = 0
    cur_start = -1
    cur_end = -1
    cur_len = 0
    last_pos = -1

    for start, stop, name in chromosome_bins:
        is_contiguous = (start == last_pos) or (last_pos == -1)
        compatible = _compatible(name)

        if compatible and is_contiguous:
            if cur_start == -1:
                cur_start = start
            cur_end = stop
            cur_len += stop - start
        else:
            if cur_len > max_len:
                max_len = cur_len
                best_start = cur_start
                best_end = cur_end
            if compatible:
                cur_start = start
                cur_end = stop
                cur_len = stop - start
            else:
                cur_start = -1
                cur_end = -1
                cur_len = 0
        last_pos = stop

    if cur_len > max_len:
        max_len = cur_len
        best_start = cur_start
        best_end = cur_end

    if max_len == 0:
        total_end = max(stop for _, stop, _ in chromosome_bins)
        return 0, total_end

    return best_start, best_end


def half_region_totals(
    region_bins: Iterable[Interval],
    region_start: int,
    region_end: int,
) -> dict[str, list[int]]:
    """Split the ``[region_start, region_end)`` window in half and tally bp per category.

    Returns ``{category: [first_half_bp, second_half_bp]}`` for each of
    ``p_arm``, ``q_arm``, ``centromere``, ``telomere``, and ``novel``.
    The ambiguous ``"arm"`` category is silently dropped (the
    flip-decision logic does not consult it).
    """
    out: dict[str, list[int]] = {
        "p_arm": [0, 0],
        "q_arm": [0, 0],
        "centromere": [0, 0],
        "telomere": [0, 0],
        "novel": [0, 0],
    }
    breakpoint = region_start + (region_end - region_start) / 2

    for start, stop, name in region_bins:
        ostart = max(start, region_start)
        oend = min(stop, region_end)
        if oend <= ostart:
            continue
        simple = get_simple_region(name)
        if simple not in out:
            continue
        first_len = max(min(breakpoint, oend) - ostart, 0)
        second_len = max(oend - max(breakpoint, ostart), 0)
        out[simple][0] += int(first_len)
        out[simple][1] += int(second_len)
    return out


def scaffold_region_majority(
    region_bins: Iterable[Interval],
    region_start: int,
    region_end: int,
) -> str:
    """Return the simple-region category that dominates the window.

    Empty input → empty string. The flip-decision ladder uses this as
    a tie-breaker when totals are below the :data:`_LIKE_THRESHOLD`.
    """
    counts: dict[str, int] = defaultdict(int)
    for start, stop, name in region_bins:
        ostart = max(start, region_start)
        oend = min(stop, region_end)
        if oend <= ostart:
            continue
        counts[get_simple_region(name)] += oend - ostart
    if not counts:
        return ""
    return max(counts, key=lambda k: counts[k])


def need_to_flip(
    region_bins: list[Interval],
    region_half_totals: dict[str, list[int]],
    *,
    region_start: int,
    region_end: int,
    scaffold_length: int,
    telo: TeloFlags,
    is_acrocentric: bool,
) -> bool:
    """Decide whether the contig is reversed (q-then-p) and should be flipped.

    Direct port of the archive's ``calculate_need_to_flip`` (the newer
    version with the centroid-based both-tel rule for chrY-like
    cases). The body is a boolean ladder over arm/centromere
    composition and telomere placement; see the inline comments for
    the intuition behind each branch.
    """
    simple_region = scaffold_region_majority(region_bins, region_start, region_end)

    has_start_tel = telo.start
    has_stop_tel = telo.stop

    continuous_start_tel = has_start_tel and (region_start == 0)
    continuous_stop_tel = has_stop_tel and (region_end == scaffold_length)

    first_half_p, second_half_p = region_half_totals["p_arm"]
    first_half_q, second_half_q = region_half_totals["q_arm"]
    p_total = first_half_p + second_half_p
    q_total = first_half_q + second_half_q
    cen_total = sum(region_half_totals["centromere"])

    is_p_like = (p_total >= _LIKE_THRESHOLD) or (simple_region == "p_arm")
    is_q_like = (q_total >= _LIKE_THRESHOLD) or (simple_region == "q_arm")
    is_c_like = (cen_total >= _LIKE_THRESHOLD) or (simple_region == "centromere")

    is_p_or_c_not_q_like = (is_p_like or is_c_like) and (not is_q_like)
    is_q_not_c_not_p_like = is_q_like and (not is_c_like) and (not is_p_like)
    is_q_and_c_not_p_like = is_q_like and is_c_like and (not is_p_like)
    is_p_and_q_like = is_p_like and is_q_like
    is_pure_telomere = (not is_p_like) and (not is_q_like) and (not is_c_like)

    if continuous_start_tel and continuous_stop_tel:
        # Both ends telomere-capped: pick the orientation that puts
        # the p-arm centroid before the q-arm centroid. Robust against
        # heterochromatic blocks (e.g. chrY's Yq12) that throw off
        # the half-total heuristic.
        p_weighted = 0.0
        q_weighted = 0.0
        for start, stop, name in region_bins:
            ostart = max(start, region_start)
            oend = min(stop, region_end)
            if oend <= ostart:
                continue
            midpoint = (ostart + oend) / 2
            length = oend - ostart
            simple = get_simple_region(name)
            if simple == "p_arm":
                p_weighted += midpoint * length
            elif simple == "q_arm":
                q_weighted += midpoint * length
        if p_total > 0:
            p_weighted /= p_total
        if q_total > 0:
            q_weighted /= q_total
        if p_total > 0 and q_total > 0:
            return p_weighted > q_weighted
        return first_half_q > second_half_q

    if continuous_start_tel or continuous_stop_tel:
        if is_pure_telomere and is_acrocentric and continuous_stop_tel:
            return True
        if is_p_or_c_not_q_like and continuous_stop_tel:
            return True
        if is_q_not_c_not_p_like and continuous_start_tel:
            return True
        if is_q_and_c_not_p_like and (first_half_q > second_half_q):
            return True
        if is_p_and_q_like:
            p_score = first_half_p - second_half_p
            q_score = second_half_q - first_half_q
            if (p_score + q_score) < 0:
                return True
        return False

    # No contiguous-region-touching telomere on either end: use
    # combined p/q-score logic.
    p_score = first_half_p - second_half_p
    q_score = second_half_q - first_half_q
    return (p_score + q_score) < 0


def flip_bins(bins: list[Interval], length: int) -> list[Interval]:
    """Mirror BED intervals to ``[length - stop, length - start)``, reversing order."""
    return [(length - stop, length - start, name) for start, stop, name in reversed(bins)]


def category_index(
    *,
    p_total: int,
    cen_total: int,
    q_total: int,
    has_start_tel: bool,
    has_stop_tel: bool,
) -> int:
    """The 8-bucket category ordering used to sort contigs within a (chrom, hap) cell.

    0 means "starts with a telomere", 7 means "ends with a telomere",
    and 1-6 are the middle buckets keyed off relative p/c/q content.
    Lower indices print first, so the karyotype reads top-down from
    p-arm-ish toward q-arm-ish.
    """
    if has_start_tel:
        return 0
    if has_stop_tel:
        return 7
    if p_total > q_total:
        return 1 if cen_total == 0 else 2
    if p_total < q_total:
        if cen_total == 0:
            return 6
        if p_total > 0:
            return 4
        return 5
    return 3


# --- orchestration --------------------------------------------------


def _stats_string(
    bins: list[Interval],
    region_start: int,
    region_end: int,
    has_start_tel: bool,
    has_stop_tel: bool,
) -> str:
    """Build the TPCQT-style summary string for one contig.

    Walks the oriented region bins in coordinate order. Adjacent
    same-letter runs are collapsed (so a long centromere prints as
    one ``C``, not 30).
    """
    chars: list[str] = []
    if has_start_tel:
        chars.append("T")
    for start, stop, name in bins:
        ostart = max(start, region_start)
        oend = min(stop, region_end)
        if oend <= ostart:
            continue
        simple = get_simple_region(name)
        if simple == "p_arm":
            chars.append("P")
        elif simple == "centromere":
            chars.append("C")
        elif simple == "q_arm":
            chars.append("Q")
    if not chars or chars == ["T"]:
        chars.append("-")
    if has_stop_tel:
        chars.append("T")
    # Collapse adjacent duplicates.
    collapsed: list[str] = []
    for ch, _ in groupby(chars):
        collapsed.append(ch)
    return "".join(collapsed)


@dataclass
class _OrientedContig:
    """Intermediate per-contig state after orientation, before ordering."""

    contig: ContigInput
    chromosome: str
    region_start: int
    region_end: int
    flipped: bool
    flipped_region_bins: list[Interval]
    # Telomere flags AFTER applying the orientation update.
    oriented_telo: TeloFlags
    half_totals: dict[str, list[int]]


def _orient(
    contig: ContigInput,
    main_chromosome: str,
    chromosome_leaves: set[str],
    acrocentrics: set[str],
) -> _OrientedContig:
    region_start, region_end = find_largest_contiguous_region(
        contig.chromosome_bins, main_chromosome, chromosome_leaves
    )
    halfs = half_region_totals(contig.region_bins, region_start, region_end)
    flipped = need_to_flip(
        contig.region_bins,
        halfs,
        region_start=region_start,
        region_end=region_end,
        scaffold_length=contig.length,
        telo=contig.telo,
        is_acrocentric=main_chromosome in acrocentrics,
    )

    if flipped:
        # Update telomere flags: swap start/stop unless both are set.
        if contig.telo.start and contig.telo.stop:
            oriented_telo = contig.telo
        else:
            oriented_telo = TeloFlags(start=contig.telo.stop, stop=contig.telo.start)
        flipped_bins = flip_bins(contig.region_bins, contig.length)
        # The region window also flips for stats-string emission.
        flipped_region_start = contig.length - region_end
        flipped_region_end = contig.length - region_start
        # Recompute half-totals on the flipped data so downstream
        # consumers (category_index, stats string) see the oriented view.
        halfs = half_region_totals(flipped_bins, flipped_region_start, flipped_region_end)
        return _OrientedContig(
            contig=contig,
            chromosome=main_chromosome,
            region_start=flipped_region_start,
            region_end=flipped_region_end,
            flipped=True,
            flipped_region_bins=flipped_bins,
            oriented_telo=oriented_telo,
            half_totals=halfs,
        )

    return _OrientedContig(
        contig=contig,
        chromosome=main_chromosome,
        region_start=region_start,
        region_end=region_end,
        flipped=False,
        flipped_region_bins=contig.region_bins,
        oriented_telo=contig.telo,
        half_totals=halfs,
    )


def classify_and_orient(
    contigs: list[ContigInput],
    *,
    chromosome_leaves: set[str],
    acrocentrics: set[str] = DEFAULT_HUMAN_ACROCENTRICS,
    min_scaffold_length: int = DEFAULT_MIN_SCAFFOLD_LENGTH,
) -> list[MapRow]:
    """Classify, orient, and order all contigs across all inputs.

    Returns one :class:`MapRow` per kept contig, in the canonical
    chromosome x haplotype x category x descending-length order that
    the downstream karyotype renderer expects.

    Contigs shorter than ``min_scaffold_length`` without a telomere
    are dropped. Contigs with no leaf-chromosome hits in their
    chromosome BED are also dropped (we can't classify them).
    """
    # Filter by length / telomere.
    kept: list[ContigInput] = [
        c for c in contigs if c.length >= min_scaffold_length or c.telo.start or c.telo.stop
    ]

    # Assign each surviving contig to a chromosome.
    main_chrom: dict[int, str | None] = {}
    for i, c in enumerate(kept):
        main_chrom[i] = assign_main_chromosome(c.chromosome_bins, chromosome_leaves)

    # Group by (chromosome, hap). Skip contigs without a chromosome.
    cells: dict[tuple[str, str], list[int]] = defaultdict(list)
    for i, c in enumerate(kept):
        chrom = main_chrom[i]
        if chrom is None:
            logger.info(
                "dropping contig %r: no leaf-chromosome features in its chromosome BED",
                c.contig_name,
            )
            continue
        cells[(chrom, c.input_name)].append(i)

    # Emit rows in chromosome order, then hap order, then category x length within each cell.
    rows: list[MapRow] = []
    for chrom, hap in sorted(cells.keys(), key=lambda k: (chromosome_sort_key(k[0]), k[1])):
        cell_indices = cells[(chrom, hap)]
        oriented = [_orient(kept[i], chrom, chromosome_leaves, acrocentrics) for i in cell_indices]
        # Build sort keys per oriented contig: (category, -length, original_name).
        decorated: list[tuple[tuple[int, int, str], _OrientedContig]] = []
        for o in oriented:
            p_total = sum(o.half_totals["p_arm"])
            q_total = sum(o.half_totals["q_arm"])
            c_total = sum(o.half_totals["centromere"])
            cat = category_index(
                p_total=p_total,
                cen_total=c_total,
                q_total=q_total,
                has_start_tel=o.oriented_telo.start,
                has_stop_tel=o.oriented_telo.stop,
            )
            decorated.append(((cat, -o.contig.length, o.contig.contig_name), o))
        decorated.sort(key=lambda t: t[0])

        for _, o in decorated:
            new_name = f"{o.chromosome}_{o.contig.input_name}_{o.contig.contig_name}"
            if o.flipped:
                new_name += "_rc"
            stats = _stats_string(
                o.flipped_region_bins,
                o.region_start,
                o.region_end,
                o.oriented_telo.start,
                o.oriented_telo.stop,
            )
            rows.append(
                MapRow(
                    new_name=new_name,
                    original_name=o.contig.contig_name,
                    input_file=o.contig.input_file,
                    hap=o.contig.input_name,
                    chromosome=o.chromosome,
                    flipped=o.flipped,
                    length=o.contig.length,
                    stats=stats,
                )
            )
    return rows


# --- BED rewriting --------------------------------------------------


def _open_bed_in(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def _open_bed_out(path: Path, *, gzip_out: bool) -> IO[str]:
    if gzip_out:
        return gzip.open(path, "wt")
    return path.open("w")


def rewrite_bed(
    input_path: Path,
    output_path: Path,
    *,
    map_rows: list[MapRow],
    gzip_out: bool | None = None,
) -> None:
    """Rewrite a per-feature-set BED to use the scaffolded contig names.

    Reads ``input_path`` (the smoothed BED that ``annotate`` produced
    for one feature set), then for each entry in ``map_rows`` -- which
    is assumed already in the canonical chromosome x hap x category x
    length order produced by :func:`classify_and_orient` for one input
    -- emits the corresponding original contig's records under its
    new encoded name. Contigs flagged ``flipped=True`` have their
    intervals mirrored to ``[length - stop, length - start)`` and
    written in reverse order, so the output remains coordinate-sorted
    within each renamed contig.

    Contigs that the input BED knows about but the map doesn't list
    (because they were dropped by the length filter or had no leaf
    chromosome) are silently omitted from the output. Contigs in the
    map but absent from the input BED are also silently skipped (an
    input might legitimately have produced no records for some
    feature set -- e.g. all-novel sequences with no smoothing pass).
    """
    if gzip_out is None:
        gzip_out = str(output_path).endswith(".gz")

    by_contig: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    with _open_bed_in(input_path) as h:
        for i, raw in enumerate(h, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise ScaffoldError(
                    f"{input_path}:{i}: expected 3+ tab-separated columns, got {len(parts)}"
                )
            seq = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError as e:
                raise ScaffoldError(f"{input_path}:{i}: non-integer coordinates: {raw!r}") from e
            rest = "\t".join(parts[3:])
            by_contig[seq].append((start, end, rest))

    with _open_bed_out(output_path, gzip_out=gzip_out) as out:
        for row in map_rows:
            recs = by_contig.get(row.original_name)
            if recs is None:
                continue
            if row.flipped:
                # Reverse and mirror so the output stays coordinate-sorted
                # under the new (oriented) coordinate system.
                emit = [
                    (row.length - stop, row.length - start, rest)
                    for start, stop, rest in reversed(recs)
                ]
            else:
                emit = recs
            for start, end, rest in emit:
                if rest:
                    out.write(f"{row.new_name}\t{start}\t{end}\t{rest}\n")
                else:
                    out.write(f"{row.new_name}\t{start}\t{end}\n")
