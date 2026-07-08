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
import shutil
import tempfile
from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field, replace
from itertools import groupby
from pathlib import Path
from typing import IO

from karyoscope.core.io.agp import AgpComponent, AgpGap, AgpObject
from karyoscope.core.io.fasta import (
    read_fasta_records,
    reverse_complement,
    write_fasta_records,
)
from karyoscope.core.io.scaffold_map import MapRow, read_map
from karyoscope.core.io.telo import TeloFlags
from karyoscope.exceptions import ScaffoldError

#: Default number of N bases inserted between concatenated contigs when
#: ``--combine-chromosomes`` is on. Matches the scaffold subcommand
#: default; large enough that the per-contig smoothing's ``max_gap``
#: (1000) never bridges it, so each gap stays a clean ``novel`` run.
DEFAULT_SCAFFOLD_GAP_SIZE = 100_000

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

    Ported from the archive's ``calculate_need_to_flip`` (the newer
    version with the centroid-based both-tel rule for chrY-like
    cases), with one addition: acrocentric short-arm fragments are
    forced telomere-first (p-ter at the top) off the raw telomere
    flags, since their satellite/novel content breaks the region-block
    continuity the other branches rely on. The body is a boolean
    ladder over arm/centromere composition and telomere placement; see
    the inline comments for the intuition behind each branch.
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

    # Acrocentric short arms (chr13/14/15/21/22 p-arms: satellite, stalk,
    # rDNA, plus the p-ter telomere) must be drawn telomere-first, i.e.
    # with p-ter at the top. These short-arm fragments are not q-arm
    # bodies (``not is_q_like``); their satellite/novel content routinely
    # breaks region-block continuity, so -- unlike the branches below --
    # we key off the raw telomere flags rather than the ``continuous_*``
    # variants (gating on continuity is exactly what left these contigs
    # inconsistently oriented). We only force the single-telomere case:
    # the lone telomere is the p-ter and belongs at the top. Contigs with
    # telomeres on both ends fall through to the centroid logic so a
    # genuinely reversed full chromosome is still detected, and q-arm
    # bodies keep their q-ter telomere at the bottom via the ladder below.
    if is_acrocentric and (not is_q_like) and (has_start_tel != has_stop_tel):
        return has_stop_tel

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


# --- streaming per-contig access ------------------------------------
#
# Both BED rewriters need per-contig random access (the output is in
# scaffold-map order, which is a permutation of the input's contig
# order) while the input arrives grouped by contig in a *different*
# order and may be gzip-compressed (so we can't seek it). Loading the
# whole BED into a dict is O(total intervals) memory -- for a
# whole-genome region BED (~21M intervals) that is multiple GB and can
# OOM. Instead we stream the input once, spilling only the contigs the
# map actually places to a per-contig temp file, and read each back on
# demand. Peak memory is then O(one contig): a flipped contig is
# buffered to reverse it, an un-flipped one streams straight through.


class _SpilledContigs:
    """Per-contig temp files for the contigs a scaffold map places.

    Built by :func:`_spill_needed_contigs`. Acts as a context manager
    that removes the temp directory on exit.
    """

    __slots__ = ("_paths", "_tmpdir")

    def __init__(self, tmpdir: Path, paths: dict[str, Path]) -> None:
        self._tmpdir = tmpdir
        self._paths = paths

    def has(self, contig: str) -> bool:
        return contig in self._paths

    def read(self, contig: str) -> Iterator[tuple[int, int, str]]:
        """Yield ``(start, end, rest)`` for ``contig`` in input order.

        ``rest`` is the tab-joined 4th-and-later columns (``""`` when
        the input had exactly three columns).
        """
        path = self._paths[contig]
        with path.open("r") as h:
            for line in h:
                start_s, end_s, rest = line.rstrip("\n").split("\t", 2)
                yield int(start_s), int(end_s), rest

    def __enter__(self) -> _SpilledContigs:
        return self

    def __exit__(self, *exc: object) -> None:
        shutil.rmtree(self._tmpdir, ignore_errors=True)


def _spill_needed_contigs(input_path: Path, needed: set[str]) -> _SpilledContigs:
    """Stream ``input_path`` once, spilling each needed contig to a temp file.

    Contigs absent from ``needed`` are validated and skipped (never
    buffered). Malformed rows raise :class:`ScaffoldError` with the
    line number, matching the previous whole-file loader's behaviour.
    Each spilled line is ``start\\tend\\trest`` so the reader can round-
    trip the record without re-splitting the original name column.
    """
    tmpdir = Path(tempfile.mkdtemp(prefix="ks_scaffold_"))
    handles: dict[str, IO[str]] = {}
    paths: dict[str, Path] = {}
    try:
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
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError as e:
                    raise ScaffoldError(
                        f"{input_path}:{i}: non-integer coordinates: {raw!r}"
                    ) from e
                seq = parts[0]
                if seq not in needed:
                    continue
                fh = handles.get(seq)
                if fh is None:
                    # Index temp files by a counter so an odd contig name
                    # (slashes, spaces) can't produce an unsafe path.
                    path = tmpdir / f"{len(paths)}.bed"
                    paths[seq] = path
                    fh = path.open("w")
                    handles[seq] = fh
                rest = "\t".join(parts[3:])
                fh.write(f"{start}\t{end}\t{rest}\n")
    finally:
        for fh in handles.values():
            fh.close()
    return _SpilledContigs(tmpdir, paths)


def _write_bed_row(out: IO[str], name: str, start: int, end: int, rest: str) -> None:
    """Write one BED row, omitting the trailing tab when ``rest`` is empty."""
    if rest:
        out.write(f"{name}\t{start}\t{end}\t{rest}\n")
    else:
        out.write(f"{name}\t{start}\t{end}\n")


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

    Streaming: the input is read once, spilling only the placed contigs
    to per-contig temp files; each is then read back in map order. Peak
    memory is one contig, not the whole BED.
    """
    if gzip_out is None:
        gzip_out = str(output_path).endswith(".gz")

    needed = {row.original_name for row in map_rows}
    with (
        _spill_needed_contigs(input_path, needed) as spilled,
        _open_bed_out(output_path, gzip_out=gzip_out) as out,
    ):
        for row in map_rows:
            if not spilled.has(row.original_name):
                continue
            if row.flipped:
                # Reverse and mirror so the output stays coordinate-sorted
                # under the new (oriented) coordinate system. Only a
                # flipped contig needs buffering (to reverse it).
                recs = list(spilled.read(row.original_name))
                for start, stop, rest in reversed(recs):
                    _write_bed_row(
                        out, row.new_name, row.length - stop, row.length - start, rest
                    )
            else:
                for start, end, rest in spilled.read(row.original_name):
                    _write_bed_row(out, row.new_name, start, end, rest)


# --- applying an existing map to a foreign BED ----------------------
#
# ``rewrite_bed`` above is the in-pipeline remap: scaffold builds the map
# and immediately rewrites the BEDs it just annotated. The helper below
# is the *standalone* remap -- apply an already-built ``scaffold_map.tsv``
# to a BED that was annotated separately, possibly against a *different*
# database than the one used to derive the map (e.g. a cytoband-database
# annotation remapped with a map built from a region/roles database).
# This is what the ``karyoscope remap-bed`` command exposes.
#
# Because the map and the BED can come from independent runs, we validate
# they actually belong together before rewriting -- a mismatched pair
# would otherwise produce a silently corrupt BED (e.g. a stale contig
# length mirrors a flipped contig to the wrong coordinates).

#: Recognised FASTA suffixes (longest first), for deriving the input
#: "stem" the map's ``input_file`` was named from. Mirrors
#: ``scaffold_run._FASTA_EXTS`` (kept local to avoid importing the heavy
#: ``scaffold_run`` module into ``scaffold`` and creating an import cycle).
_FASTA_EXTS: tuple[str, ...] = (
    ".fasta.gz",
    ".fa.gz",
    ".fna.gz",
    ".fasta",
    ".fa",
    ".fna",
)


def _fasta_stem(name: str) -> str:
    """The output stem a FASTA basename produces (drops a known FASTA suffix)."""
    lower = name.lower()
    for ext in _FASTA_EXTS:
        if lower.endswith(ext):
            return name[: -len(ext)]
    return Path(name).stem


@dataclass(frozen=True)
class RemapStats:
    """Summary of a standalone BED remap (see :func:`remap_bed_with_map`)."""

    #: Distinct contig names in the input BED.
    bed_contigs: int
    #: Contigs listed in the scaffold map.
    map_contigs: int
    #: BED contigs that the map places (these are emitted, renamed).
    mapped_contigs: int
    #: BED contigs absent from the map (dropped from the output -- these are
    #: the short / unplaced contigs scaffolding legitimately excludes).
    dropped_contigs: int
    #: Map contigs with no records in the BED (skipped -- e.g. a feature set
    #: that produced nothing for that contig).
    map_contigs_absent_from_bed: int


def remap_bed_with_map(
    input_bed: Path,
    output_bed: Path,
    map_path: Path,
    *,
    gzip_out: bool | None = None,
    strict: bool = False,
) -> RemapStats:
    """Apply an existing ``scaffold_map.tsv`` to a separately-annotated BED.

    Reads the map, validates that the BED and map plausibly describe the
    same assembly, then delegates to :func:`rewrite_bed` (the single source
    of truth for the actual coordinate rewrite, including flipped-contig
    mirroring).

    Validation -- the BED carries no provenance pointer back to its source
    FASTA, so we rely on what the map records (``original_name`` + per-contig
    ``length``):

    * **Hard error** when *no* BED contig name appears in the map -- the two
      files describe different assemblies (or an argument was swapped).
    * **Hard error** when a contig present in *both* has a BED interval whose
      end exceeds the contig ``length`` recorded in the map -- the BED was
      annotated against a different (longer) sequence, which would corrupt
      the flip math.
    * **Warning** when the BED filename's stem doesn't match the map's source
      FASTA stem -- advisory only, since files are routinely renamed.

    A BED contig that is *absent* from the map is **not** an error: the map
    only lists contigs that survived scaffolding (length filter + leaf
    chromosome), while the original-coordinate BED contains every contig.
    Those are dropped from the output (as :func:`rewrite_bed` already does);
    their count is reported in :class:`RemapStats`.

    ``strict`` promotes the advisory conditions (stem mismatch; any map contig
    that has no records in the BED) to hard errors.
    """
    rows = read_map(map_path)
    if not rows:
        raise ScaffoldError(f"scaffold map is empty: {map_path}")

    map_length: dict[str, int] = {}
    for r in rows:
        # A well-formed map lists each original contig once; if a name repeats,
        # keep the largest length (the safest bound for the end<=length check).
        map_length[r.original_name] = max(map_length.get(r.original_name, 0), r.length)

    # Filename-stem advisory: the map's source FASTA stem should prefix the
    # BED basename (annotate names outputs '<stem>.<db>.<fs>.smoothed.bed[.gz]').
    map_stem = _fasta_stem(rows[0].input_file)
    if map_stem and not input_bed.name.startswith(f"{map_stem}."):
        msg = (
            f"BED {input_bed.name!r} does not look like it came from the same "
            f"assembly as the map (expected a name starting with {map_stem + '.'!r}, "
            f"the stem of the map's source FASTA {rows[0].input_file!r})"
        )
        if strict:
            raise ScaffoldError(msg + "; pass without --strict to treat this as a warning")
        logger.warning("%s", msg)

    # Single pass over the BED: distinct contig names + per-contig max end
    # (only the contigs the map knows about can violate the length bound).
    bed_names: set[str] = set()
    max_end: dict[str, int] = {}
    with _open_bed_in(input_bed) as h:
        for i, raw in enumerate(h, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise ScaffoldError(
                    f"{input_bed}:{i}: expected 3+ tab-separated columns, got {len(parts)}"
                )
            seq = parts[0]
            bed_names.add(seq)
            if seq in map_length:
                try:
                    end = int(parts[2])
                except ValueError as e:
                    raise ScaffoldError(f"{input_bed}:{i}: non-integer coordinates: {raw!r}") from e
                if end > max_end.get(seq, 0):
                    max_end[seq] = end

    overlap = bed_names & set(map_length)
    if not overlap:
        raise ScaffoldError(
            f"no contig in {input_bed.name} matches the scaffold map {map_path.name} "
            f"({len(bed_names)} BED contigs, {len(map_length)} map contigs); "
            "these do not appear to describe the same assembly"
        )

    violations = [
        (c, max_end[c], map_length[c]) for c in overlap if max_end.get(c, 0) > map_length[c]
    ]
    if violations:
        sample = ", ".join(f"{c} (end {e} > length {ln})" for c, e, ln in sorted(violations)[:5])
        raise ScaffoldError(
            f"{input_bed.name} has intervals beyond the contig lengths recorded in "
            f"{map_path.name} for {len(violations)} contig(s): {sample}"
            + ("" if len(violations) <= 5 else ", ...")
            + ". The BED was annotated against a different assembly than the map."
        )

    map_absent = [r.original_name for r in rows if r.original_name not in bed_names]
    if map_absent and strict:
        raise ScaffoldError(
            f"{len(map_absent)} map contig(s) have no records in {input_bed.name} "
            f"(e.g. {', '.join(sorted(map_absent)[:5])}); pass without --strict to allow"
        )

    stats = RemapStats(
        bed_contigs=len(bed_names),
        map_contigs=len(map_length),
        mapped_contigs=len(overlap),
        dropped_contigs=len(bed_names - set(map_length)),
        map_contigs_absent_from_bed=len(map_absent),
    )
    logger.info(
        "remap %s: %d/%d BED contigs placed by map; %d dropped (unscaffolded); "
        "%d map contigs had no BED records",
        input_bed.name,
        stats.mapped_contigs,
        stats.bed_contigs,
        stats.dropped_contigs,
        stats.map_contigs_absent_from_bed,
    )
    rewrite_bed(input_bed, output_bed, map_rows=rows, gzip_out=gzip_out)
    return stats


# --- FASTA rewriting ------------------------------------------------


def rewrite_fasta(
    input_path: Path,
    output_path: Path,
    *,
    map_rows: list[MapRow],
    keep_unscaffolded: bool = True,
    gzip_out: bool | None = None,
    line_width: int | None = None,
) -> None:
    """Write a scaffolded FASTA from ``input_path``.

    Walks ``map_rows`` in order, emitting each scaffolded contig
    under its encoded name (reverse-complementing when ``flipped`` is
    True). When ``keep_unscaffolded`` is True (the default), any
    contig present in the source FASTA but absent from ``map_rows``
    is appended at the end under its original name. This matches the
    archive's ``scaffold_hap_assembly.py`` behaviour and keeps the
    output assembly complete; pass ``False`` to drop unscaffolded
    contigs entirely.

    ``gzip_out=None`` gzips iff ``output_path`` ends in ``.gz``.
    ``line_width=None`` writes each sequence on a single line.

    Contigs in ``map_rows`` whose ``original_name`` is missing from
    the source FASTA are silently skipped -- the same forgiving
    semantics as :func:`rewrite_bed`.
    """
    records = read_fasta_records(input_path)

    placed: set[str] = set()
    out_records: dict[str, str] = {}

    for row in map_rows:
        seq = records.get(row.original_name)
        if seq is None:
            continue
        if row.flipped:
            seq = reverse_complement(seq)
        out_records[row.new_name] = seq
        placed.add(row.original_name)

    if keep_unscaffolded:
        # Preserve source order for the appended contigs.
        for name, seq in records.items():
            if name in placed:
                continue
            out_records[name] = seq

    write_fasta_records(out_records, output_path, gzip_out=gzip_out, line_width=line_width)


# --- combined-chromosome scaffolding --------------------------------
#
# When ``--combine-chromosomes`` is on, every contig of one
# ``(chromosome, haplotype)`` group is concatenated into a single output
# sequence named ``<chrom>_<hap>``, with ``gap_size`` N bases between
# adjacent contigs. The FASTA, the per-feature-set BEDs, and the AGP all
# derive from one shared layout (:func:`plan_combined_layout`) so their
# coordinate systems agree exactly.
#
# Coordinate model (see the scaffold subcommand's design notes): each
# per-contig BED tiles ``[0, E)`` where ``E = L - k + 1`` (``L`` = true
# sequence length, ``k`` = k-mer size). The map's ``length`` field is
# ``E``. We never hardcode ``k``: object offsets come from true
# sequence lengths (read from the FASTA), and each contig's tiling end
# comes from its own ``E``. The ``novel`` gap interval between two
# contigs therefore spans ``[offset_i + E_i, offset_{i+1})`` -- the
# literal N gap plus the ``k-1`` untiled tail of contig ``i`` whose
# k-mers, in the concatenated assembly, would overlap the Ns. This is
# byte-identical to re-annotating the combined FASTA, and stays correct
# for a future variable-k database with no special case.


@dataclass(frozen=True)
class PlacedComponent:
    """One contig placed inside a combined object.

    ``object_start`` is the 0-based offset (in true bp) of the contig
    within its output object. ``true_length`` is the contig's full
    sequence length; ``bed_extent`` is ``E`` (where the per-contig BED
    tiling stops, i.e. the map's ``length``).
    """

    row: MapRow
    object_start: int
    true_length: int
    bed_extent: int


@dataclass(frozen=True)
class ScaffoldObject:
    """One output FASTA record: a renamed singleton or a combined chromosome.

    ``combined`` is True when the object concatenates a whole
    ``(chromosome, hap)`` group under a simplified ``<chrom>_<hap>``
    name; False for a singleton emitted under a ``<chrom>_<hap>_<A|B|C...>``
    name (an acrocentric group left uncombined: its contigs stay as
    separate records but are renamed in canonical order). ``gap_size`` is
    the N-run length between components (irrelevant when there is only one).
    """

    name: str
    components: list[PlacedComponent]
    gap_size: int
    combined: bool


def _column_label(index: int) -> str:
    """Spreadsheet-style bijective base-26 label: 0->A, 25->Z, 26->AA, ...

    Used to suffix the contigs of an uncombined acrocentric group; a
    group never has enough contigs to reach ``AA`` in practice, but the
    bijective scheme keeps the labels unambiguous if it ever does.
    """
    label = ""
    index += 1
    while index:
        index, rem = divmod(index - 1, 26)
        label = chr(ord("A") + rem) + label
    return label


def _acrocentric_singleton_name(chrom: str, hap: str, group_index: int) -> str:
    """Name for one contig of an uncombined acrocentric ``(chrom, hap)`` group.

    Acrocentric groups left uncombined keep their contigs as separate
    records, but -- like the combined groups -- drop the original contig
    name and any ``_rc`` orientation suffix in favour of a clean
    ``<chrom>_<hap>_<A|B|C...>`` label. ``group_index`` is the contig's
    position within its ``(chrom, hap)`` group in canonical map order, so
    the same letter is produced by both :func:`plan_combined_layout` (the
    FASTA/BED/AGP object name) and :func:`combined_map_rows` (the renderer's
    join key).
    """
    return f"{chrom}_{hap}_{_column_label(group_index)}"


def plan_combined_layout(
    map_rows: list[MapRow],
    true_lengths: dict[str, int],
    *,
    gap_size: int = DEFAULT_SCAFFOLD_GAP_SIZE,
    acrocentrics: set[str] = DEFAULT_HUMAN_ACROCENTRICS,
    combine_acrocentrics: bool = False,
) -> list[ScaffoldObject]:
    """Plan the combined-chromosome output layout for one input.

    Groups ``map_rows`` by ``(chromosome, hap)`` in map order (which is
    already the canonical chromosome x hap x category x length order),
    and turns each group into one or more :class:`ScaffoldObject`:

    * A non-acrocentric group -- or an acrocentric group when
      ``combine_acrocentrics`` is True -- becomes a single combined
      object named ``<chrom>_<hap>`` whose components are placed at
      cumulative offsets ``Σ(true_length + gap_size)``.
    * An acrocentric group left uncombined becomes one singleton object
      per contig, each named ``<chrom>_<hap>_<A|B|C...>`` by canonical
      order (dropping the original contig name and ``_rc`` suffix, just
      like the combined groups do).

    Rows whose ``original_name`` is absent from ``true_lengths`` (the
    contig was not in the source FASTA) are skipped, mirroring the
    forgiving semantics of :func:`rewrite_fasta` / :func:`rewrite_bed`.
    """
    # Group preserving first-seen order; map_rows already arrive grouped.
    groups: dict[tuple[str, str], list[MapRow]] = {}
    for row in map_rows:
        groups.setdefault((row.chromosome, row.hap), []).append(row)

    objects: list[ScaffoldObject] = []
    for (chrom, hap), rows in groups.items():
        present = [r for r in rows if r.original_name in true_lengths]
        if not present:
            continue
        is_acro = chrom in acrocentrics
        combine = not (is_acro and not combine_acrocentrics)

        if combine:
            components: list[PlacedComponent] = []
            offset = 0
            for r in present:
                length = true_lengths[r.original_name]
                components.append(
                    PlacedComponent(
                        row=r,
                        object_start=offset,
                        true_length=length,
                        bed_extent=r.length,
                    )
                )
                offset += length + gap_size
            objects.append(
                ScaffoldObject(
                    name=f"{chrom}_{hap}",
                    components=components,
                    gap_size=gap_size,
                    combined=True,
                )
            )
        else:
            # Index over the full group (not ``present``) so a contig's
            # letter is stable even if an earlier sibling was dropped for
            # being absent from ``true_lengths``; that matches the index
            # ``combined_map_rows`` uses, which has no length filter.
            for idx, r in enumerate(rows):
                if r.original_name not in true_lengths:
                    continue
                length = true_lengths[r.original_name]
                objects.append(
                    ScaffoldObject(
                        name=_acrocentric_singleton_name(chrom, hap, idx),
                        components=[
                            PlacedComponent(
                                row=r,
                                object_start=0,
                                true_length=length,
                                bed_extent=r.length,
                            )
                        ],
                        gap_size=gap_size,
                        combined=False,
                    )
                )
    return objects


def combined_map_rows(
    map_rows: list[MapRow],
    *,
    acrocentrics: set[str] = DEFAULT_HUMAN_ACROCENTRICS,
    combine_acrocentrics: bool = False,
) -> list[MapRow]:
    """Synthetic map rows describing the combined-chromosome objects.

    Mirrors the object naming and the ``(chromosome, hap)`` grouping of
    :func:`plan_combined_layout`, but needs no FASTA: it produces one
    :class:`MapRow` per output object so the karyotype renderer can join
    the combined binned BED (keyed by object name) to chromosome / hap /
    telomere metadata.

    A combined object gets a synthetic ``new_name`` of ``<chrom>_<hap>``
    (matching the combined BED's sequence name) and a ``stats`` proxy
    carrying only the end-telomere flags the renderer reads: a leading
    ``T`` iff the first component starts with a telomere, a trailing
    ``T`` iff the last component ends with one. Acrocentric groups left
    uncombined emit one row per contig, each renamed to the same
    ``<chrom>_<hap>_<A|B|C...>`` label the layout gives its singleton
    object, so the renderer's join key matches.

    The ``length`` field is the sum of component ``bed_extent`` values
    (informational only -- the renderer derives sequence lengths from the
    binned BED, not from this field). ``flipped`` is always ``False`` on
    a combined object: orientation was already baked into the component
    coordinates when the combined BED was written.
    """
    groups: dict[tuple[str, str], list[MapRow]] = {}
    for row in map_rows:
        groups.setdefault((row.chromosome, row.hap), []).append(row)

    out: list[MapRow] = []
    for (chrom, hap), rows in groups.items():
        is_acro = chrom in acrocentrics
        combine = not (is_acro and not combine_acrocentrics)
        if not combine:
            # One synthetic row per contig, renamed to the clean
            # <chrom>_<hap>_<A|B|C...> label the layout emits for the
            # singleton object; all other fields (orientation, telomere
            # stats) are preserved so the renderer sees them.
            for idx, r in enumerate(rows):
                out.append(replace(r, new_name=_acrocentric_singleton_name(chrom, hap, idx)))
            continue
        start_t = rows[0].stats.startswith("T")
        stop_t = rows[-1].stats.endswith("T")
        stats = ("T" if start_t else "x") + ("T" if stop_t else "x")
        out.append(
            MapRow(
                new_name=f"{chrom}_{hap}",
                original_name=f"{chrom}_{hap}",
                input_file=rows[0].input_file,
                hap=hap,
                chromosome=chrom,
                flipped=False,
                length=sum(r.length for r in rows),
                stats=stats,
            )
        )
    return out


def _to_agp_objects(
    objects: list[ScaffoldObject],
    leftovers: list[tuple[str, int]],
) -> list[AgpObject]:
    """Build AGP objects from the layout plus any unscaffolded leftovers.

    ``leftovers`` is ``[(name, length)]`` for contigs kept in the output
    FASTA under their original names (``keep_unscaffolded``). Each is a
    one-component object so the AGP fully describes the output FASTA.
    """
    agp: list[AgpObject] = []
    for obj in objects:
        parts: list[AgpComponent | AgpGap] = []
        for i, comp in enumerate(obj.components):
            parts.append(
                AgpComponent(
                    component_id=comp.row.original_name,
                    object_start=comp.object_start,
                    object_end=comp.object_start + comp.true_length,
                    length=comp.true_length,
                    orientation="-" if comp.row.flipped else "+",
                )
            )
            if i + 1 < len(obj.components):
                gap_start = comp.object_start + comp.true_length
                gap_end = obj.components[i + 1].object_start
                parts.append(
                    AgpGap(
                        object_start=gap_start,
                        object_end=gap_end,
                        length=gap_end - gap_start,
                    )
                )
        agp.append(AgpObject(name=obj.name, parts=parts))

    for name, length in leftovers:
        agp.append(
            AgpObject(
                name=name,
                parts=[
                    AgpComponent(
                        component_id=name,
                        object_start=0,
                        object_end=length,
                        length=length,
                        orientation="+",
                    )
                ],
            )
        )
    return agp


def write_combined_fasta(
    records: dict[str, str],
    objects: list[ScaffoldObject],
    output_path: Path,
    *,
    keep_unscaffolded: bool = True,
    gzip_out: bool | None = None,
    line_width: int | None = None,
) -> list[tuple[str, int]]:
    """Write the combined-chromosome FASTA and return the leftover list.

    Each :class:`ScaffoldObject` is emitted under its name as the
    concatenation of its (oriented) component sequences joined by
    ``gap_size`` Ns. Singleton objects join one sequence, so no gap is
    inserted. Contigs kept by ``keep_unscaffolded`` are appended under
    their original names.

    Returns ``[(name, length)]`` for the appended leftovers, so the
    caller can hand them to :func:`_to_agp_objects` and keep the AGP a
    complete description of the FASTA.
    """
    placed: set[str] = set()
    out_records: dict[str, str] = {}

    for obj in objects:
        gap = "N" * obj.gap_size
        seqs: list[str] = []
        for comp in obj.components:
            seq = records.get(comp.row.original_name)
            if seq is None:
                continue
            if comp.row.flipped:
                seq = reverse_complement(seq)
            seqs.append(seq)
            placed.add(comp.row.original_name)
        if not seqs:
            continue
        out_records[obj.name] = gap.join(seqs)

    leftovers: list[tuple[str, int]] = []
    if keep_unscaffolded:
        for name, seq in records.items():
            if name in placed:
                continue
            out_records[name] = seq
            leftovers.append((name, len(seq)))

    write_fasta_records(out_records, output_path, gzip_out=gzip_out, line_width=line_width)
    return leftovers


class _ObjectCoalescer:
    """Coalesce and write one combined object's BED rows as a stream.

    Concatenating per-contig BEDs only creates new adjacencies at the
    junctions (e.g. a contig ending in ``novel``, the inserted ``novel``
    gap, and the next contig starting in ``novel``). Merging them keeps
    the combined BED identical to what annotate would emit for the
    concatenated sequence.

    Holds a single pending ``(start, end, rest)`` and extends it when
    the next interval abuts and shares its label, else flushes it under
    the object's name -- so peak memory is one row, not the whole object.
    """

    __slots__ = ("_end", "_name", "_out", "_rest", "_start")

    def __init__(self, out: IO[str], name: str) -> None:
        self._out = out
        self._name = name
        self._start = -1
        self._end = -1
        self._rest: str | None = None

    def push(self, start: int, end: int, rest: str) -> None:
        if self._rest is not None and self._end == start and self._rest == rest:
            self._end = end
            return
        if self._rest is not None:
            _write_bed_row(self._out, self._name, self._start, self._end, self._rest)
        self._start = start
        self._end = end
        self._rest = rest

    def flush(self) -> None:
        if self._rest is not None:
            _write_bed_row(self._out, self._name, self._start, self._end, self._rest)
        self._rest = None


def rewrite_bed_combined(
    input_path: Path,
    output_path: Path,
    *,
    objects: list[ScaffoldObject],
    gzip_out: bool | None = None,
) -> None:
    """Rewrite a per-feature-set BED into combined-chromosome coordinates.

    For each :class:`ScaffoldObject`, every component's intervals are
    shifted by the component's ``object_start`` (flipped contigs are
    mirrored within ``[0, E)`` first, as in :func:`rewrite_bed`), and a
    ``novel`` interval is inserted between consecutive components to
    fill the N gap plus the untiled ``k-1`` boundary. The result tiles
    ``[0, object_extent)`` per object with no gaps, then adjacent
    same-label intervals are coalesced.

    A component whose contig is absent from the input BED has its whole
    ``[object_start, object_start + E)`` extent filled with ``novel`` so
    the tiling stays complete (this does not happen for real annotate
    output, which tiles every contig fully).

    Streaming: the input is read once, spilling only the placed contigs
    to per-contig temp files. Each object's tiling is coalesced on the
    fly (a running ``_merge_adjacent`` fold) and written directly, so
    peak memory is one contig rather than the whole BED.
    """
    if gzip_out is None:
        gzip_out = str(output_path).endswith(".gz")

    needed = {comp.row.original_name for obj in objects for comp in obj.components}
    with (
        _spill_needed_contigs(input_path, needed) as spilled,
        _open_bed_out(output_path, gzip_out=gzip_out) as out,
    ):
        for obj in objects:
            # Coalesce the object's tiling on the fly (a running
            # _merge_adjacent fold) and write directly, so we never hold
            # the whole object's intervals -- only one pending row.
            coalescer = _ObjectCoalescer(out, obj.name)
            n = len(obj.components)
            for idx, comp in enumerate(obj.components):
                off = comp.object_start
                if not spilled.has(comp.row.original_name):
                    # No records for this contig in this feature set: fill
                    # its extent with novel to keep the tiling complete.
                    logger.warning(
                        "contig %r has no records in %s; filling its %d bp "
                        "with novel in the combined BED",
                        comp.row.original_name,
                        input_path.name,
                        comp.bed_extent,
                    )
                    coalescer.push(off, off + comp.bed_extent, "novel")
                elif comp.row.flipped:
                    # Mirror within [0, E) and reverse; only a flipped
                    # contig needs buffering.
                    recs = list(spilled.read(comp.row.original_name))
                    ext = comp.bed_extent
                    for start, stop, rest in reversed(recs):
                        coalescer.push(off + ext - stop, off + ext - start, rest)
                else:
                    for start, end, rest in spilled.read(comp.row.original_name):
                        coalescer.push(off + start, off + end, rest)
                if idx + 1 < n:
                    gap_start = off + comp.bed_extent
                    gap_end = obj.components[idx + 1].object_start
                    if gap_end > gap_start:
                        coalescer.push(gap_start, gap_end, "novel")
            coalescer.flush()
