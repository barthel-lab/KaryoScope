"""Render an assembly karyotype as SVG.

Faithful port of the archive's ``KaryoScope_assembly.py`` adapted for
the v0.1 data model: instead of a single combined per-sample BED +
two-column ``scaffold_stats.tsv``, the renderer consumes a list of
per-input :class:`RenderInput` records, each carrying the per-input
scaffold map (the source-of-truth file Stage 5d-1 introduced), the
per-input binned scaffolded BED, and optionally per-input centromere
coordinates for the centromere mode.

Three render modes:

* ``"genome"`` -- whole-chromosome view. Pixels-per-bp tuned so a
  human chromosome fits comfortably; 10 Mb scale bar.
* ``"subtelomere"`` -- zoomed view of the p-arm and q-arm telomeric
  ends. Only contigs flagged with a telomere are drawn. Scale 1 px
  per 300 bp; 10 kb scale bar.
* ``"centromere"`` -- zoomed view of each contig's centromere. Only
  contigs that appear in the centromere coordinates table are drawn.
  Scale 1 px per 25 kb; 1 Mb scale bar.

Sex-determination filtering follows the archive: per chromosome,
``get_expected_haps`` returns the haps that should be drawn given
the sample's sex and the chosen sex-determination system (XY / X0 /
ZW / ZO). Combined with the actual haps present in the data, this
decides which (chrom, hap) cells get an empty column vs a populated
column vs are skipped entirely.

The renderer writes a single SVG; layout is computed top-down across
chromosomes left-to-right, with haplotypes within each chromosome
laid out side by side.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from math import floor
from pathlib import Path
from typing import Literal

import drawsvg as draw

from karyoscope.core.io.scaffold_map import MapRow
from karyoscope.core.scaffold import Interval, chromosome_sort_key
from karyoscope.exceptions import KaryotypeError

logger = logging.getLogger(__name__)


#: The three render modes.
Mode = Literal["genome", "subtelomere", "centromere"]
_VALID_MODES: tuple[Mode, ...] = ("genome", "subtelomere", "centromere")


#: Recognised sex-determination systems. Each declares the
#: heterogametic sex (the one with two different sex chromosomes), the
#: homogametic sex chromosome, and the full list of sex chromosomes
#: (in the order ``[heterogametic_first, homogametic]``).
PREDEFINED_SEX_SYSTEMS: dict[str, dict[str, object]] = {
    "XY": {
        "heterogametic_sex": "male",
        "homogametic_chromosome": "chrX",
        "sex_chromosomes": ["chrY", "chrX"],
    },
    "X0": {
        "heterogametic_sex": "male",
        "homogametic_chromosome": "chrX",
        "sex_chromosomes": ["chrX"],
    },
    "ZW": {
        "heterogametic_sex": "female",
        "homogametic_chromosome": "chrZ",
        "sex_chromosomes": ["chrZ", "chrW"],
    },
    "ZO": {
        "heterogametic_sex": "female",
        "homogametic_chromosome": "chrZ",
        "sex_chromosomes": ["chrZ"],
    },
}


#: Default human chromosome list. We seed the renderer's CHROMOSOMES
#: with these so that human assemblies show an empty column for any
#: chromosome that's missing in the actual data (rather than just
#: silently omitting it). For non-human assemblies the chromosomes
#: are discovered from the map files.
DEFAULT_HUMAN_CHROMOSOMES: tuple[str, ...] = tuple(
    [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
)


# --- per-mode visual parameters -----------------------------------------


@dataclass(frozen=True)
class _ModeParams:
    pixels_per_pos: float
    scale_bar_length: int
    scale_bar_label: str


_MODE_PARAMS: dict[str, _ModeParams] = {
    "genome": _ModeParams(
        pixels_per_pos=4 / 1_000_000, scale_bar_length=10_000_000, scale_bar_label="10 Mbp"
    ),
    "subtelomere": _ModeParams(
        pixels_per_pos=1 / 300, scale_bar_length=10_000, scale_bar_label="10 Kbp"
    ),
    "centromere": _ModeParams(
        pixels_per_pos=1 / 25_000, scale_bar_length=1_000_000, scale_bar_label="1 Mbp"
    ),
}


# --- inputs --------------------------------------------------------------


@dataclass
class RenderInput:
    """One input file's worth of renderer-ready data.

    ``map_rows`` provides the chromosome / hap / stats triples per
    contig (the map file is the source of truth for those, rather
    than parsing the encoded ``<chrom>_<hap>_<contig>`` name).
    ``binned_bed`` is the per-contig list of binned intervals at the
    mode-appropriate bin size. ``centromere_ranges`` is required
    only for centromere mode.
    """

    map_rows: list[MapRow]
    binned_bed: dict[str, list[Interval]]
    centromere_ranges: dict[str, tuple[int, int]] | None = None


# --- sex-aware hap filtering --------------------------------------------


def get_expected_haps(
    chromosome: str,
    sex: str | None,
    haplotypes: list[str],
    sex_determination_system: str | dict,
) -> list[str]:
    """Return the haplotypes that should appear for ``chromosome``.

    Direct port of the archive's logic. When ``sex`` is ``None`` the
    sex chromosomes get no expectations (so they only appear if there
    is actual data); autosomes get the full haplotype list.
    """
    if isinstance(sex_determination_system, str):
        system = PREDEFINED_SEX_SYSTEMS.get(sex_determination_system.upper())
        if not system:
            raise KaryotypeError(
                f"unknown sex-determination system {sex_determination_system!r}; "
                f"expected one of {sorted(PREDEFINED_SEX_SYSTEMS)}"
            )
    elif isinstance(sex_determination_system, dict):
        system = sex_determination_system
    else:
        raise KaryotypeError("sex_determination_system must be a string or a dict")

    sex_chromosomes = list(system["sex_chromosomes"])  # type: ignore[arg-type]
    heterogametic_sex: str = system["heterogametic_sex"]  # type: ignore[assignment]
    homogametic_chromosome: str = system["homogametic_chromosome"]  # type: ignore[assignment]

    if sex is None:
        # No expectations on sex chromosomes -- only actual data shows up.
        if chromosome in sex_chromosomes:
            return []
        return haplotypes

    sex_lower = sex.lower()
    if sex_lower == "reference":
        return [haplotypes[0]] if haplotypes else []
    if chromosome not in sex_chromosomes:
        return haplotypes
    is_heterogametic = sex_lower == heterogametic_sex.lower()
    if is_heterogametic:
        try:
            idx = sex_chromosomes.index(chromosome)
        except ValueError:
            return []
        return [haplotypes[idx]] if idx < len(haplotypes) else []
    return haplotypes if chromosome == homogametic_chromosome else []


# --- helpers -------------------------------------------------------------


def _telomere_flags_from_stats(stats: str) -> tuple[bool, bool]:
    """Read TPCQT-style stats string; first char 'T' means start telomere,
    last char 'T' means stop telomere."""
    has_start = stats.startswith("T")
    has_stop = stats.endswith("T") and len(stats) > 0
    return has_start, has_stop


def _haps_natural_sort_key(hap: str) -> tuple[int, int, str]:
    """Order haplotypes naturally: hap1, hap2, ..., then alphabetical, then 'unassigned' last."""
    if hap == "unassigned":
        return (2, 0, "")
    if hap.startswith("hap") and hap[3:].isdigit():
        return (0, int(hap[3:]), "")
    return (1, 0, hap)


# --- main renderer -------------------------------------------------------


@dataclass
class _CellLayout:
    chrom_block_widths: dict[str, int] = field(default_factory=dict)
    chrom_start_x: dict[str, int] = field(default_factory=dict)
    hap_block_widths: dict[str, dict[str, int]] = field(default_factory=dict)
    hap_start_x: dict[str, dict[str, int]] = field(default_factory=dict)
    drawable_haps_per_chrom: dict[str, list[str]] = field(default_factory=dict)
    image_width: int = 0


def render_karyotype(
    inputs: list[RenderInput],
    *,
    colors: dict[str, str],
    mode: Mode = "genome",
    sex: str | None = None,
    sex_determination_system: str | dict = "XY",
    background_color: str = "white",
    subtelomere_boundary: int = 250_000,
    max_num_sequences: int = 400,
    seed_human_chromosomes: bool = True,
    output_path: Path,
    sample_label: str | None = None,
    database_id: str | None = None,
    feature_set_label: str | None = None,
    smoothed: bool = True,
    show_title: bool = True,
    show_legend: bool = True,
) -> None:
    """Render and save a karyotype SVG.

    Parameters
    ----------
    inputs
        Per-input :class:`RenderInput` records. Each contig in
        ``binned_bed`` must have a corresponding row in ``map_rows``
        (matched on ``new_name`` == BED first column).
    colors
        ``{feature: hex_color}`` for the requested feature set. The
        ``novel`` sentinel and the no-colour-fallback (white) are
        handled internally.
    mode
        Which view to render. See module docstring.
    sex
        Sample sex (``"male"``, ``"female"``, ``"reference"``) or
        ``None`` for unknown.
    sex_determination_system
        ``"XY"`` (default), ``"X0"``, ``"ZW"``, ``"ZO"``, or a
        custom dict matching :data:`PREDEFINED_SEX_SYSTEMS`'s shape.
    background_color
        ``"white"`` (draws sequence outlines) or ``"black"`` (omits
        outlines, uses white text and scale bar for contrast).
    subtelomere_boundary
        Subtelomere window size in bp (subtelomere mode only).
    max_num_sequences
        Cap on number of contigs drawn; surplus contigs are dropped
        with a warning. Matches the archive default.
    seed_human_chromosomes
        When True (default), seed the CHROMOSOMES list with the
        standard human set so missing chromosomes still show empty
        columns. Disable for non-human assemblies.
    output_path
        Where to write the SVG.
    sample_label, database_id, feature_set_label, smoothed
        Metadata fields rendered in the title band at the top of the
        SVG when ``show_title=True``. Each is optional; missing fields
        are omitted from the title rather than left as placeholders.
    show_title
        Draw the title band at the top. Default True.
    show_legend
        Draw the color legend in the right margin. Only features that
        actually appear in the rendered data are listed (so the
        legend stays compact for large feature sets). Default True.
    """
    if mode not in _VALID_MODES:
        raise KaryotypeError(f"unknown mode {mode!r}; expected one of {_VALID_MODES}")

    mode_params = _MODE_PARAMS[mode]
    pixels_per_pos = mode_params.pixels_per_pos

    # --- aggregate all inputs into rendered-set views ----------------

    # New_name -> MapRow (across all inputs).
    map_by_name: dict[str, MapRow] = {}
    for ri in inputs:
        for row in ri.map_rows:
            map_by_name[row.new_name] = row

    # Telomere flags from the stats string.
    tel_start_sequences: set[str] = set()
    tel_stop_sequences: set[str] = set()
    for name, row in map_by_name.items():
        has_start, has_stop = _telomere_flags_from_stats(row.stats)
        if has_start:
            tel_start_sequences.add(name)
        if has_stop:
            tel_stop_sequences.add(name)

    # Per-sequence (contig) length from the binned BED -- max stop seen.
    sequence_lengths: dict[str, int] = {}
    for ri in inputs:
        for name, intervals in ri.binned_bed.items():
            if not intervals:
                continue
            seq_max = max(stop for _, stop, _ in intervals)
            sequence_lengths[name] = max(sequence_lengths.get(name, 0), seq_max)

    # CHROMOSOMES and HAPLOTYPES.
    chroms_seen: set[str] = set()
    haps_seen: set[str] = set()
    for row in map_by_name.values():
        chroms_seen.add(row.chromosome)
        haps_seen.add(row.hap)

    CHROMOSOMES: list[str] = []
    if seed_human_chromosomes:
        CHROMOSOMES = list(DEFAULT_HUMAN_CHROMOSOMES)
    for c in chroms_seen:
        if c not in CHROMOSOMES:
            CHROMOSOMES.append(c)
    CHROMOSOMES.sort(key=chromosome_sort_key)

    HAPLOTYPES: list[str] = sorted(haps_seen, key=_haps_natural_sort_key)
    if not HAPLOTYPES:
        HAPLOTYPES = ["hap1", "hap2"]

    # Centromere coordinates (centromere mode only).
    centromere_coords: dict[str, tuple[int, int]] = {}
    if mode == "centromere":
        for ri in inputs:
            if ri.centromere_ranges is None:
                continue
            for name, (cstart, cend) in ri.centromere_ranges.items():
                if name in centromere_coords:
                    existing_start, existing_end = centromere_coords[name]
                    centromere_coords[name] = (
                        min(existing_start, cstart),
                        max(existing_end, cend),
                    )
                else:
                    centromere_coords[name] = (cstart, cend)
        if not centromere_coords:
            raise KaryotypeError(
                "centromere mode requires centromere coordinates; none were provided "
                "across any of the inputs"
            )
        max_centromere_length = max(end - start for start, end in centromere_coords.values())
    else:
        max_centromere_length = 0

    # --- layout constants (from archive) -----------------------------

    min_label_width = 15 if len(HAPLOTYPES) > 1 else 25
    initial_x = 50
    sequence_gap = 15
    hap_gap = 10
    chrom_gap = 24
    circle_radius = 3
    x_border = initial_x + sequence_gap - 1
    y_border = 25
    # Title band (optional, drawn above the karyotype). All
    # below-title-band y constants get offset by ``title_band_height``
    # via :data:`title_offset`.
    title_band_height = 35 if show_title else 0
    title_offset = title_band_height
    chrom_label_y = 23 + title_offset
    chrom_line_y = 28 + title_offset
    hap_label_y = 43 + title_offset
    hap_line_y = 48 + title_offset
    initial_y = (60 if len(HAPLOTYPES) > 1 else 40) + title_offset
    min_q_arm_offset = 5
    text_color = "#000000" if background_color == "white" else "#FFFFFF"

    # Legend band (optional, drawn at the right margin). The width is
    # added to ``image_width`` below. Per-row height matches the small
    # text size we use for feature labels.
    legend_band_width = 170 if show_legend else 0
    legend_row_height = 16
    legend_swatch_size = 12
    legend_text_size = 11
    legend_pad_x = 18  # gap between karyotype right edge and legend

    P_Q_ARM_GAP = 50
    if mode == "subtelomere":
        p_arm_height = subtelomere_boundary * pixels_per_pos
        q_arm_start_y = initial_y + p_arm_height + P_Q_ARM_GAP
        q_arm_height = subtelomere_boundary * pixels_per_pos
        final_image_height = int(q_arm_start_y + q_arm_height + initial_y)
    elif mode == "centromere":
        final_image_height = int(
            initial_y + (max_centromere_length * pixels_per_pos) + y_border + initial_y
        )
    else:  # genome
        final_image_height = 0  # computed below from max_stop_y

    # --- helpers depending on layout consts --------------------------

    def pos_to_y(pos: int) -> int:
        return int(initial_y + floor(pos * pixels_per_pos))

    # --- Pass 1: collect sequences to plot ---------------------------

    telomeric_sequences_set = tel_start_sequences | tel_stop_sequences
    sequences_to_plot: list[str] = []
    seen_sequences: set[str] = set()
    max_stop_y = 0

    # We iterate per-input in input order; within each input, map order
    # (canonical chromosome x hap x category x length).
    for ri in inputs:
        for row in ri.map_rows:
            seq = row.new_name
            if mode == "subtelomere" and seq not in telomeric_sequences_set:
                continue
            if mode == "centromere" and seq not in centromere_coords:
                continue
            if seq in seen_sequences:
                continue
            if len(sequences_to_plot) >= max_num_sequences:
                logger.warning(
                    "max_num_sequences=%d reached; dropping contig %r",
                    max_num_sequences,
                    seq,
                )
                continue
            sequences_to_plot.append(seq)
            seen_sequences.add(seq)
            if mode == "genome":
                stop = sequence_lengths.get(seq, 0)
                if stop:
                    max_stop_y = max(max_stop_y, pos_to_y(stop))

    if mode == "genome":
        final_image_height = int(max_stop_y + y_border)

    # --- group by chromosome and haplotype ---------------------------

    sequences_per_chrom_hap: dict[str, dict[str, list[str]]] = {
        c: {h: [] for h in HAPLOTYPES} for c in CHROMOSOMES
    }
    intra_hap_indices: dict[str, int] = {}
    for seq in sequences_to_plot:
        row = map_by_name[seq]
        chrom, hap = row.chromosome, row.hap
        sequences_per_chrom_hap.setdefault(chrom, {}).setdefault(hap, []).append(seq)
        intra_hap_indices[seq] = len(sequences_per_chrom_hap[chrom][hap]) - 1

    # --- Pass 2: x-axis layout ---------------------------------------

    layout = _CellLayout()
    current_x = initial_x

    for chromosome in CHROMOSOMES:
        expected_haps = get_expected_haps(chromosome, sex, HAPLOTYPES, sex_determination_system)
        actual_haps_with_data = [
            h for h, seqs in sequences_per_chrom_hap.get(chromosome, {}).items() if seqs
        ]
        expected_filtered = [h for h in expected_haps if h != "unassigned"]
        combined = set(expected_filtered) | set(actual_haps_with_data)
        haps_to_draw = [h for h in HAPLOTYPES if h in combined]
        if "unassigned" in combined and "unassigned" not in haps_to_draw:
            haps_to_draw.append("unassigned")
        if not haps_to_draw:
            continue

        layout.hap_block_widths[chromosome] = {}
        layout.drawable_haps_per_chrom[chromosome] = haps_to_draw
        total_chrom_width = 0
        for hap in haps_to_draw:
            seqs = sequences_per_chrom_hap.get(chromosome, {}).get(hap, [])
            n = len(seqs)
            width_from_seqs = (2 * circle_radius) + (n - 1) * sequence_gap if n > 0 else 0
            block_width = max(width_from_seqs, min_label_width)
            layout.hap_block_widths[chromosome][hap] = block_width
            total_chrom_width += block_width
        if len(haps_to_draw) > 1:
            total_chrom_width += hap_gap * (len(haps_to_draw) - 1)

        layout.chrom_block_widths[chromosome] = total_chrom_width
        layout.chrom_start_x[chromosome] = current_x
        hap_current_x = current_x
        layout.hap_start_x[chromosome] = {}
        for hap in haps_to_draw:
            layout.hap_start_x[chromosome][hap] = hap_current_x
            hap_current_x += layout.hap_block_widths[chromosome][hap] + hap_gap
        current_x += total_chrom_width + chrom_gap

    karyotype_right_edge = current_x - chrom_gap + x_border
    layout.image_width = karyotype_right_edge + legend_band_width
    image_height = final_image_height

    # --- prepare the SVG canvas -------------------------------------

    d = draw.Drawing(layout.image_width, image_height, id_prefix="k")
    d.append(draw.Rectangle(0, 0, layout.image_width, image_height, fill=background_color))

    # --- title band (top) ------------------------------------------

    if show_title:
        title_parts: list[str] = []
        if sample_label:
            title_parts.append(sample_label)
        if database_id:
            title_parts.append(database_id)
        title_parts.append(f"{mode} view")
        if feature_set_label:
            title_parts.append(feature_set_label)
        if smoothed:
            title_parts.append("smoothed")
        title_text = "  |  ".join(title_parts)
        # Centre the title over the karyotype area (excluding the
        # legend band) so it visually balances the chromosome columns.
        d.append(
            draw.Text(
                title_text,
                14,
                karyotype_right_edge / 2,
                title_band_height - 12,
                text_anchor="middle",
                fill=text_color,
                font_weight="bold",
            )
        )

    x_coords: dict[str, float] = {}
    for chromosome, haps in layout.drawable_haps_per_chrom.items():
        for hap in haps:
            seqs = sequences_per_chrom_hap.get(chromosome, {}).get(hap, [])
            if not seqs:
                continue
            block_start_x = layout.hap_start_x[chromosome][hap]
            block_width = layout.hap_block_widths[chromosome][hap]
            width_from_seqs = (2 * circle_radius) + (len(seqs) - 1) * sequence_gap
            padding = (block_width - width_from_seqs) / 2
            for seq in seqs:
                if seq in intra_hap_indices:
                    seq_index = intra_hap_indices[seq]
                    x_coords[seq] = (
                        block_start_x + padding + circle_radius + (seq_index * sequence_gap)
                    )

    # --- Pass 3: header labels --------------------------------------

    for chromosome, start_x in layout.chrom_start_x.items():
        block_width = layout.chrom_block_widths[chromosome]
        d.append(draw.Rectangle(start_x, chrom_line_y, block_width, 1, fill=text_color))
        chromosome_text = f"{chromosome}p" if mode == "subtelomere" else chromosome
        d.append(
            draw.Text(
                chromosome_text,
                14,
                start_x + (block_width / 2),
                chrom_label_y,
                text_anchor="middle",
                fill=text_color,
            )
        )
        if len(HAPLOTYPES) > 1:
            for hap in layout.drawable_haps_per_chrom.get(chromosome, []):
                hap_start = layout.hap_start_x[chromosome][hap]
                hap_width = layout.hap_block_widths[chromosome][hap]
                hap_text = f"h{hap[3:]}" if hap.startswith("hap") else hap[:1]
                d.append(draw.Rectangle(hap_start, hap_line_y, hap_width, 1, fill=text_color))
                d.append(
                    draw.Text(
                        hap_text,
                        14,
                        hap_start + (hap_width / 2),
                        hap_label_y,
                        text_anchor="middle",
                        fill=text_color,
                    )
                )

    # --- Final pass: draw colored rectangles per feature ------------

    missing_features_warned: set[str] = set()
    features_drawn: set[str] = set()  # tracked for the legend

    def _color_for(feature: str) -> str:
        if feature in colors:
            return colors[feature]
        if feature not in missing_features_warned:
            logger.warning("feature %r not found in colors map; using white", feature)
            missing_features_warned.add(feature)
        return "#ffffff"

    for ri in inputs:
        for seq, intervals in ri.binned_bed.items():
            if seq not in seen_sequences:
                continue
            x = x_coords.get(seq)
            if x is None:
                continue
            for start, stop, feature in intervals:
                color = _color_for(feature)
                features_drawn.add(feature)
                if mode == "subtelomere":
                    seq_len = sequence_lengths.get(seq)
                    if not seq_len:
                        continue
                    # P-arm chunk
                    if seq in tel_start_sequences and start < subtelomere_boundary:
                        draw_start = max(0, start)
                        draw_stop = min(subtelomere_boundary, stop)
                        start_y = initial_y + (draw_start * pixels_per_pos)
                        stop_y = initial_y + (draw_stop * pixels_per_pos)
                        if stop_y > start_y:
                            d.append(
                                draw.Rectangle(
                                    x - circle_radius,
                                    start_y,
                                    2 * circle_radius,
                                    stop_y - start_y,
                                    fill=color,
                                )
                            )
                    # Q-arm chunk
                    if seq in tel_stop_sequences and stop > seq_len - subtelomere_boundary:
                        draw_start = max(seq_len - subtelomere_boundary, start)
                        draw_stop = min(seq_len, stop)
                        q_rel_start = draw_start - (seq_len - subtelomere_boundary)
                        q_rel_stop = draw_stop - (seq_len - subtelomere_boundary)
                        start_y = q_arm_start_y + (q_rel_start * pixels_per_pos)
                        stop_y = q_arm_start_y + (q_rel_stop * pixels_per_pos)
                        if stop_y > start_y:
                            d.append(
                                draw.Rectangle(
                                    x - circle_radius,
                                    start_y,
                                    2 * circle_radius,
                                    stop_y - start_y,
                                    fill=color,
                                )
                            )
                elif mode == "centromere":
                    cen_start, cen_stop = centromere_coords[seq]
                    draw_start = max(start, cen_start)
                    draw_stop = min(stop, cen_stop)
                    if draw_stop > draw_start:
                        y_rel_start = draw_start - cen_start
                        y_rel_stop = draw_stop - cen_start
                        start_y = initial_y + (y_rel_start * pixels_per_pos)
                        stop_y = initial_y + (y_rel_stop * pixels_per_pos)
                        d.append(
                            draw.Rectangle(
                                x - circle_radius,
                                start_y,
                                2 * circle_radius,
                                stop_y - start_y,
                                fill=color,
                            )
                        )
                else:  # genome
                    start_y = pos_to_y(start)
                    stop_y = pos_to_y(stop)
                    d.append(
                        draw.Rectangle(
                            x - circle_radius,
                            start_y,
                            2 * circle_radius,
                            stop_y - start_y,
                            fill=color,
                        )
                    )

    # --- sequence outlines (white background only) -------------------

    if background_color == "white":
        for seq, x in x_coords.items():
            if mode == "subtelomere":
                if seq in tel_start_sequences:
                    height = subtelomere_boundary * pixels_per_pos
                    d.append(
                        draw.Rectangle(
                            x - circle_radius,
                            initial_y,
                            2 * circle_radius,
                            height,
                            fill="none",
                            stroke="black",
                            stroke_width=1,
                        )
                    )
                if seq in tel_stop_sequences:
                    height = subtelomere_boundary * pixels_per_pos
                    d.append(
                        draw.Rectangle(
                            x - circle_radius,
                            q_arm_start_y,
                            2 * circle_radius,
                            height,
                            fill="none",
                            stroke="black",
                            stroke_width=1,
                        )
                    )
            elif mode == "centromere":
                cen_start, cen_stop = centromere_coords[seq]
                height = (cen_stop - cen_start) * pixels_per_pos
                d.append(
                    draw.Rectangle(
                        x - circle_radius,
                        initial_y,
                        2 * circle_radius,
                        height,
                        fill="none",
                        stroke="black",
                        stroke_width=1,
                    )
                )
            else:  # genome
                seq_len = sequence_lengths.get(seq)
                if seq_len:
                    height = pos_to_y(seq_len) - initial_y
                    d.append(
                        draw.Rectangle(
                            x - circle_radius,
                            initial_y,
                            2 * circle_radius,
                            height,
                            fill="none",
                            stroke="black",
                            stroke_width=1,
                        )
                    )

    # --- telomere indicator circles (genome mode only) ----------------

    if mode == "genome":
        telomere_color = "#006884"
        for seq, x in x_coords.items():
            if seq in tel_start_sequences:
                d.append(
                    draw.Circle(
                        x,
                        initial_y,
                        2 * circle_radius,
                        fill=telomere_color,
                        stroke="black",
                        stroke_width=1,
                    )
                )
            if seq in tel_stop_sequences:
                seq_len = sequence_lengths.get(seq)
                if seq_len:
                    q_arm_y = max(pos_to_y(seq_len), initial_y + min_q_arm_offset)
                    d.append(
                        draw.Circle(
                            x,
                            q_arm_y,
                            2 * circle_radius,
                            fill=telomere_color,
                            stroke="black",
                            stroke_width=1,
                        )
                    )

    # --- scale bar ---------------------------------------------------

    scale_bar_pixel_height = mode_params.scale_bar_length * pixels_per_pos
    scale_bar_x = 40
    scale_bar_width = 2
    label1_x = 35
    label2_x = 20
    pos_per_pixel = 1 / pixels_per_pos
    if pos_per_pixel >= 1000:
        resolution_label_text = f"1px = {pos_per_pixel / 1000:.1f} kbp".replace(".0", "")
    else:
        resolution_label_text = f"1px = {int(pos_per_pixel)} bp"
    label_y_center = initial_y + (scale_bar_pixel_height / 2)
    d.append(
        draw.Rectangle(
            scale_bar_x, initial_y, scale_bar_width, scale_bar_pixel_height, fill=text_color
        )
    )
    d.append(
        draw.Text(
            mode_params.scale_bar_label,
            10,
            label1_x,
            label_y_center,
            text_anchor="middle",
            fill=text_color,
            transform=f"rotate(-90, {label1_x}, {label_y_center})",
        )
    )
    d.append(
        draw.Text(
            resolution_label_text,
            10,
            label2_x,
            label_y_center,
            text_anchor="middle",
            fill=text_color,
            transform=f"rotate(-90, {label2_x}, {label_y_center})",
        )
    )
    if mode == "subtelomere":
        q_arm_end_y = q_arm_start_y + q_arm_height
        scale_bar_start_y = q_arm_end_y - scale_bar_pixel_height
        q_arm_label_y_center = scale_bar_start_y + (scale_bar_pixel_height / 2)
        d.append(
            draw.Rectangle(
                scale_bar_x,
                scale_bar_start_y,
                scale_bar_width,
                scale_bar_pixel_height,
                fill=text_color,
            )
        )
        d.append(
            draw.Text(
                mode_params.scale_bar_label,
                10,
                label1_x,
                q_arm_label_y_center,
                text_anchor="middle",
                fill=text_color,
                transform=f"rotate(-90, {label1_x}, {q_arm_label_y_center})",
            )
        )
        d.append(
            draw.Text(
                resolution_label_text,
                10,
                label2_x,
                q_arm_label_y_center,
                text_anchor="middle",
                fill=text_color,
                transform=f"rotate(-90, {label2_x}, {q_arm_label_y_center})",
            )
        )

    # --- footer labels (subtelomere mode only) -----------------------

    if mode == "subtelomere":
        drawing_end_y = q_arm_start_y + (subtelomere_boundary * pixels_per_pos)
        for chromosome, start_x in layout.chrom_start_x.items():
            block_width = layout.chrom_block_widths[chromosome]
            chromosome_text = f"{chromosome}q"
            if len(HAPLOTYPES) == 1:
                bottom_chrom_line_y = drawing_end_y + 12
                bottom_chrom_label_y = drawing_end_y + 27
                d.append(
                    draw.Rectangle(start_x, bottom_chrom_line_y, block_width, 1, fill=text_color)
                )
                d.append(
                    draw.Text(
                        chromosome_text,
                        14,
                        start_x + (block_width / 2),
                        bottom_chrom_label_y,
                        text_anchor="middle",
                        fill=text_color,
                    )
                )
            else:
                bottom_hap_line_y = drawing_end_y + 12
                bottom_hap_label_y = drawing_end_y + 27
                bottom_chrom_line_y = drawing_end_y + 32
                bottom_chrom_label_y = drawing_end_y + 47
                for hap in layout.drawable_haps_per_chrom.get(chromosome, []):
                    hap_start = layout.hap_start_x[chromosome][hap]
                    hap_width = layout.hap_block_widths[chromosome][hap]
                    hap_text = f"h{hap[3:]}" if hap.startswith("hap") else hap[:1]
                    d.append(
                        draw.Rectangle(hap_start, bottom_hap_line_y, hap_width, 1, fill=text_color)
                    )
                    d.append(
                        draw.Text(
                            hap_text,
                            14,
                            hap_start + (hap_width / 2),
                            bottom_hap_label_y,
                            text_anchor="middle",
                            fill=text_color,
                        )
                    )
                d.append(
                    draw.Rectangle(start_x, bottom_chrom_line_y, block_width, 1, fill=text_color)
                )
                d.append(
                    draw.Text(
                        chromosome_text,
                        14,
                        start_x + (block_width / 2),
                        bottom_chrom_label_y,
                        text_anchor="middle",
                        fill=text_color,
                    )
                )

    # --- legend (right margin) --------------------------------------

    if show_legend and features_drawn:
        # One row per feature actually rendered (small font, tight
        # rows). Sorted with a stable convention: numbered chromosome
        # leaves naturally (chr1, chr2, ...), then everything else
        # alphabetical. Falls back cleanly when the feature names
        # don't match the chromosome pattern.
        def _legend_sort_key(name: str) -> tuple[int, int, str]:
            if name.startswith("chr"):
                suffix = name[3:]
                if suffix.isdigit():
                    return (0, int(suffix), "")
                return (1, 0, suffix)
            return (2, 0, name)

        legend_x = karyotype_right_edge + legend_pad_x
        legend_y = initial_y
        for i, feature in enumerate(sorted(features_drawn, key=_legend_sort_key)):
            row_y = legend_y + i * legend_row_height
            # Bail out if the legend would overflow the SVG height
            # (rare on tall karyotypes but possible for small ones).
            if row_y + legend_swatch_size > image_height - 5:
                logger.info(
                    "legend truncated at %d entries; remaining features fit outside the SVG height",
                    i,
                )
                break
            d.append(
                draw.Rectangle(
                    legend_x,
                    row_y,
                    legend_swatch_size,
                    legend_swatch_size,
                    fill=_color_for(feature),
                    stroke="black",
                    stroke_width=0.5,
                )
            )
            d.append(
                draw.Text(
                    feature,
                    legend_text_size,
                    legend_x + legend_swatch_size + 4,
                    row_y + legend_swatch_size - 2,
                    text_anchor="start",
                    fill=text_color,
                )
            )

    d.save_svg(str(output_path))


def convert_svg(svg_path: Path, target_path: Path) -> None:
    """Convert an existing SVG file to another format (PDF / PNG).

    Format is inferred from ``target_path``'s extension. PDF / PNG
    use cairosvg, which in turn requires the native ``libcairo``
    library at runtime; if cairo isn't installed in the active
    environment, this function raises :class:`KaryotypeError` with
    an actionable install hint. SVG-to-SVG is a plain file copy and
    has no dependency on cairo.
    """
    ext = target_path.suffix.lower()
    if ext == ".svg":
        if svg_path.resolve() != target_path.resolve():
            target_path.write_bytes(svg_path.read_bytes())
        return
    if ext not in (".pdf", ".png"):
        raise KaryotypeError(f"unsupported output format {ext!r}; expected one of .svg, .pdf, .png")

    try:
        import cairosvg
    except OSError as e:
        # cairosvg loads libcairo at import time; missing libcairo
        # surfaces as OSError ("no library called 'cairo'... was found").
        raise KaryotypeError(
            f"cannot convert to {ext}: native libcairo is not installed. "
            "Install in the active conda env with `conda install -c conda-forge cairo`, "
            "or use `--format svg` only."
        ) from e

    svg_bytes = svg_path.read_bytes()
    if ext == ".pdf":
        cairosvg.svg2pdf(bytestring=svg_bytes, write_to=str(target_path))
    else:  # .png
        cairosvg.svg2png(bytestring=svg_bytes, write_to=str(target_path))
