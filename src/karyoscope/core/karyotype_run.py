"""Orchestrator for ``karyoscope karyotype``.

Glues the algorithm in :mod:`karyoscope.core.karyotype` to the rest of
the pipeline. Auto-derives missing prerequisites the same way
``scaffold`` and ``centromeres`` do:

1. Resolve the database, parse its manifest + colors.
2. For each ``-i`` input, make sure scaffolded BEDs exist for the
   requested feature sets (cascade through :func:`scaffold_run`).
3. Bin each per-input scaffolded BED at the mode-appropriate bin
   size if not already binned.
4. For centromere mode, ensure centromere coordinates exist for each
   input (cascade through :func:`centromeres_run`).
5. Load the per-input binned BEDs + map rows + (centromere mode)
   centromere coordinates into :class:`RenderInput` records.
6. Per requested feature set, render one SVG.

Output paths follow the convention
``<base>.<dbid>.<mode>.<feature_set>.karyotype.svg`` per feature set.
``<base>`` is either an explicit ``--output PATH`` (then we use the
basename without ``.svg``) or, when no explicit path is given, a
sample-identifying base derived from the input stems: the single
input's stem, or for multi-input runs the longest common prefix of the
input stems (e.g. ``GM04890.haplotype1`` + ``GM04890.haplotype2`` ->
``GM04890``), falling back to the first stem when the stems share no
separator-delimited prefix. See :func:`_common_base`.
"""

from __future__ import annotations

import gzip
import logging
import shutil
import tempfile
import time
from collections import OrderedDict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from karyoscope.core.annotate import annotate_batch
from karyoscope.core.bin import bin_features, leaves_for
from karyoscope.core.centromeres import (
    DEFAULT_COARSE_BIN_SIZE,
    DEFAULT_FINE_BIN_SIZE,
    _resolve_centromere_role,
    centromeres_run,
    find_centromere_ranges,
)
from karyoscope.core.io.bgzip import bgzip_file
from karyoscope.core.io.colors import (
    colors_for_set,
    parse_colors_and_groups,
    validate_colors,
    validate_legend_groups,
)
from karyoscope.core.io.hierarchy import parse_hierarchy
from karyoscope.core.io.scaffold_map import MapRow, map_signature, read_map
from karyoscope.core.karyotype import (
    DEFAULT_HUMAN_CHROMOSOMES,
    RenderInput,
    convert_svg,
    render_karyotype,
)
from karyoscope.core.scaffold import (
    DEFAULT_HUMAN_ACROCENTRICS,
    DEFAULT_MIN_SCAFFOLD_LENGTH,
    DEFAULT_SCAFFOLD_GAP_SIZE,
    Interval,
    combined_map_rows,
    rewrite_bed,
)
from karyoscope.core.scaffold_run import InputSpec, _resolve_roles, scaffold_run
from karyoscope.exceptions import KaryotypeError
from karyoscope.installed import resolve_database
from karyoscope.manifest import validate_database_layout
from karyoscope.progress import SILENT, Progress

logger = logging.getLogger(__name__)


#: Default bin size (bp) per render mode. The archive's per-mode
#: pixel scales are tuned to these defaults; deviating substantially
#: produces under- or over-detailed SVGs.
DEFAULT_BIN_SIZE_BY_MODE: dict[str, int] = {
    "genome": 1_000_000,
    "centromere": 100_000,
    "subtelomere": 100,
}

#: Target number of render bins across the longest sequence when the genome-view
#: bin size is chosen automatically (no ``--bin-size`` override). Human 250 Mb /
#: 250 -> ~1 Mb (the old default); Arabidopsis ~32 Mb / 250 -> ~130 kb, which
#: restores feature diversity that 1 Mb bins wash out via the plurality-per-bin
#: rule. Only ``genome`` scales; ``centromere``/``subtelomere`` keep fixed
#: fine-grained bins (scaling them to the whole-genome target would coarsen them).
_GENOME_TARGET_BIN_COUNT = 250
_MIN_AUTO_BIN_SIZE = 10_000


def _scaffolding_prereq_note(
    *,
    requested: list[str],
    scaffold_manifest_roles: dict[str, str],
    scaffold_available: list[str],
    centromere_fs: str | None,
    scaffold_db_id: str | None,
) -> str:
    """Explain feature sets the cascade needs but the user didn't ask to plot.

    Laying out a karyotype needs more than the sets being drawn: the
    chromosome-assignment set says which chromosome each contig belongs to,
    and the region-assignment set orients it (and, in centromere mode,
    locates the centromere). :func:`scaffold_run` therefore appends those
    role sets to whatever ``--feature-set`` asked for.

    Without this note the progress output is quietly baffling — asking for
    two feature sets and watching ``annotate`` report three, with no hint
    which extra one appeared or why. Returns ``""`` when the requested sets
    already cover the roles, so the common case stays uncluttered.

    Deliberately phrased as a requirement rather than an action: the BEDs
    may already exist from an earlier run, in which case the cascade
    short-circuits and nothing is annotated at all.
    """
    try:
        chromosome_fs, region_fs = _resolve_roles(scaffold_manifest_roles, scaffold_available)
    except Exception:
        # Role resolution failing is reported properly by the cascade; a
        # missing progress line is not worth turning into an error here.
        return ""

    already = set(requested)
    extras: list[str] = []
    for fs in (chromosome_fs, region_fs, centromere_fs):
        if fs is not None and fs not in already and fs not in extras:
            extras.append(fs)
    if not extras:
        return ""

    names = ", ".join(extras)
    source = f" from {scaffold_db_id}" if scaffold_db_id else ""
    return f"scaffolding also needs {names}{source} (annotated if missing; not rendered)"


def _auto_bin_size(mode: str, max_seq_len: int) -> int:
    """Genome-view bin size scaled to the longest sequence; fixed otherwise."""
    if mode != "genome" or max_seq_len <= 0:
        return DEFAULT_BIN_SIZE_BY_MODE[mode]
    return max(_MIN_AUTO_BIN_SIZE, round(max_seq_len / _GENOME_TARGET_BIN_COUNT))


#: All available render modes, in canonical order. Used when the
#: user doesn't specify ``--mode`` -- the default is to render every
#: mode for every requested feature set.
ALL_MODES: tuple[str, ...] = ("genome", "centromere", "subtelomere")


_FASTA_EXTS: tuple[str, ...] = (
    ".fasta.gz",
    ".fa.gz",
    ".fna.gz",
    ".fasta",
    ".fa",
    ".fna",
)


def _input_stem(path: Path) -> str:
    name = path.name
    lower = name.lower()
    for ext in _FASTA_EXTS:
        if lower.endswith(ext):
            return name[: -len(ext)]
    return path.stem


#: Separator characters used to find a token boundary in input stems.
_STEM_SEPARATORS = "._-"


def _common_base(stems: list[str]) -> str:
    """The filename base shared by one or more input stems.

    For a single input this is just its stem. For multi-input runs --
    e.g. separate ``hap1``/``hap2`` FASTAs -- the output filename should
    name the *sample*, not just the first input, so it isn't misleading
    (a both-hap plot named ``...haplotype1...``). We take the
    character-wise longest common prefix of the stems, trim it back to
    the last separator (so a partially-shared trailing token like
    ``haplotype`` is dropped rather than left dangling), and strip
    trailing separators::

        ["GM04890.haplotype1", "GM04890.haplotype2"]  -> "GM04890"
        ["BJ.hap1", "BJ.hap2"]                         -> "BJ"
        ["HG002.maternal", "HG002.paternal"]           -> "HG002"

    Falls back to the first stem when the stems are identical or share
    no separator-delimited prefix (so the result is never empty).
    """
    if not stems:
        return ""
    first = stems[0]
    if len(stems) == 1 or all(s == first for s in stems):
        return first

    lcp_len = len(first)
    for s in stems[1:]:
        limit = min(lcp_len, len(s))
        i = 0
        while i < limit and first[i] == s[i]:
            i += 1
        lcp_len = i
    lcp = first[:lcp_len]

    cut = max((lcp.rfind(c) for c in _STEM_SEPARATORS), default=-1)
    if cut == -1:
        return first  # no shared separator-delimited prefix
    base = lcp[:cut].rstrip(_STEM_SEPARATORS)
    return base or first


def _scaffolded_bed_path(
    out_dir: Path, stem: str, db_id: str, fs: str, variant: str = "smoothed"
) -> Path:
    gz = out_dir / f"{stem}.{db_id}.{fs}.{variant}.scaffolded.bed.gz"
    if gz.is_file():
        return gz
    plain = out_dir / f"{stem}.{db_id}.{fs}.{variant}.scaffolded.bed"
    if plain.is_file():
        return plain
    return gz


def _annotation_bed_path(
    out_dir: Path, stem: str, db_id: str, fs: str, variant: str = "smoothed"
) -> Path:
    """The annotate-produced BED path (smoothed or presmoothed).

    Used on the ``--no-scaffolding`` codepath where scaffolding skipped
    writing per-FS scaffolded BEDs. We bin this file and apply the
    scaffold map at bin time via :func:`rewrite_bed`.
    """
    gz = out_dir / f"{stem}.{db_id}.{fs}.{variant}.bed.gz"
    if gz.is_file():
        return gz
    plain = out_dir / f"{stem}.{db_id}.{fs}.{variant}.bed"
    if plain.is_file():
        return plain
    return gz


def _binned_scaffolded_bed_path(
    out_dir: Path,
    stem: str,
    db_id: str,
    fs: str,
    bin_size: int,
    variant: str = "smoothed",
) -> Path:
    return out_dir / f"{stem}.{db_id}.{fs}.{variant}.scaffolded.binned{bin_size}.bed.gz"


#: Filename tag distinguishing combined-chromosome outputs from the
#: per-contig scaffolded outputs. Must match
#: ``scaffold_run._COMBINED_TAG`` so karyotype reads the combined BEDs
#: that scaffold writes.
_COMBINED_TAG = "combined_chromosomes"


def _colors_filename_tag(colors_path: Path | None) -> str:
    """Output-filename tag for a custom ``--colors`` file (empty for default colours).

    Appending the colour file's stem (e.g. ``.colors_chromosome``) keeps a
    custom-colour render from overwriting the default-colour one in the same
    directory.
    """
    return f".{colors_path.stem}" if colors_path is not None else ""


def _combined_scaffolded_bed_path(
    out_dir: Path, stem: str, db_id: str, fs: str, variant: str = "smoothed"
) -> Path:
    """The combined-chromosome scaffolded BED scaffold writes (mode 'both')."""
    gz = out_dir / f"{stem}.{db_id}.{fs}.{variant}.scaffolded.{_COMBINED_TAG}.bed.gz"
    if gz.is_file():
        return gz
    plain = out_dir / f"{stem}.{db_id}.{fs}.{variant}.scaffolded.{_COMBINED_TAG}.bed"
    if plain.is_file():
        return plain
    return gz


def _binned_combined_scaffolded_bed_path(
    out_dir: Path,
    stem: str,
    db_id: str,
    fs: str,
    bin_size: int,
    variant: str = "smoothed",
) -> Path:
    return (
        out_dir
        / f"{stem}.{db_id}.{fs}.{variant}.scaffolded.{_COMBINED_TAG}.binned{bin_size}.bed.gz"
    )


def _binned_mapsig_path(binned: Path) -> Path:
    """Sidecar path recording the scaffold-map signature a binned BED was built from."""
    return binned.with_name(binned.name + ".mapsig")


def _binned_bed_is_current(binned: Path, map_rows: list[MapRow] | None) -> bool:
    """True if ``binned`` was built from the same scaffold map as ``map_rows``.

    The binned-scaffolded BED bakes in the scaffold map's rename and
    orientation (its sequence names are the encoded
    ``<chrom>_<hap>_<contig>`` names, and coordinates may be flipped).
    Reusing it after the map changed -- e.g. after hap inference was
    corrected so a contig moved from ``hap1`` to ``hap2`` -- would serve
    a stale haplotype layout. We compare the current map's signature
    against the one recorded in the ``.mapsig`` sidecar when the binned
    BED was written.

    Fails closed: a missing or unreadable sidecar (e.g. a binned BED
    from before this guard existed) is treated as not-current, so it is
    rebuilt once. ``map_rows is None`` means the caller cannot supply a
    map to check against; we preserve the legacy reuse-if-present
    behaviour in that case rather than force an un-checkable rebuild.
    """
    if map_rows is None:
        return True
    try:
        recorded = _binned_mapsig_path(binned).read_text().strip()
    except OSError:
        return False
    return recorded == map_signature(map_rows)


def _write_binned_mapsig(binned: Path, map_rows: list[MapRow] | None) -> None:
    """Record the scaffold-map signature next to a freshly built binned BED."""
    if map_rows is None:
        return
    _binned_mapsig_path(binned).write_text(map_signature(map_rows) + "\n")


def _load_binned_bed(path: Path) -> OrderedDict[str, list[Interval]]:
    opener = gzip.open if path.suffix == ".gz" else open
    out: OrderedDict[str, list[Interval]] = OrderedDict()
    with opener(path, "rt") as h:
        for raw in h:
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            out.setdefault(parts[0], []).append((start, end, parts[3]))
    return out


def _load_centromeres_bed(path: Path) -> dict[str, tuple[int, int]]:
    opener = gzip.open if path.suffix == ".gz" else open
    out: dict[str, tuple[int, int]] = {}
    with opener(path, "rt") as h:
        for raw in h:
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            out[parts[0]] = (start, end)
    return out


def _ensure_binned_scaffolded(
    *,
    out_dir: Path,
    stem: str,
    db_id: str,
    fs: str,
    bin_size: int,
    leaf_set: set[str],
    auto: bool,
    input_name: str,
    threads: int,
    map_rows: list[MapRow] | None = None,
    variant: str = "smoothed",
    combined: bool = False,
) -> Path:
    """Return the binned scaffolded BED path, building it if missing or stale.

    When ``combined`` is True, operates on the combined-chromosome
    scaffolded BED (sequence names ``<chrom>_<hap>``) that
    ``scaffold_run(combine_chromosomes=True, mode="both")`` always
    materialises. The combined BED is binned directly; there is no
    ``--no-scaffolding`` annotation fallback on this path, because the
    combined coordinate transform only exists in the full-resolution
    combined BED. ``map_rows`` here should be the *combined* map rows so
    the ``.mapsig`` guard invalidates when the underlying scaffold map
    changes.

    An existing binned BED is reused only when it is also *current* with
    respect to the scaffold map (see :func:`_binned_bed_is_current`): the
    binned BED bakes the map's rename + orientation into its contents, so
    a map that changed since it was built (e.g. corrected hap inference
    moving a contig to a different haplotype) forces a rebuild rather
    than serving a stale layout. Each freshly built binned BED records
    the map signature it was built from in a ``.mapsig`` sidecar.

    Two construction paths:

    * **Scaffolded BED on disk** (the historical path): bin it directly.
    * **Scaffolded BED missing but a map is available** (the
      ``--no-scaffolding`` path): bin the annotation (smoothed or
      presmoothed) BED at the requested ``bin_size``, then stream the
      binned output through :func:`rewrite_bed` to apply the map
      (rename contigs + mirror coordinates for flipped contigs).

    ``variant`` selects whether we read the smoothed or presmoothed
    annotation BEDs and produce correspondingly named intermediates.
    """
    if combined:
        out = _binned_combined_scaffolded_bed_path(
            out_dir, stem, db_id, fs, bin_size, variant=variant
        )
    else:
        out = _binned_scaffolded_bed_path(out_dir, stem, db_id, fs, bin_size, variant=variant)
    if out.is_file() and _binned_bed_is_current(out, map_rows):
        return out
    if not auto:
        # File absent, or present but built from a superseded scaffold
        # map (e.g. hap inference was corrected since). We were told not
        # to auto-derive, so we can't fix it here -- but silently
        # serving a stale haplotype layout is exactly the bug this guard
        # exists to prevent, so fail loudly instead.
        stale = out.is_file()
        raise KaryotypeError(
            f"{'stale' if stale else 'missing'} binned scaffolded BED for "
            f"{input_name}, feature set {fs!r}, bin size {bin_size} "
            f"(at {out}). "
            + ("The scaffold map changed since it was built. " if stale else "")
            + "Re-run with auto-derive enabled"
            + (f" (or delete {out.name} and its .mapsig)." if stale else ".")
        )
    if combined:
        # The combined coordinate transform lives only in the
        # full-resolution combined BED, which scaffold always writes in
        # combine mode. Bin it directly; there is no annotation fallback.
        combined_src = _combined_scaffolded_bed_path(out_dir, stem, db_id, fs, variant=variant)
        if not combined_src.is_file():
            raise KaryotypeError(
                f"cannot bin combined {fs!r} for {input_name}: combined "
                f"scaffolded BED missing at {combined_src} (scaffold should "
                f"have produced it in combine mode)"
            )
        bin_features(
            combined_src,
            out,
            bin_size=bin_size,
            leaf_set=leaf_set or None,
            threads=threads,
        )
        _write_binned_mapsig(out, map_rows)
        return out
    scaffolded_src = _scaffolded_bed_path(out_dir, stem, db_id, fs, variant=variant)
    if scaffolded_src.is_file():
        bin_features(
            scaffolded_src,
            out,
            bin_size=bin_size,
            leaf_set=leaf_set or None,
            threads=threads,
        )
        _write_binned_mapsig(out, map_rows)
        return out

    # Fallback: bin the annotation BED, then apply the scaffold map.
    if map_rows is None:
        raise KaryotypeError(
            f"cannot bin {fs!r} for {input_name}: scaffolded BED missing at "
            f"{scaffolded_src} and no scaffold map provided for post-bin "
            f"renaming."
        )
    annotation_src = _annotation_bed_path(out_dir, stem, db_id, fs, variant=variant)
    if not annotation_src.is_file():
        raise KaryotypeError(
            f"cannot bin {fs!r} for {input_name}: {variant} BED missing at "
            f"{annotation_src} (and scaffolded BED also missing)"
        )
    tmpdir = Path(tempfile.mkdtemp(prefix="ks_karyo_bin_", dir=out_dir))
    try:
        tmp_binned = tmpdir / "binned.bed.gz"
        bin_features(
            annotation_src,
            tmp_binned,
            bin_size=bin_size,
            leaf_set=leaf_set or None,
            threads=threads,
        )
        rewrite_bed(tmp_binned, out, map_rows=map_rows)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    _write_binned_mapsig(out, map_rows)
    return out


def _centromeres_bed_path(out_dir: Path, stem: str, db_id: str) -> Path:
    gz = out_dir / f"{stem}.{db_id}.centromeres.bed.gz"
    if gz.is_file():
        return gz
    plain = out_dir / f"{stem}.{db_id}.centromeres.bed"
    if plain.is_file():
        return plain
    return gz


def _combined_centromeres_bed_path(out_dir: Path, stem: str, db_id: str) -> Path:
    gz = out_dir / f"{stem}.{db_id}.centromeres.{_COMBINED_TAG}.bed.gz"
    if gz.is_file():
        return gz
    plain = out_dir / f"{stem}.{db_id}.centromeres.{_COMBINED_TAG}.bed"
    if plain.is_file():
        return plain
    return gz


def _ensure_combined_centromeres(
    *,
    out_dir: Path,
    stem: str,
    db_id: str,
    centromere_fs: str,
    centromere_leaves: set[str],
    map_rows: list[MapRow],
    auto: bool,
    input_name: str,
    threads: int,
    bgzip: bool,
    variant: str = "smoothed",
) -> dict[str, tuple[int, int]]:
    """Detect centromere ranges in the combined coordinate system.

    The non-combined cascade runs ``centromeres_run`` per contig; on the
    combine path we detect against the combined-chromosome BEDs instead,
    so the ranges are keyed by ``<chrom>_<hap>`` and expressed in
    combined coordinates -- exactly what the centromere-mode renderer
    needs to join against the combined binned BED.

    Reuses :func:`find_centromere_ranges` and the combined binnings
    (coarse 1 Mb + fine 100 kb) that :func:`_ensure_binned_scaffolded`
    produces and caches. Writes the result to a
    ``centromeres.combined_chromosomes.bed`` for parity with the
    per-contig output, then returns it in memory.
    """
    coarse_path = _ensure_binned_scaffolded(
        out_dir=out_dir,
        stem=stem,
        db_id=db_id,
        fs=centromere_fs,
        bin_size=DEFAULT_COARSE_BIN_SIZE,
        leaf_set=centromere_leaves,
        auto=auto,
        input_name=input_name,
        threads=threads,
        map_rows=map_rows,
        variant=variant,
        combined=True,
    )
    coarse_bins = _load_binned_bed(coarse_path)
    fine_bins: OrderedDict[str, list[Interval]] | None = None
    if DEFAULT_FINE_BIN_SIZE:
        fine_path = _ensure_binned_scaffolded(
            out_dir=out_dir,
            stem=stem,
            db_id=db_id,
            fs=centromere_fs,
            bin_size=DEFAULT_FINE_BIN_SIZE,
            leaf_set=centromere_leaves,
            auto=auto,
            input_name=input_name,
            threads=threads,
            map_rows=map_rows,
            variant=variant,
            combined=True,
        )
        fine_bins = _load_binned_bed(fine_path)
    ranges = find_centromere_ranges(coarse_bins, fine_bins)

    plain = out_dir / f"{stem}.{db_id}.centromeres.{_COMBINED_TAG}.bed"
    with plain.open("w") as h:
        for contig, (cstart, cend) in ranges.items():
            h.write(f"{contig}\t{cstart}\t{cend}\n")
    if bgzip:
        bgzip_file(plain, threads=threads)
    return dict(ranges)


#: Output formats supported by :func:`karyotype_run`. SVG is the
#: native drawsvg output; PDF and PNG are produced by converting the
#: SVG via :mod:`cairosvg`.
SUPPORTED_FORMATS: tuple[str, ...] = ("svg", "pdf", "png")


@dataclass
class KaryotypeResult:
    """One render produced by :func:`karyotype_run`.

    ``output_paths`` is the list of files written for this
    (mode, feature_set) combination -- one per requested format. The
    first entry is always the SVG (the canonical output); other
    formats are derived from it.
    """

    feature_set: str
    mode: str
    output_paths: list[Path]

    @property
    def svg_path(self) -> Path:
        """The SVG output (always the first entry)."""
        return self.output_paths[0]


def karyotype_run(
    inputs: list[InputSpec],
    *,
    db_root: Path,
    db_id: str | None = None,
    scaffold_db_id: str | None = None,
    scaffold_db_root: Path | None = None,
    feature_sets: list[str] | None = None,
    modes: list[str] | None = None,
    sex: str | None = None,
    sex_determination_system: str = "XY",
    background_color: str = "white",
    bin_size: int | None = None,
    pixels_per_mb: float | None = None,
    subtelomere_boundary: int = 250_000,
    min_scaffold_length: int = DEFAULT_MIN_SCAFFOLD_LENGTH,
    telo_motif: str | None = None,
    acrocentrics: set[str] | None = None,
    combine_chromosomes: bool = False,
    scaffold_gap_size: int = DEFAULT_SCAFFOLD_GAP_SIZE,
    combine_acrocentrics: bool = False,
    split_haps_regex: str | None = None,
    threads: int = 0,
    auto: bool = True,
    bgzip: bool = True,
    scaffolding: bool = True,
    annotation_variant: str = "smoothed",
    output_dir: Path | None = None,
    output_path: Path | None = None,
    seed_human_chromosomes: bool = True,
    formats: list[str] | None = None,
    progress: Progress = SILENT,
    sample_label: str | None = None,
    show_title: bool = True,
    show_legend: bool = True,
    colors_path: Path | None = None,
) -> list[KaryotypeResult]:
    """Render one SVG per (mode, feature_set) combination.

    With no overrides, the default behaviour is to render every
    available mode (``genome``, ``centromere``, ``subtelomere``) for
    every feature set declared in the database's manifest. Either
    axis can be restricted via ``modes=`` or ``feature_sets=``; the
    Cartesian product of what's specified is what gets rendered.

    Most parameters mirror the CLI flags; see
    :mod:`karyoscope.commands.karyotype` for the user-facing
    documentation.
    """
    if not inputs:
        raise KaryotypeError("at least one --input is required")

    t_karyo_start = time.perf_counter()
    requested_modes: list[str] = list(modes) if modes else list(ALL_MODES)
    unknown_modes = [m for m in requested_modes if m not in ALL_MODES]
    if unknown_modes:
        raise KaryotypeError(f"unknown mode(s) {unknown_modes!r}; expected from {list(ALL_MODES)}")

    # Output formats: default to just SVG. Preserve user-supplied
    # order so the listed formats appear in result.output_paths in
    # the same order on disk (cosmetic). SVG is always produced
    # first (everything else is derived from it).
    requested_formats: list[str] = [f.lower() for f in formats] if formats else ["svg"]
    unknown_formats = [f for f in requested_formats if f not in SUPPORTED_FORMATS]
    if unknown_formats:
        raise KaryotypeError(
            f"unsupported format(s) {unknown_formats!r}; expected from {list(SUPPORTED_FORMATS)}"
        )
    if "svg" not in requested_formats:
        # We always need an SVG (it's the source for PDF/PNG); insert
        # at the front so the file is written first.
        requested_formats = ["svg", *requested_formats]

    db_id_resolved, db_dir = resolve_database(db_root, db_id)
    manifest = validate_database_layout(db_dir)
    available = list(manifest.feature_sets)

    # Layout / scaffold database. When --scaffold-db is given, chromosome
    # ordering + region orientation + centromere detection are derived from
    # THIS database (which must carry those roles), while the feature set(s)
    # in --feature-set are plotted + coloured from --db. This lets a
    # plot-only database (e.g. a cytoband database with no chromosome/region
    # feature sets) borrow the layout from a roles-bearing database such as
    # KS_human_CHM13_v2. Default (unset): the --db database supplies both
    # layout and plotting -- the original single-database behaviour.
    use_scaffold_db = scaffold_db_id is not None and scaffold_db_id != db_id_resolved
    if use_scaffold_db:
        if combine_chromosomes:
            raise KaryotypeError(
                "--scaffold-db cannot be combined with --combine-chromosomes "
                "(combined-chromosome layout requires the plotted feature set's "
                "scaffolded BEDs in the layout database)."
            )
        scaffold_db_root_eff = scaffold_db_root if scaffold_db_root is not None else db_root
        scaffold_db_id_resolved, scaffold_db_dir = resolve_database(
            scaffold_db_root_eff, scaffold_db_id
        )
        scaffold_manifest = validate_database_layout(scaffold_db_dir)
        scaffold_available = list(scaffold_manifest.feature_sets)
        # Fail early with a clear message if the layout DB lacks the
        # chromosome/region roles -- it can't drive scaffold otherwise.
        _resolve_roles(scaffold_manifest.roles, scaffold_available)
    else:
        scaffold_db_root_eff = db_root
        scaffold_db_id_resolved = db_id_resolved
        scaffold_manifest = manifest
        scaffold_available = available

    requested: list[str] = list(feature_sets) if feature_sets else list(available)
    unknown = [fs for fs in requested if fs not in available]
    if unknown:
        raise KaryotypeError(
            f"requested feature set(s) {unknown!r} not declared in manifest "
            f"(available: {available!r})"
        )
    if not requested:
        raise KaryotypeError("no feature sets to render; manifest is empty?")

    # bin_size override only makes sense with a single mode (different
    # modes have different natural bin sizes).
    if bin_size is not None:
        if len(requested_modes) != 1:
            raise KaryotypeError(
                "--bin-size can only be set when exactly one --mode is "
                f"requested (got modes={requested_modes!r}). Drop --bin-size "
                "or restrict to one mode."
            )
        if bin_size < 1:
            raise KaryotypeError(f"--bin-size must be a positive integer, got {bin_size}")

    # Colours: a user-supplied --colors file overrides the database default.
    # Its stem tags the output filename so a custom-colour render never clobbers
    # the default-colour one (e.g. '...smoothed.colors_chromosome.karyotype.svg').
    colors_source = colors_path if colors_path is not None else db_dir / manifest.colors
    colors, legend_groups_all = parse_colors_and_groups(colors_source)
    colors_tag = _colors_filename_tag(colors_path)
    hierarchy = parse_hierarchy(db_dir / manifest.hierarchy)
    # Chromosomes to seed the karyotype layout with (so one missing from the
    # sample still gets an empty column): the chromosome-role feature set's
    # leaves. Replaces the old hardcoded human list; identical for the human DB
    # (chr1..chrY), and organism-correct elsewhere. Empty if unresolvable ->
    # only chromosomes present in the data are drawn.
    chromosome_role_fs = _resolve_roles(scaffold_manifest.roles, scaffold_available)[0]
    expected_chromosomes = list(leaves_for(hierarchy, chromosome_role_fs))

    # Hard check that every hierarchy node has a colour. ``download``
    # runs the same check at install time; we re-run it here so
    # locally-installed databases (which bypass download validation)
    # still get caught before any SVG is rendered.
    color_issues = validate_colors(hierarchy, colors)
    color_issues += validate_legend_groups(colors, legend_groups_all)
    if color_issues:
        raise KaryotypeError(
            f"database {db_id_resolved!r} is missing colour entries; refusing "
            "to render. Fix colors.tsv (or use `karyoscope info` to see the "
            "full list of issues):\n  - " + "\n  - ".join(color_issues)
        )

    if combine_chromosomes and scaffold_gap_size < 0:
        raise KaryotypeError(f"scaffold_gap_size must be >= 0, got {scaffold_gap_size}")

    # The acrocentric set drives which (chromosome, hap) groups get
    # combined. Resolve it the same way scaffold does so the synthetic
    # combined map rows agree with the combined BED's object names.
    acros_set = acrocentrics if acrocentrics is not None else set(DEFAULT_HUMAN_ACROCENTRICS)

    # Centromere detection feature set (combine path only -- the
    # non-combine path lets centromeres_run resolve it internally).
    want_centromere = "centromere" in requested_modes
    centromere_fs: str | None = None
    centromere_leaves: set[str] = set()
    if combine_chromosomes and want_centromere:
        centromere_fs = _resolve_centromere_role(manifest.roles, available)
        centromere_leaves = leaves_for(hierarchy, centromere_fs)

    # Announce the run before the cascade starts. Everything from here --
    # annotate, scaffold, bin, render -- can take tens of minutes, and the
    # nested annotate reports through a child of this reporter, so the user
    # sees the expensive middle of the pipeline rather than a blank terminal.
    #
    # Placed after role resolution rather than at the top of the function so
    # the announcement can name the feature sets the cascade will pull in
    # beyond the ones being rendered. Nothing above this point is expensive
    # (manifest, hierarchy, and colour validation), so nothing is hidden by
    # announcing here.
    sample_names = ", ".join(spec.path.name for spec in inputs)
    progress.start(
        f"Rendering karyotypes for {sample_names} against {db_id_resolved}",
        f"{len(requested_modes)} mode(s) x {len(requested)} feature set(s) "
        f"= {len(requested_modes) * len(requested)} render(s)",
        _scaffolding_prereq_note(
            requested=requested,
            scaffold_manifest_roles=scaffold_manifest.roles,
            scaffold_available=scaffold_available,
            centromere_fs=centromere_fs,
            scaffold_db_id=scaffold_db_id_resolved if use_scaffold_db else None,
        ),
    )
    tracker = progress.track([f"{mode}/{fs}" for mode in requested_modes for fs in requested])
    # The cascade's own reporters nest one level in, so their headlines read
    # as steps of this run rather than as separate commands that started on
    # their own.
    cascade_progress = progress.child()

    logger.info(
        "rendering karyotype(s): %d input(s), modes=%s, feature_sets=%s "
        "(= %d SVG(s) x %d format(s))%s",
        len(inputs),
        requested_modes,
        requested,
        len(requested_modes) * len(requested),
        len(requested_formats),
        " [combine-chromosomes]" if combine_chromosomes else "",
    )

    # Make sure scaffolded BEDs exist for every input + requested
    # feature set. scaffold_run short-circuits on existing files.
    if combine_chromosomes:
        # Combine path: cascade scaffold in 'both' mode so it writes the
        # combined-chromosome BEDs (and the combined FASTA + AGP). The
        # centromere-detection feature set is added to the set when
        # centromere mode is requested so its combined BED exists for
        # in-coordinate centromere detection below.
        scaffold_feature_sets = list(requested)
        if centromere_fs is not None and centromere_fs not in scaffold_feature_sets:
            scaffold_feature_sets.append(centromere_fs)
        scaffold_run(
            inputs,
            progress=cascade_progress,
            db_root=db_root,
            db_id=db_id_resolved,
            feature_sets=scaffold_feature_sets,
            mode="both",
            min_scaffold_length=min_scaffold_length,
            telo_motif=telo_motif,
            acrocentrics=acros_set,
            combine_chromosomes=True,
            scaffold_gap_size=scaffold_gap_size,
            combine_acrocentrics=combine_acrocentrics,
            split_haps_regex=split_haps_regex,
            threads=threads,
            bgzip=bgzip,
            auto=auto,
            output_dir=output_dir,
            annotation_variant=annotation_variant,
        )
    elif use_scaffold_db:
        # Two-database path: derive the scaffold map (chromosome ordering +
        # region orientation) from the LAYOUT database, then ensure the
        # plotted feature set's annotation BED exists against the PLOT
        # database. The map is applied to the plot BED at bin time in the
        # render loop -- the same code path --no-scaffolding uses.
        #
        # Layout pass: feature_sets=[] makes scaffold annotate only the
        # layout DB's role sets (chromosome + region) and write
        # {stem}.{scaffold_db_id}.scaffold_map.tsv; no per-feature-set
        # scaffolded BEDs are produced.
        scaffold_run(
            inputs,
            progress=cascade_progress,
            db_root=scaffold_db_root_eff,
            db_id=scaffold_db_id_resolved,
            feature_sets=[],
            mode="bed",
            min_scaffold_length=min_scaffold_length,
            telo_motif=telo_motif,
            acrocentrics=acrocentrics,
            split_haps_regex=split_haps_regex,
            threads=threads,
            bgzip=bgzip,
            auto=auto,
            output_dir=output_dir,
            write_scaffolded_beds=False,
            annotation_variant=annotation_variant,
        )
        # Plot pass: ensure each requested feature set's annotation BED
        # exists against the PLOT database. Skipped when already present;
        # honours --no-auto. Inputs are grouped by (output dir, missing
        # feature sets) so each group is annotated in a single batched call
        # — on the HKS backend that loads each feature set's index once for
        # the whole group (e.g. both haplotypes together) instead of once
        # per input.
        plot_groups: dict[tuple[Path, frozenset[str]], list[Path]] = {}
        for spec in inputs:
            stem = _input_stem(spec.path)
            out_dir = output_dir if output_dir is not None else spec.path.parent
            missing = [
                fs
                for fs in requested
                if not _annotation_bed_path(
                    out_dir, stem, db_id_resolved, fs, variant=annotation_variant
                ).is_file()
            ]
            if not missing:
                continue
            if not auto:
                raise KaryotypeError(
                    f"missing {annotation_variant} annotation BED(s) for "
                    f"{spec.path.name} feature set(s) {missing!r} in plot database "
                    f"{db_id_resolved!r}; re-run with auto-derive enabled or run "
                    f"`karyoscope annotate` first."
                )
            out_dir.mkdir(parents=True, exist_ok=True)
            plot_groups.setdefault((out_dir, frozenset(missing)), []).append(spec.path)

        for (out_dir, missing_fs), paths in plot_groups.items():
            annotate_batch(
                input_paths=paths,
                output_dir=out_dir,
                db_root=db_root,
                db_id=db_id_resolved,
                # Preserve the manifest/requested order for the batched call.
                feature_sets=[fs for fs in requested if fs in missing_fs],
                threads=threads,
                smooth=(annotation_variant == "smoothed"),
                bgzip=bgzip,
                progress=cascade_progress,
            )
    else:
        scaffold_run(
            inputs,
            progress=cascade_progress,
            db_root=db_root,
            db_id=db_id_resolved,
            feature_sets=requested,
            mode="bed",
            min_scaffold_length=min_scaffold_length,
            telo_motif=telo_motif,
            acrocentrics=acrocentrics,
            split_haps_regex=split_haps_regex,
            threads=threads,
            bgzip=bgzip,
            auto=auto,
            output_dir=output_dir,
            write_scaffolded_beds=scaffolding,
            annotation_variant=annotation_variant,
        )

    # For centromere mode, also ensure the centromere coordinates file
    # exists per input. centromeres_run cascades through scaffold and
    # bin internally. Only call when centromere mode is actually
    # requested -- it's expensive (an extra bin pass) and unnecessary
    # for genome / subtelomere outputs. On the combine path the
    # per-input combined centromere ranges are detected lazily in the
    # render loop instead (keyed by <chrom>_<hap>, in combined coords).
    if want_centromere and not combine_chromosomes:
        centromeres_run(
            inputs,
            db_root=scaffold_db_root_eff,
            db_id=scaffold_db_id_resolved,
            min_scaffold_length=min_scaffold_length,
            telo_motif=telo_motif,
            acrocentrics=acrocentrics,
            split_haps_regex=split_haps_regex,
            bgzip=bgzip,
            threads=threads,
            auto=auto,
            output_dir=output_dir,
            write_scaffolded_beds=scaffolding,
        )

    # Per (input, feature_set): bin scaffolded BED at mode-appropriate
    # size, load into memory.
    per_input_state: list[tuple[InputSpec, Path, str]] = []
    for spec in inputs:
        stem = _input_stem(spec.path)
        out_dir = output_dir if output_dir is not None else spec.path.parent
        out_dir.mkdir(parents=True, exist_ok=True)
        per_input_state.append((spec, out_dir, stem))

    # Sane default output base. Used only when output_path is None.
    if output_path is not None:
        # Strip .svg if present; we add it per-(mode, fs) below.
        base_name = output_path.name
        if base_name.lower().endswith(".svg"):
            base_name = base_name[:-4]
        results_dir = output_path.parent
    else:
        # Default the on-disk filename to a sample-identifying base so it
        # tells the user at a glance which sample produced the SVGs.
        # Multi-input runs (e.g. separate hap1/hap2 FASTAs) collapse to
        # the common prefix of the input stems (GM04890.haplotype1 +
        # GM04890.haplotype2 -> GM04890), matching what the title band
        # already shows and avoiding a misleading "...haplotype1..." name
        # for a both-hap plot. Pass --output to override.
        base_name = _common_base([stem for _, _, stem in per_input_state])
        results_dir = output_dir if output_dir is not None else per_input_state[0][1]

    # Sample label for the SVG title band. Defaults to the first
    # input's stem, joined with " + " for multi-input runs so a quick
    # glance at the title tells the reader which assembly produced it.
    if sample_label is None:
        stems = [stem for _, _, stem in per_input_state]
        sample_label = " + ".join(stems) if stems else None

    # Combine path: precompute the per-input combined map rows (one
    # synthetic <chrom>_<hap> row per combined object) and, when
    # centromere mode is requested, the combined centromere ranges.
    # Both are reused across every (mode, feature_set) render for that
    # input. Keyed by the input FASTA path.
    combined_rows_by_input: dict[Path, list[MapRow]] = {}
    combined_cen_by_input: dict[Path, dict[str, tuple[int, int]]] = {}
    if combine_chromosomes:
        for spec, out_dir, stem in per_input_state:
            per_contig_rows = read_map(out_dir / f"{stem}.{db_id_resolved}.scaffold_map.tsv")
            crows = combined_map_rows(
                per_contig_rows,
                acrocentrics=acros_set,
                combine_acrocentrics=combine_acrocentrics,
            )
            combined_rows_by_input[spec.path] = crows
            if want_centromere and centromere_fs is not None:
                combined_cen_by_input[spec.path] = _ensure_combined_centromeres(
                    out_dir=out_dir,
                    stem=stem,
                    db_id=db_id_resolved,
                    centromere_fs=centromere_fs,
                    centromere_leaves=centromere_leaves,
                    map_rows=crows,
                    auto=auto,
                    input_name=spec.path.name,
                    threads=threads,
                    bgzip=bgzip,
                    variant=annotation_variant,
                )

    results: list[KaryotypeResult] = []
    results_dir.mkdir(parents=True, exist_ok=True)

    # Longest scaffolded sequence across inputs -- drives the auto (data-
    # driven) genome-view bin size below (and, in the renderer, the zoom).
    max_seq_len = 0
    for _spec, out_dir, stem in per_input_state:
        map_path = out_dir / f"{stem}.{scaffold_db_id_resolved}.scaffold_map.tsv"
        if map_path.is_file():
            for row in read_map(map_path):
                max_seq_len = max(max_seq_len, row.length)

    # Outer loop: mode (sets bin_size). Inner loop: feature set.
    # Binned BEDs are cached per (input, fs, bin_size) on disk so the
    # second time we hit the same fs at a different bin_size, we just
    # do a fresh bin call -- not re-running scaffold/annotate.
    for current_mode in requested_modes:
        current_bin_size = (
            bin_size if bin_size is not None else _auto_bin_size(current_mode, max_seq_len)
        )

        for fs in requested:
            t_view = time.perf_counter()
            leaves = leaves_for(hierarchy, fs)
            # Legend sort: hand the renderer the list of child names
            # in hierarchy.tsv file order. Internal nodes appear first
            # (defined as children of "categorized"), then their
            # subtrees, mirroring the order the user wrote the file.
            fs_feature_order = [row.child for row in hierarchy.rows_in(fs)]

            render_inputs: list[RenderInput] = []
            for spec, out_dir, stem in per_input_state:
                # Read map first; the binner needs it on the
                # ``--no-scaffolding`` path to apply rename + flip at
                # bin time when no on-disk scaffolded BED exists. On the
                # combine path the renderer instead consumes the
                # synthetic combined map rows keyed by <chrom>_<hap>.
                per_contig_rows = read_map(
                    out_dir / f"{stem}.{scaffold_db_id_resolved}.scaffold_map.tsv"
                )
                render_map_rows = (
                    combined_rows_by_input[spec.path] if combine_chromosomes else per_contig_rows
                )
                binned_path = _ensure_binned_scaffolded(
                    out_dir=out_dir,
                    stem=stem,
                    db_id=db_id_resolved,
                    fs=fs,
                    bin_size=current_bin_size,
                    leaf_set=leaves,
                    auto=auto,
                    input_name=spec.path.name,
                    threads=threads,
                    map_rows=render_map_rows,
                    variant=annotation_variant,
                    combined=combine_chromosomes,
                )
                binned_bed = _load_binned_bed(binned_path)

                centromere_ranges: dict[str, tuple[int, int]] | None = None
                if current_mode == "centromere":
                    if combine_chromosomes:
                        centromere_ranges = combined_cen_by_input[spec.path]
                    else:
                        cpath = _centromeres_bed_path(out_dir, stem, scaffold_db_id_resolved)
                        if not cpath.is_file():
                            raise KaryotypeError(
                                f"missing centromeres BED for {spec.path.name} "
                                f"(expected at {cpath})"
                            )
                        centromere_ranges = _load_centromeres_bed(cpath)

                render_inputs.append(
                    RenderInput(
                        map_rows=render_map_rows,
                        binned_bed=binned_bed,
                        centromere_ranges=centromere_ranges,
                    )
                )

            flat_colors = colors_for_set(colors, fs)
            # Legend grouping is per feature set and optional; a database
            # without the column yields {} and the renderer keeps its
            # per-feature legend unchanged.
            fs_groups = legend_groups_all.get(fs, {})
            fs_group_order = list(dict.fromkeys(fs_groups.values())) if fs_groups else None
            combined_tag = f".{_COMBINED_TAG}" if combine_chromosomes else ""
            stem_for_paths = (
                f"{base_name}.{db_id_resolved}.{current_mode}.{fs}."
                f"{annotation_variant}{combined_tag}{colors_tag}.karyotype"
            )
            svg_path = results_dir / f"{stem_for_paths}.svg"
            logger.info(
                "rendering karyotype: mode=%s, feature_set=%s -> %s",
                current_mode,
                fs,
                svg_path.name,
            )
            t_render = time.perf_counter()
            render_karyotype(
                render_inputs,
                colors=flat_colors,
                legend_groups=fs_groups or None,
                legend_group_order=fs_group_order,
                mode=current_mode,  # type: ignore[arg-type]
                sex=sex,
                sex_determination_system=sex_determination_system,
                background_color=background_color,
                subtelomere_boundary=subtelomere_boundary,
                seed_human_chromosomes=seed_human_chromosomes,
                expected_chromosomes=expected_chromosomes,
                pixels_per_mb=pixels_per_mb,
                output_path=svg_path,
                sample_label=sample_label,
                database_id=db_id_resolved,
                feature_set_label=fs,
                smoothed=(annotation_variant == "smoothed"),
                show_title=show_title,
                show_legend=show_legend,
                feature_order=fs_feature_order,
            )
            logger.info("rendered %s in %.1fs", svg_path.name, time.perf_counter() - t_render)

            # Convert SVG to additional formats as requested. The SVG
            # itself stays on disk only if the user asked for it; for
            # PDF/PNG-only runs we delete the SVG after conversion.
            output_paths: list[Path] = []
            keep_svg = "svg" in requested_formats
            for fmt in requested_formats:
                if fmt == "svg":
                    output_paths.append(svg_path)
                    continue
                target = results_dir / f"{stem_for_paths}.{fmt}"
                logger.info("converting %s -> %s", svg_path.name, target.name)
                t_conv = time.perf_counter()
                convert_svg(svg_path, target)
                logger.info("converted %s in %.1fs", target.name, time.perf_counter() - t_conv)
                output_paths.append(target)
            if not keep_svg:
                svg_path.unlink(missing_ok=True)

            results.append(
                KaryotypeResult(
                    feature_set=fs,
                    mode=current_mode,
                    output_paths=output_paths,
                )
            )
            # Timed from the top of the inner loop, so the first view of a
            # mode carries its bin pass (and, on a cold run, the cascade
            # work its feature set needed) rather than reporting a
            # misleadingly fast render.
            tracker.step(f"{current_mode}/{fs}", time.perf_counter() - t_view)

    logger.info(
        "karyotype complete in %.1fs (%d render(s) -> %d file(s))",
        time.perf_counter() - t_karyo_start,
        len(results),
        sum(len(r.output_paths) for r in results),
    )
    return results


# Silence "unused" warnings on imports kept for future extension.
_ = (Iterable, DEFAULT_HUMAN_CHROMOSOMES, DEFAULT_HUMAN_ACROCENTRICS)
