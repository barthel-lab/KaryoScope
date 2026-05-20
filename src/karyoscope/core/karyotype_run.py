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
basename without ``.svg``) or, when no explicit path is given,
the literal string ``"karyotype"`` so the filename is sample-agnostic
(the sample is implied by ``<dbid>`` and the ``--outdir``).
"""

from __future__ import annotations

import gzip
import logging
import time
from collections import OrderedDict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from karyoscope.core.annotate import resolve_database
from karyoscope.core.bin import bin_features, leaves_for
from karyoscope.core.centromeres import centromeres_run
from karyoscope.core.io.colors import colors_for_set, parse_colors, validate_colors
from karyoscope.core.io.hierarchy import parse_hierarchy
from karyoscope.core.io.scaffold_map import read_map
from karyoscope.core.karyotype import (
    DEFAULT_HUMAN_CHROMOSOMES,
    RenderInput,
    convert_svg,
    render_karyotype,
)
from karyoscope.core.scaffold import (
    DEFAULT_HUMAN_ACROCENTRICS,
    DEFAULT_MIN_SCAFFOLD_LENGTH,
    Interval,
)
from karyoscope.core.scaffold_run import InputSpec, scaffold_run
from karyoscope.exceptions import KaryotypeError
from karyoscope.manifest import validate_database_layout

logger = logging.getLogger(__name__)


#: Default bin size (bp) per render mode. The archive's per-mode
#: pixel scales are tuned to these defaults; deviating substantially
#: produces under- or over-detailed SVGs.
DEFAULT_BIN_SIZE_BY_MODE: dict[str, int] = {
    "genome": 1_000_000,
    "centromere": 100_000,
    "subtelomere": 100,
}


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


def _scaffolded_bed_path(out_dir: Path, stem: str, db_id: str, fs: str) -> Path:
    gz = out_dir / f"{stem}.{db_id}.{fs}.smoothed.scaffolded.bed.gz"
    if gz.is_file():
        return gz
    plain = out_dir / f"{stem}.{db_id}.{fs}.smoothed.scaffolded.bed"
    if plain.is_file():
        return plain
    return gz


def _binned_scaffolded_bed_path(
    out_dir: Path,
    stem: str,
    db_id: str,
    fs: str,
    bin_size: int,
) -> Path:
    return out_dir / f"{stem}.{db_id}.{fs}.smoothed.scaffolded.binned{bin_size}.bed.gz"


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
) -> Path:
    out = _binned_scaffolded_bed_path(out_dir, stem, db_id, fs, bin_size)
    if out.is_file():
        return out
    if not auto:
        raise KaryotypeError(
            f"missing binned scaffolded BED for {input_name}, feature set "
            f"{fs!r}, bin size {bin_size} (expected at {out}). Re-run with "
            "auto-derive enabled."
        )
    src = _scaffolded_bed_path(out_dir, stem, db_id, fs)
    if not src.is_file():
        raise KaryotypeError(f"cannot bin {fs!r} for {input_name}: scaffolded BED missing at {src}")
    # bin_features logs its own start (with leaf_set + threads) and
    # completion lines; no need for a redundant announcement here.
    bin_features(src, out, bin_size=bin_size, leaf_set=leaf_set or None, threads=threads)
    return out


def _centromeres_bed_path(out_dir: Path, stem: str, db_id: str) -> Path:
    gz = out_dir / f"{stem}.{db_id}.centromeres.bed.gz"
    if gz.is_file():
        return gz
    plain = out_dir / f"{stem}.{db_id}.centromeres.bed"
    if plain.is_file():
        return plain
    return gz


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
    feature_sets: list[str] | None = None,
    modes: list[str] | None = None,
    sex: str | None = None,
    sex_determination_system: str = "XY",
    background_color: str = "white",
    bin_size: int | None = None,
    subtelomere_boundary: int = 250_000,
    min_scaffold_length: int = DEFAULT_MIN_SCAFFOLD_LENGTH,
    acrocentrics: set[str] | None = None,
    split_haps_regex: str | None = None,
    threads: int = 0,
    auto: bool = True,
    output_dir: Path | None = None,
    output_path: Path | None = None,
    seed_human_chromosomes: bool = True,
    formats: list[str] | None = None,
    sample_label: str | None = None,
    show_title: bool = True,
    show_legend: bool = True,
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

    colors = parse_colors(db_dir / manifest.colors)
    hierarchy = parse_hierarchy(db_dir / manifest.hierarchy)

    # Hard check that every hierarchy node has a colour. ``download``
    # runs the same check at install time; we re-run it here so
    # locally-installed databases (which bypass download validation)
    # still get caught before any SVG is rendered.
    color_issues = validate_colors(hierarchy, colors)
    if color_issues:
        raise KaryotypeError(
            f"database {db_id_resolved!r} is missing colour entries; refusing "
            "to render. Fix colors.tsv (or use `karyoscope info` to see the "
            "full list of issues):\n  - " + "\n  - ".join(color_issues)
        )

    logger.info(
        "rendering karyotype(s): %d input(s), modes=%s, feature_sets=%s "
        "(= %d SVG(s) x %d format(s))",
        len(inputs),
        requested_modes,
        requested,
        len(requested_modes) * len(requested),
        len(requested_formats),
    )

    # Make sure scaffolded BEDs exist for every input + requested
    # feature set. scaffold_run short-circuits on existing files.
    scaffold_run(
        inputs,
        db_root=db_root,
        db_id=db_id_resolved,
        feature_sets=requested,
        mode="bed",
        min_scaffold_length=min_scaffold_length,
        acrocentrics=acrocentrics,
        split_haps_regex=split_haps_regex,
        threads=threads,
        auto=auto,
        output_dir=output_dir,
    )

    # For centromere mode, also ensure the centromere coordinates file
    # exists per input. centromeres_run cascades through scaffold and
    # bin internally. Only call when centromere mode is actually
    # requested -- it's expensive (an extra bin pass) and unnecessary
    # for genome / subtelomere outputs.
    if "centromere" in requested_modes:
        centromeres_run(
            inputs,
            db_root=db_root,
            db_id=db_id_resolved,
            min_scaffold_length=min_scaffold_length,
            acrocentrics=acrocentrics,
            split_haps_regex=split_haps_regex,
            threads=threads,
            auto=auto,
            output_dir=output_dir,
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
        base_name = "karyotype"
        results_dir = output_dir if output_dir is not None else per_input_state[0][1]

    # Sample label for the SVG title band. Defaults to the first
    # input's stem, joined with " + " for multi-input runs so a quick
    # glance at the title tells the reader which assembly produced it.
    if sample_label is None:
        stems = [stem for _, _, stem in per_input_state]
        sample_label = " + ".join(stems) if stems else None

    results: list[KaryotypeResult] = []
    results_dir.mkdir(parents=True, exist_ok=True)

    # Outer loop: mode (sets bin_size). Inner loop: feature set.
    # Binned BEDs are cached per (input, fs, bin_size) on disk so the
    # second time we hit the same fs at a different bin_size, we just
    # do a fresh bin call -- not re-running scaffold/annotate.
    for current_mode in requested_modes:
        current_bin_size = (
            bin_size if bin_size is not None else DEFAULT_BIN_SIZE_BY_MODE[current_mode]
        )

        for fs in requested:
            leaves = leaves_for(hierarchy, fs)
            # Legend sort: hand the renderer the list of child names
            # in hierarchy.tsv file order. Internal nodes appear first
            # (defined as children of "categorized"), then their
            # subtrees, mirroring the order the user wrote the file.
            fs_feature_order = [row.child for row in hierarchy.rows_in(fs)]

            render_inputs: list[RenderInput] = []
            for spec, out_dir, stem in per_input_state:
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
                )
                binned_bed = _load_binned_bed(binned_path)
                map_rows = read_map(out_dir / f"{stem}.{db_id_resolved}.scaffold_map.tsv")

                centromere_ranges: dict[str, tuple[int, int]] | None = None
                if current_mode == "centromere":
                    cpath = _centromeres_bed_path(out_dir, stem, db_id_resolved)
                    if not cpath.is_file():
                        raise KaryotypeError(
                            f"missing centromeres BED for {spec.path.name} (expected at {cpath})"
                        )
                    centromere_ranges = _load_centromeres_bed(cpath)

                render_inputs.append(
                    RenderInput(
                        map_rows=map_rows,
                        binned_bed=binned_bed,
                        centromere_ranges=centromere_ranges,
                    )
                )

            flat_colors = colors_for_set(colors, fs)
            stem_for_paths = f"{base_name}.{db_id_resolved}.{current_mode}.{fs}.karyotype"
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
                mode=current_mode,  # type: ignore[arg-type]
                sex=sex,
                sex_determination_system=sex_determination_system,
                background_color=background_color,
                subtelomere_boundary=subtelomere_boundary,
                seed_human_chromosomes=seed_human_chromosomes,
                output_path=svg_path,
                sample_label=sample_label,
                database_id=db_id_resolved,
                feature_set_label=fs,
                smoothed=True,
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

    logger.info(
        "karyotype complete in %.1fs (%d render(s) -> %d file(s))",
        time.perf_counter() - t_karyo_start,
        len(results),
        sum(len(r.output_paths) for r in results),
    )
    return results


# Silence "unused" warnings on imports kept for future extension.
_ = (Iterable, DEFAULT_HUMAN_CHROMOSOMES, DEFAULT_HUMAN_ACROCENTRICS)
