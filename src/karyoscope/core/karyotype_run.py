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
from collections import OrderedDict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from karyoscope.core.annotate import resolve_database
from karyoscope.core.bin import bin_features, leaves_for
from karyoscope.core.centromeres import centromeres_run
from karyoscope.core.io.colors import colors_for_set, parse_colors
from karyoscope.core.io.hierarchy import parse_hierarchy
from karyoscope.core.io.scaffold_map import read_map
from karyoscope.core.karyotype import (
    DEFAULT_HUMAN_CHROMOSOMES,
    RenderInput,
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
    "full": 1_000_000,
    "centromere": 100_000,
    "subtelomere": 100,
}


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
    logger.info("binning %s -> %s (bin_size=%d)", src, out, bin_size)
    bin_features(src, out, bin_size=bin_size, leaf_set=leaf_set or None)
    return out


def _centromeres_bed_path(out_dir: Path, stem: str, db_id: str) -> Path:
    gz = out_dir / f"{stem}.{db_id}.centromeres.bed.gz"
    if gz.is_file():
        return gz
    plain = out_dir / f"{stem}.{db_id}.centromeres.bed"
    if plain.is_file():
        return plain
    return gz


@dataclass
class KaryotypeResult:
    """One SVG written by :func:`karyotype_run`."""

    feature_set: str
    mode: str
    svg_path: Path


def karyotype_run(
    inputs: list[InputSpec],
    *,
    db_root: Path,
    db_id: str | None = None,
    feature_sets: list[str] | None = None,
    mode: str = "full",
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
) -> list[KaryotypeResult]:
    """Render one SVG per requested feature set.

    Most parameters mirror the CLI flags; see
    :mod:`karyoscope.commands.karyotype` for the user-facing
    documentation.
    """
    if not inputs:
        raise KaryotypeError("at least one --input is required")
    if mode not in ("full", "subtelomere", "centromere"):
        raise KaryotypeError(
            f"unknown mode {mode!r}; expected 'full', 'subtelomere', or 'centromere'"
        )

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

    if bin_size is None:
        bin_size = DEFAULT_BIN_SIZE_BY_MODE[mode]
    if bin_size < 1:
        raise KaryotypeError(f"--bin-size must be a positive integer, got {bin_size}")

    colors = parse_colors(db_dir / manifest.colors)
    hierarchy = parse_hierarchy(db_dir / manifest.hierarchy)

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
    # bin internally.
    if mode == "centromere":
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

    results: list[KaryotypeResult] = []
    for fs in requested:
        leaves = leaves_for(hierarchy, fs)

        render_inputs: list[RenderInput] = []
        for spec, out_dir, stem in per_input_state:
            binned_path = _ensure_binned_scaffolded(
                out_dir=out_dir,
                stem=stem,
                db_id=db_id_resolved,
                fs=fs,
                bin_size=bin_size,
                leaf_set=leaves,
                auto=auto,
                input_name=spec.path.name,
            )
            binned_bed = _load_binned_bed(binned_path)
            map_rows = read_map(out_dir / f"{stem}.{db_id_resolved}.scaffold_map.tsv")

            centromere_ranges: dict[str, tuple[int, int]] | None = None
            if mode == "centromere":
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
        svg_path = results_dir / f"{base_name}.{db_id_resolved}.{mode}.{fs}.karyotype.svg"
        results_dir.mkdir(parents=True, exist_ok=True)
        render_karyotype(
            render_inputs,
            colors=flat_colors,
            mode=mode,  # type: ignore[arg-type]
            sex=sex,
            sex_determination_system=sex_determination_system,
            background_color=background_color,
            subtelomere_boundary=subtelomere_boundary,
            seed_human_chromosomes=seed_human_chromosomes,
            output_path=svg_path,
        )
        results.append(KaryotypeResult(feature_set=fs, mode=mode, svg_path=svg_path))

    return results


# Silence "unused" warnings on imports kept for future extension.
_ = (Iterable, DEFAULT_HUMAN_CHROMOSOMES, DEFAULT_HUMAN_ACROCENTRICS)
