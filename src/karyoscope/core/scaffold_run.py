"""High-level orchestration for ``karyoscope scaffold``.

This is the auto-derive cascade and per-input fan-out that the CLI
layer sits on top of. It is intentionally separate from the algorithm
in :mod:`karyoscope.core.scaffold` so the algorithm stays testable in
isolation without file system or subprocess concerns.

Flow for one ``karyoscope scaffold`` invocation:

1. Resolve the database (id + dir, validated layout).
2. Read manifest roles to pick the chromosome-assignment and
   region-assignment feature sets, falling back to literal
   ``"chromosome"`` / ``"region"`` with a warning when the manifest
   omits them.
3. Per ``-i NAME=PATH`` input, in order, assign a hap label
   (explicit > filename-stem inference > positional ``inputN``).
4. Per input, derive any missing artefacts: annotation BEDs (via
   :func:`karyoscope.core.annotate.annotate`), telomere file (via
   ``seqtk telo``), 1 Mb binned BEDs (via
   :func:`karyoscope.core.bin.bin_features`). Existing files are
   reused. ``--no-auto`` (passed as ``auto=False``) disables this
   and turns missing inputs into hard errors.
5. Build :class:`ContigInput` records by reading the per-contig hap
   classifications (:func:`karyoscope.core.hap_inference.classify_contigs`),
   loading the binned BEDs, and parsing the telomere file.
6. Run :func:`karyoscope.core.scaffold.classify_and_orient` once
   across all inputs to get the global ordering.
7. Per input, write ``<stem>.<dbid>.scaffold_map.tsv``,
   ``<stem>.<dbid>.scaffold_stats.tsv``, and the per-feature-set
   ``<stem>.<dbid>.<fs>.smoothed.scaffolded.bed[.gz]`` files.

The orchestrator runs all of step 4 sequentially per input. The
design committed to this for Stage 5d-1 with the structural caveat
that the future move to per-input parallelism should be mechanical;
that's why the cascade is expressed as a list of independent
``(input, tool)`` units rather than a deeply nested function.
"""

from __future__ import annotations

import logging
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from karyoscope.core.annotate import _bgzip_file, annotate, resolve_database
from karyoscope.core.bin import bin_features, leaves_for
from karyoscope.core.hap_inference import (
    assign_per_input_labels,
    classify_contigs,
    read_fasta_contig_names,
)
from karyoscope.core.io.hierarchy import parse_hierarchy
from karyoscope.core.io.scaffold_map import MapRow, write_legacy_stats, write_map
from karyoscope.core.io.telo import TeloFlags, parse_telo_file, run_seqtk_telo
from karyoscope.core.scaffold import (
    DEFAULT_HUMAN_ACROCENTRICS,
    DEFAULT_MIN_SCAFFOLD_LENGTH,
    ContigInput,
    classify_and_orient,
    rewrite_bed,
    rewrite_fasta,
)
from karyoscope.exceptions import ScaffoldError
from karyoscope.manifest import validate_database_layout

logger = logging.getLogger(__name__)


#: Recognised FASTA suffixes (longest first). Used to derive the
#: input "stem" for output naming.
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


@dataclass
class InputSpec:
    """One ``-i`` argument plus its optional ``--telo`` companion."""

    name: str | None  # explicit hap label, or None for auto-inference
    path: Path
    telo_path: Path | None = None


@dataclass
class _ResolvedInput:
    """Post-name-resolution: stem, hap label, and (optional) telo override."""

    spec: InputSpec
    hap_label: str
    stem: str
    out_dir: Path


#: Output modes for :func:`scaffold_run`. ``"fasta"`` writes only the
#: scaffolded FASTA per input; ``"bed"`` writes only the per-feature-set
#: scaffolded BEDs (used internally by ``karyoscope karyotype``);
#: ``"both"`` writes both. The map and legacy stats files are always
#: written.
ScaffoldMode = str  # Literal["fasta", "bed", "both"]; loose for back-compat
_VALID_MODES: tuple[str, ...] = ("fasta", "bed", "both")


@dataclass
class ScaffoldResult:
    """What :func:`scaffold_run` wrote for one input."""

    input_path: Path
    hap_label: str
    map_path: Path
    stats_path: Path
    scaffolded_beds: dict[str, Path] = field(default_factory=dict)
    scaffolded_fasta: Path | None = None


# --- helpers --------------------------------------------------------


def _resolve_roles(manifest_roles: dict[str, str], available: list[str]) -> tuple[str, str]:
    """Pick the chromosome- and region-assignment feature sets.

    Falls back to literal names with a warning when the manifest omits
    them. Errors when the resolved set is not declared in
    ``available``.
    """
    chrom_set = manifest_roles.get("chromosome_assignment")
    if chrom_set is None:
        chrom_set = "chromosome"
        logger.warning(
            "manifest has no roles.chromosome_assignment; falling back to "
            "feature set %r. Add `roles: { chromosome_assignment: <name> }` "
            "to manifest.yaml to silence this warning.",
            chrom_set,
        )
    if chrom_set not in available:
        raise ScaffoldError(
            f"chromosome-assignment feature set {chrom_set!r} not declared in the "
            f"database's manifest. Available feature sets: {available!r}"
        )

    region_set = manifest_roles.get("region_assignment")
    if region_set is None:
        region_set = "region"
        logger.warning(
            "manifest has no roles.region_assignment; falling back to "
            "feature set %r. Add `roles: { region_assignment: <name> }` "
            "to manifest.yaml to silence this warning.",
            region_set,
        )
    if region_set not in available:
        raise ScaffoldError(
            f"region-assignment feature set {region_set!r} not declared in the "
            f"database's manifest. Available feature sets: {available!r}"
        )
    return chrom_set, region_set


def _smoothed_bed_path(out_dir: Path, stem: str, db_id: str, fs: str) -> Path:
    """The canonical path for ``annotate``'s smoothed BED output."""
    # annotate writes .bed.gz by default; we look for either.
    gz = out_dir / f"{stem}.{db_id}.{fs}.smoothed.bed.gz"
    if gz.is_file():
        return gz
    plain = out_dir / f"{stem}.{db_id}.{fs}.smoothed.bed"
    if plain.is_file():
        return plain
    return gz  # default to expecting .gz


def _binned_bed_path(out_dir: Path, stem: str, db_id: str, fs: str, bin_size: int) -> Path:
    return out_dir / f"{stem}.{db_id}.{fs}.smoothed.binned{bin_size}.bed.gz"


def _telo_path(out_dir: Path, stem: str) -> Path:
    return out_dir / f"{stem}.telo"


def _ensure_annotated(
    spec: _ResolvedInput,
    *,
    db_root: Path,
    db_id: str,
    feature_sets: list[str],
    threads: int,
    auto: bool,
    bgzip: bool,
) -> None:
    """Run annotate if any required smoothed BED is missing.

    ``bgzip`` controls whether the auto-derived ``annotate`` outputs
    are compressed; it should inherit the scaffold-level setting so
    a user running ``karyoscope scaffold --no-bgzip`` on a fresh input
    gets uncompressed annotation BEDs from the cascade too (matching
    what a manual ``karyoscope annotate --no-bgzip`` would have
    produced). Downstream readers handle both ``.bed`` and ``.bed.gz``
    transparently via ``chunked_seq_reader``, so the choice only
    affects on-disk storage, not behaviour.
    """
    missing = [
        fs
        for fs in feature_sets
        if not _smoothed_bed_path(spec.out_dir, spec.stem, db_id, fs).is_file()
    ]
    if not missing:
        return
    if not auto:
        raise ScaffoldError(
            f"missing annotation BEDs for {spec.spec.path.name} "
            f"(feature sets: {missing}). Re-run with auto-derive enabled, or run "
            f"`karyoscope annotate` first."
        )
    logger.info(
        "annotating %s for missing feature set(s) %s",
        spec.spec.path,
        missing,
    )
    annotate(
        input_path=spec.spec.path,
        output_dir=spec.out_dir,
        db_root=db_root,
        db_id=db_id,
        feature_sets=missing,
        threads=threads,
        smooth=True,
        keep_presmoothed=True,
        keep_intermediates=False,
        bgzip=bgzip,
    )


def _ensure_telo(spec: _ResolvedInput, *, auto: bool) -> Path:
    """Return the telo file path, running seqtk if necessary."""
    if spec.spec.telo_path is not None:
        if not spec.spec.telo_path.is_file():
            raise ScaffoldError(f"--telo file does not exist: {spec.spec.telo_path}")
        return spec.spec.telo_path

    out = _telo_path(spec.out_dir, spec.stem)
    if out.is_file():
        return out
    if not auto:
        raise ScaffoldError(
            f"missing telomere file for {spec.spec.path.name} (expected at "
            f"{out}). Re-run with auto-derive enabled, pass --telo explicitly, "
            f"or run `seqtk telo` first."
        )
    run_seqtk_telo(spec.spec.path, out)
    return out


def _ensure_binned(
    spec: _ResolvedInput,
    *,
    db_id: str,
    fs: str,
    bin_size: int,
    leaf_set: set[str],
    auto: bool,
    threads: int,
) -> Path:
    """Return the binned-BED path, running the binner if necessary."""
    out = _binned_bed_path(spec.out_dir, spec.stem, db_id, fs, bin_size)
    if out.is_file():
        return out
    if not auto:
        raise ScaffoldError(
            f"missing binned BED for {spec.spec.path.name}, feature set {fs!r}, "
            f"bin size {bin_size} (expected at {out}). Re-run with auto-derive "
            f"enabled."
        )
    src = _smoothed_bed_path(spec.out_dir, spec.stem, db_id, fs)
    if not src.is_file():
        raise ScaffoldError(
            f"can't bin {fs!r} for {spec.spec.path.name}: smoothed BED missing "
            f"at {src} (annotate should have produced it)"
        )
    # bin_features logs its own start (with leaf_set + threads) and
    # completion lines; no need for a redundant announcement here.
    bin_features(src, out, bin_size=bin_size, leaf_set=leaf_set or None, threads=threads)
    return out


def _load_binned_bed(path: Path) -> dict[str, list[tuple[int, int, str]]]:
    import gzip

    opener = gzip.open if path.suffix == ".gz" else open
    out: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
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
            out[parts[0]].append((start, end, parts[3]))
    return out


# --- main entry point -----------------------------------------------


def scaffold_run(
    inputs: list[InputSpec],
    *,
    db_root: Path,
    db_id: str | None = None,
    feature_sets: list[str] | None = None,
    mode: ScaffoldMode = "fasta",
    bin_size: int = 1_000_000,
    min_scaffold_length: int = DEFAULT_MIN_SCAFFOLD_LENGTH,
    acrocentrics: set[str] | None = None,
    split_haps_regex: str | None = None,
    threads: int = 0,
    bgzip: bool = True,
    keep_unscaffolded: bool = True,
    auto: bool = True,
    output_dir: Path | None = None,
) -> dict[str, ScaffoldResult]:
    """Run the full ``karyoscope scaffold`` pipeline.

    Parameters
    ----------
    inputs
        One :class:`InputSpec` per ``-i`` argument. Order matters
        only for the positional ``inputN`` fallback labels.
    db_root, db_id
        Database location. ``db_id=None`` picks the unique installed
        database (errors when zero or many are installed).
    feature_sets
        Which per-feature-set scaffolded BEDs to emit. ``None``
        defaults to every feature set in the manifest.
    bin_size
        Bin size in bp for the orientation BEDs. Matches the
        manuscript's 1 Mb benchmark by default.
    min_scaffold_length
        Drop contigs shorter than this that have no telomere.
    acrocentrics
        Chromosome names that count as acrocentric in the flip
        decision. ``None`` falls back to the human default with a
        single warning (handled at the CLI layer).
    split_haps_regex
        Optional user regex applied per-contig (overrides the
        built-in hap inference patterns). See
        :mod:`karyoscope.core.hap_inference`.
    threads
        Threads for any auto-run ``annotate`` invocations.
    bgzip
        bgzip the scaffolded output BEDs (default: yes).
    auto
        Auto-derive missing inputs (annotate, seqtk telo, bin).
        Disable to require every input to exist.
    output_dir
        Where the scaffolded BEDs and the map file land. ``None``
        means "next to each input FASTA" (the directory the
        scaffolded outputs for that input share with the existing
        annotation BEDs).
    """
    if not inputs:
        raise ScaffoldError("at least one --input is required")
    if mode not in _VALID_MODES:
        raise ScaffoldError(f"unknown mode {mode!r}; expected one of {_VALID_MODES}")

    t_scaffold_start = time.perf_counter()
    db_id_resolved, db_dir = resolve_database(db_root, db_id)
    manifest = validate_database_layout(db_dir)
    available = list(manifest.feature_sets)

    chromosome_fs, region_fs = _resolve_roles(manifest.roles, available)
    if mode == "fasta":
        # FASTA-only mode never writes per-feature-set scaffolded BEDs,
        # so there's no point requesting them from annotate. The role
        # sets (used for classify+orient) are still annotated.
        requested: list[str] = []
    elif feature_sets is None:
        requested = list(available)
    else:
        unknown = [fs for fs in feature_sets if fs not in available]
        if unknown:
            raise ScaffoldError(
                f"requested feature set(s) {unknown!r} not declared in manifest "
                f"(available: {available!r})"
            )
        requested = list(feature_sets)

    # We always need the role sets present in annotate output, even when
    # the user didn't ask for them in --feature-set (and even in
    # mode='fasta' -- the role sets drive classify_and_orient).
    # Preserve the user's / manifest's order so progress messages and
    # output files appear in the order users expect; only append role
    # sets at the end if they weren't already requested. (A previous
    # version used ``sorted(set(...))`` which alphabetised, surfacing
    # in dogfooding as "acrocentric" arriving first instead of the
    # manifest's "chromosome".)
    seen = set(requested)
    annotate_sets = list(requested) + [fs for fs in (chromosome_fs, region_fs) if fs not in seen]

    hierarchy = parse_hierarchy(db_dir / manifest.hierarchy)
    chromosome_leaves = leaves_for(hierarchy, chromosome_fs)
    region_leaves = leaves_for(hierarchy, region_fs)
    if not chromosome_leaves:
        raise ScaffoldError(
            f"chromosome feature set {chromosome_fs!r} has no leaf features; "
            "cannot classify contigs to chromosomes"
        )

    # Per-input resolution: hap labels and output directories.
    name_path_pairs = [(spec.name, spec.path) for spec in inputs]
    hap_labels = assign_per_input_labels(name_path_pairs)
    resolved: list[_ResolvedInput] = []
    for spec, hap_label in zip(inputs, hap_labels, strict=True):
        stem = _input_stem(spec.path)
        out_dir = output_dir if output_dir is not None else spec.path.parent
        out_dir.mkdir(parents=True, exist_ok=True)
        resolved.append(_ResolvedInput(spec=spec, hap_label=hap_label, stem=stem, out_dir=out_dir))

    is_only_input = len(resolved) == 1
    acros_set = acrocentrics if acrocentrics is not None else set(DEFAULT_HUMAN_ACROCENTRICS)

    logger.info(
        "scaffolding %d input(s) against %s (mode=%s)",
        len(resolved),
        db_id_resolved,
        mode,
    )

    # --- auto-derive cascade per input -------------------------------
    # Sequential for Stage 5d-1; the per-input units are independent
    # so a future batched-parallel pass would be a structural change
    # in this loop only.
    for r in resolved:
        if not r.spec.path.is_file():
            raise ScaffoldError(f"input file not found: {r.spec.path}")
        _ensure_annotated(
            r,
            db_root=db_root,
            db_id=db_id_resolved,
            feature_sets=annotate_sets,
            threads=threads,
            auto=auto,
            bgzip=bgzip,
        )
        _ensure_telo(r, auto=auto)
        _ensure_binned(
            r,
            db_id=db_id_resolved,
            fs=chromosome_fs,
            bin_size=bin_size,
            leaf_set=chromosome_leaves,
            auto=auto,
            threads=threads,
        )
        _ensure_binned(
            r,
            db_id=db_id_resolved,
            fs=region_fs,
            bin_size=bin_size,
            leaf_set=region_leaves,
            auto=auto,
            threads=threads,
        )

    # --- build ContigInput list across all inputs --------------------
    all_contigs: list[ContigInput] = []
    contigs_per_input: dict[str, list[ContigInput]] = defaultdict(list)
    for r in resolved:
        contig_names = read_fasta_contig_names(r.spec.path)
        explicit = r.spec.name is not None
        per_contig_hap = classify_contigs(
            contig_names,
            file_level_label=r.hap_label,
            split_haps_regex=split_haps_regex,
            is_only_input=is_only_input,
            explicit_name_given=explicit,
        )
        chrom_bins = _load_binned_bed(
            _binned_bed_path(r.out_dir, r.stem, db_id_resolved, chromosome_fs, bin_size)
        )
        region_bins = _load_binned_bed(
            _binned_bed_path(r.out_dir, r.stem, db_id_resolved, region_fs, bin_size)
        )
        telo_path = _ensure_telo(r, auto=auto)
        telo_flags = parse_telo_file(telo_path)

        for name in contig_names:
            length = max(
                (stop for _, stop, _ in chrom_bins.get(name, [])),
                default=0,
            )
            region_max = max(
                (stop for _, stop, _ in region_bins.get(name, [])),
                default=0,
            )
            length = max(length, region_max)
            if length == 0:
                # No coverage at all -- annotate produced no records for
                # this contig in either role feature set. Skip silently.
                continue
            ci = ContigInput(
                input_name=per_contig_hap[name],
                input_file=r.spec.path.name,
                contig_name=name,
                length=length,
                chromosome_bins=sorted(chrom_bins.get(name, [])),
                region_bins=sorted(region_bins.get(name, [])),
                telo=telo_flags.get(name, TeloFlags(False, False)),
            )
            all_contigs.append(ci)
            contigs_per_input[r.spec.path.name].append(ci)

    # --- classify + orient (joint across all inputs) -----------------
    logger.info(
        "classifying + orienting %d contig(s) across %d input(s)",
        len(all_contigs),
        len(resolved),
    )
    t_classify_start = time.perf_counter()
    rows = classify_and_orient(
        all_contigs,
        chromosome_leaves=chromosome_leaves,
        acrocentrics=acros_set,
        min_scaffold_length=min_scaffold_length,
    )
    logger.info(
        "classified %d scaffold row(s) in %.1fs",
        len(rows),
        time.perf_counter() - t_classify_start,
    )

    # Group rows by input file for the writer.
    rows_per_input: dict[str, list[MapRow]] = defaultdict(list)
    for row in rows:
        rows_per_input[row.input_file].append(row)

    # --- write outputs per input -------------------------------------
    results: dict[str, ScaffoldResult] = {}
    for r in resolved:
        per_input_rows = rows_per_input.get(r.spec.path.name, [])
        map_path = r.out_dir / f"{r.stem}.{db_id_resolved}.scaffold_map.tsv"
        stats_path = r.out_dir / f"{r.stem}.{db_id_resolved}.scaffold_stats.tsv"
        write_map(per_input_rows, map_path)
        write_legacy_stats(per_input_rows, stats_path)

        scaffolded_beds: dict[str, Path] = {}
        if mode in ("bed", "both"):
            for fs in requested:
                src = _smoothed_bed_path(r.out_dir, r.stem, db_id_resolved, fs)
                if not src.is_file():
                    logger.warning(
                        "smoothed BED for %s / %s not found at %s; skipping",
                        r.spec.path.name,
                        fs,
                        src,
                    )
                    continue
                out_plain = r.out_dir / f"{r.stem}.{db_id_resolved}.{fs}.smoothed.scaffolded.bed"
                logger.info(
                    "rewriting scaffolded BED for %s / %s -> %s",
                    r.spec.path.name,
                    fs,
                    out_plain.name,
                )
                t_rb = time.perf_counter()
                rewrite_bed(src, out_plain, map_rows=per_input_rows, gzip_out=False)
                logger.info("wrote %s in %.1fs", out_plain.name, time.perf_counter() - t_rb)
                out_final = _bgzip_file(out_plain, threads=threads) if bgzip else out_plain
                scaffolded_beds[fs] = out_final

        scaffolded_fasta: Path | None = None
        if mode in ("fasta", "both"):
            # Write plain .fa first, then bgzip (when requested) so the
            # output is samtools-faidx compatible -- matching how the
            # BED outputs are compressed above.
            fasta_plain = r.out_dir / f"{r.stem}.{db_id_resolved}.scaffolded.fa"
            logger.info(
                "writing scaffolded FASTA for %s -> %s",
                r.spec.path.name,
                fasta_plain,
            )
            t_rf = time.perf_counter()
            rewrite_fasta(
                r.spec.path,
                fasta_plain,
                map_rows=per_input_rows,
                keep_unscaffolded=keep_unscaffolded,
                gzip_out=False,
            )
            logger.info("wrote %s in %.1fs", fasta_plain.name, time.perf_counter() - t_rf)
            scaffolded_fasta = _bgzip_file(fasta_plain, threads=threads) if bgzip else fasta_plain

        results[r.spec.path.name] = ScaffoldResult(
            input_path=r.spec.path,
            hap_label=r.hap_label,
            map_path=map_path,
            stats_path=stats_path,
            scaffolded_beds=scaffolded_beds,
            scaffolded_fasta=scaffolded_fasta,
        )

    logger.info(
        "scaffold complete in %.1fs (%d input(s))",
        time.perf_counter() - t_scaffold_start,
        len(resolved),
    )
    return results
