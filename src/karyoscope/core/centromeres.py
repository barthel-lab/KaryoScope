"""Centromere coordinate extraction from binned scaffolded region BEDs.

Faithful port of the archive's ``get_centromeres.py``. Per scaffolded
contig, identify the centromere's start/end coordinates by:

1. **Coarse pass** (typically 1 Mb bins): for each contig, find the
   minimum start and maximum stop of bins whose feature classifies
   as ``"centromere"`` via
   :func:`karyoscope.core.scaffold.get_simple_region`. This gives a
   conservative outer bound for the centromere region.
2. **Fine pass** (typically 100 kb bins, optional): search within
   ``[coarse_min - 1Mb, coarse_max + 1Mb]`` (a buffer either side of
   the coarse range) for centromere-classified bins, and tighten the
   coordinates. Bins outside the buffer window are ignored so that
   stray centromeric signal far from the actual centromere doesn't
   widen the call.

Output: ``{contig_name: (start, end)}`` for contigs that have any
centromere-classified bins in the coarse pass. Contigs with no
centromere are omitted entirely.

The input bins must come from the *scaffolded* coordinate system
(post-flip), so the output coordinates are directly consumable by
``karyoscope karyotype --mode centromere``.
"""

from __future__ import annotations

import gzip
import logging
import shutil
import tempfile
import time
from collections import OrderedDict
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from karyoscope.core.annotate import _bgzip_file, resolve_database
from karyoscope.core.bin import bin_features, leaves_for
from karyoscope.core.io.hierarchy import parse_hierarchy
from karyoscope.core.io.scaffold_map import MapRow, read_map
from karyoscope.core.scaffold import (
    DEFAULT_MIN_SCAFFOLD_LENGTH,
    Interval,
    get_simple_region,
    rewrite_bed,
)
from karyoscope.core.scaffold_run import InputSpec, scaffold_run
from karyoscope.exceptions import CentromereError
from karyoscope.manifest import validate_database_layout

logger = logging.getLogger(__name__)


#: Default search-window buffer for the fine-refinement pass. Matches
#: the archive's hardcoded constant.
DEFAULT_REFINEMENT_BUFFER = 1_000_000


def find_centromere_ranges(
    coarse_bins: Mapping[str, list[Interval]],
    fine_bins: Mapping[str, list[Interval]] | None = None,
    *,
    refinement_buffer: int = DEFAULT_REFINEMENT_BUFFER,
) -> OrderedDict[str, tuple[int, int]]:
    """Return ``{contig: (start, end)}`` for every contig with centromeric content.

    Parameters
    ----------
    coarse_bins
        Per-contig list of intervals from the coarse binned scaffolded
        region BED (typically 1 Mb bins).
    fine_bins
        Optional per-contig list from a finer binned scaffolded region
        BED (typically 100 kb). When provided, the coarse range is
        tightened by considering only bins within
        ``[coarse_min - refinement_buffer, coarse_max + refinement_buffer]``.
    refinement_buffer
        How far either side of the coarse range to look during the
        fine pass. Defaults to 1 Mb. Set to 0 to require the fine
        bins to be strictly inside the coarse range.

    Insertion order of the returned dict follows ``coarse_bins``.
    Contigs absent from ``coarse_bins`` (or present but with no
    centromere-classified bin) are absent from the result.
    """
    coarse_ranges: OrderedDict[str, tuple[int, int]] = OrderedDict()
    for contig, bins in coarse_bins.items():
        cmin = -1
        cmax = 0
        for start, stop, feature in bins:
            if get_simple_region(feature) != "centromere":
                continue
            if cmin == -1 or start < cmin:
                cmin = start
            if stop > cmax:
                cmax = stop
        if cmax > 0:
            coarse_ranges[contig] = (max(cmin, 0), cmax)

    if fine_bins is None:
        return coarse_ranges

    out: OrderedDict[str, tuple[int, int]] = OrderedDict()
    for contig, coarse_range in coarse_ranges.items():
        window_start = max(0, coarse_range[0] - refinement_buffer)
        window_end = coarse_range[1] + refinement_buffer
        fmin = -1
        fmax = 0
        for start, stop, feature in fine_bins.get(contig, []):
            if start < window_start or stop > window_end:
                continue
            if get_simple_region(feature) != "centromere":
                continue
            if fmin == -1 or start < fmin:
                fmin = start
            if stop > fmax:
                fmax = stop
        if fmax > 0:
            out[contig] = (max(fmin, 0), fmax)
        else:
            # No fine signal -- keep the coarse range so the contig
            # still gets a centromere call.
            out[contig] = coarse_range

    return out


# --- orchestrator ---------------------------------------------------


#: Defaults for the orchestrator (overridable via the CLI flags).
DEFAULT_COARSE_BIN_SIZE = 1_000_000
DEFAULT_FINE_BIN_SIZE = 100_000


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


def _resolve_centromere_role(manifest_roles: dict[str, str], available: list[str]) -> str:
    """Pick the feature set to use for centromere detection.

    Fall-back chain: ``roles.centromere_detection`` ->
    ``roles.region_assignment`` -> the literal name ``"region"``. Each
    fallback step emits a warning. The final resolved set must be
    declared in the manifest's ``feature_sets``; otherwise we raise.
    """
    fs = manifest_roles.get("centromere_detection")
    if fs is None:
        fs = manifest_roles.get("region_assignment")
        if fs is None:
            fs = "region"
            logger.warning(
                "manifest has neither roles.centromere_detection nor "
                "roles.region_assignment; falling back to feature set %r. "
                "Add `roles: { centromere_detection: <name> }` to "
                "manifest.yaml to silence this warning.",
                fs,
            )
        else:
            logger.warning(
                "manifest has no roles.centromere_detection; falling back "
                "to roles.region_assignment value %r.",
                fs,
            )
    if fs not in available:
        raise CentromereError(
            f"centromere-detection feature set {fs!r} not declared in the "
            f"database's manifest. Available feature sets: {available!r}"
        )
    return fs


@dataclass
class CentromereResult:
    """What :func:`centromeres_run` produced for one input."""

    input_path: Path
    hap_label: str
    centromeres_bed: Path
    num_contigs: int


def _scaffolded_bed_path(out_dir: Path, stem: str, db_id: str, fs: str) -> Path:
    """The conventional path scaffold writes for a per-feature-set BED."""
    gz = out_dir / f"{stem}.{db_id}.{fs}.smoothed.scaffolded.bed.gz"
    if gz.is_file():
        return gz
    plain = out_dir / f"{stem}.{db_id}.{fs}.smoothed.scaffolded.bed"
    if plain.is_file():
        return plain
    return gz  # default to expecting .gz


def _smoothed_bed_path(out_dir: Path, stem: str, db_id: str, fs: str) -> Path:
    """The annotate-produced (unscaffolded) smoothed BED path.

    Looked at only on the ``write_scaffolded_beds=False`` codepath
    (when scaffolding skipped writing per-FS scaffolded BEDs). The
    binning step then runs against this file and the scaffold map is
    applied post-bin via :func:`rewrite_bed`.
    """
    gz = out_dir / f"{stem}.{db_id}.{fs}.smoothed.bed.gz"
    if gz.is_file():
        return gz
    plain = out_dir / f"{stem}.{db_id}.{fs}.smoothed.bed"
    if plain.is_file():
        return plain
    return gz  # default to expecting .gz


def _binned_scaffolded_bed_path(
    out_dir: Path,
    stem: str,
    db_id: str,
    fs: str,
    bin_size: int,
) -> Path:
    return out_dir / f"{stem}.{db_id}.{fs}.smoothed.scaffolded.binned{bin_size}.bed.gz"


def _load_binned_bed(path: Path) -> OrderedDict[str, list[Interval]]:
    """Load a binned BED into ``{contig: [(start, stop, feature), ...]}``.

    Preserves insertion order so the downstream output BED is emitted
    in the same per-contig order the scaffolded BED used (which is
    the canonical chromosome x hap x category ordering from
    :func:`classify_and_orient`).
    """
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
            contig = parts[0]
            out.setdefault(contig, []).append((start, end, parts[3]))
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
) -> Path:
    """Return the binned scaffolded BED path, building it if missing.

    Two construction paths:

    * **Scaffolded BED on disk** (the historical path): bin it directly.
    * **Scaffolded BED missing but a map is available** (the
      ``--no-scaffolding`` path): bin the *smoothed* (unscaffolded)
      BED at the requested ``bin_size``, then stream the binned output
      through :func:`rewrite_bed` to apply the map (rename contigs +
      mirror coordinates for flipped contigs). Skipping the
      full-resolution scaffold BED rewrite is the whole point of
      ``--no-scaffolding`` -- the map application on binned data is
      microseconds rather than the minutes the smoothed rewrite costs.
    """
    out = _binned_scaffolded_bed_path(out_dir, stem, db_id, fs, bin_size)
    if out.is_file():
        return out
    if not auto:
        raise CentromereError(
            f"missing binned scaffolded BED for {input_name}, feature set "
            f"{fs!r}, bin size {bin_size} (expected at {out}). Re-run with "
            f"auto-derive enabled."
        )
    scaffolded_src = _scaffolded_bed_path(out_dir, stem, db_id, fs)
    if scaffolded_src.is_file():
        # bin_features logs its own start (with leaf_set + threads) and
        # completion lines; no need for a redundant announcement here.
        bin_features(
            scaffolded_src, out, bin_size=bin_size, leaf_set=leaf_set or None, threads=threads
        )
        return out

    # Fallback: bin the smoothed BED, then apply the scaffold map.
    if map_rows is None:
        raise CentromereError(
            f"cannot bin {fs!r} for {input_name}: scaffolded BED missing at "
            f"{scaffolded_src} and no scaffold map provided for post-bin "
            f"renaming."
        )
    smoothed_src = _smoothed_bed_path(out_dir, stem, db_id, fs)
    if not smoothed_src.is_file():
        raise CentromereError(
            f"cannot bin {fs!r} for {input_name}: smoothed BED missing at "
            f"{smoothed_src} (and scaffolded BED also missing)"
        )
    tmpdir = Path(tempfile.mkdtemp(prefix="ks_centro_bin_", dir=out_dir))
    try:
        tmp_binned = tmpdir / "binned.bed.gz"
        bin_features(
            smoothed_src,
            tmp_binned,
            bin_size=bin_size,
            leaf_set=leaf_set or None,
            threads=threads,
        )
        # rewrite_bed applies the map: rename per row.new_name, mirror
        # coords for flipped contigs. Output respects the .gz suffix.
        rewrite_bed(tmp_binned, out, map_rows=map_rows)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return out


def centromeres_run(
    inputs: list[InputSpec],
    *,
    db_root: Path,
    db_id: str | None = None,
    centromere_feature_set: str | None = None,
    coarse_bin_size: int = DEFAULT_COARSE_BIN_SIZE,
    fine_bin_size: int | None = DEFAULT_FINE_BIN_SIZE,
    min_scaffold_length: int = DEFAULT_MIN_SCAFFOLD_LENGTH,
    acrocentrics: set[str] | None = None,
    split_haps_regex: str | None = None,
    threads: int = 0,
    bgzip: bool = True,
    auto: bool = True,
    output_dir: Path | None = None,
    write_scaffolded_beds: bool = True,
) -> dict[str, CentromereResult]:
    """Per input, extract per-contig centromere coordinates.

    Cascades through scaffold (which itself cascades through annotate,
    seqtk telo, and bin) so the user can run this against raw FASTAs.

    Parameters
    ----------
    inputs
        Same :class:`InputSpec` form as :func:`scaffold_run` -- one per
        FASTA, optional explicit hap name and optional precomputed
        telomere file.
    centromere_feature_set
        Override the manifest's ``roles.centromere_detection``. ``None``
        (the default) falls back through the manifest role chain.
    coarse_bin_size, fine_bin_size
        Bin sizes for the coarse pass and (optional) fine refinement
        pass. Pass ``fine_bin_size=None`` or ``0`` to disable the
        fine pass entirely.
    auto
        Auto-derive any missing prerequisites (scaffolded BED, binned
        scaffolded BEDs). Disable to require they exist already.
    output_dir
        Where to write the centromere BEDs. ``None`` means "next to
        each input FASTA" -- same convention as scaffold.

    Returns ``{input_path_basename: CentromereResult}`` per input.
    """
    if not inputs:
        raise CentromereError("at least one --input is required")
    if coarse_bin_size < 1:
        raise CentromereError(
            f"--coarse-bin-size must be a positive integer, got {coarse_bin_size}"
        )
    if fine_bin_size is not None and fine_bin_size < 0:
        raise CentromereError(f"--fine-bin-size must be >= 0 (0 disables), got {fine_bin_size}")

    t_cen_start = time.perf_counter()
    db_id_resolved, db_dir = resolve_database(db_root, db_id)
    manifest = validate_database_layout(db_dir)
    available = list(manifest.feature_sets)

    if centromere_feature_set is None:
        centromere_fs = _resolve_centromere_role(manifest.roles, available)
    elif centromere_feature_set not in available:
        raise CentromereError(
            f"--centromere-feature-set {centromere_feature_set!r} not declared "
            f"in the database's manifest. Available: {available!r}"
        )
    else:
        centromere_fs = centromere_feature_set

    hierarchy = parse_hierarchy(db_dir / manifest.hierarchy)
    centromere_leaves = leaves_for(hierarchy, centromere_fs)

    logger.info(
        "extracting centromeres from %d input(s) against %s (feature_set=%s)",
        len(inputs),
        db_id_resolved,
        centromere_fs,
    )

    # Make sure the scaffolded region BEDs exist for every input by
    # calling scaffold_run with mode='bed'. scaffold_run itself is
    # cheap when its outputs already exist (the per-step _ensure_*
    # helpers short-circuit on existing files).
    scaffold_run(
        inputs,
        db_root=db_root,
        db_id=db_id_resolved,
        feature_sets=[centromere_fs],
        mode="bed",
        min_scaffold_length=min_scaffold_length,
        acrocentrics=acrocentrics,
        split_haps_regex=split_haps_regex,
        threads=threads,
        bgzip=bgzip,
        auto=auto,
        output_dir=output_dir,
        write_scaffolded_beds=write_scaffolded_beds,
    )

    # Per input, bin the scaffolded BED at the two bin sizes (skipping
    # the fine pass when disabled), load, run find_centromere_ranges,
    # and write the output BED.
    results: dict[str, CentromereResult] = {}
    for spec in inputs:
        stem = _input_stem(spec.path)
        out_dir = output_dir if output_dir is not None else spec.path.parent
        out_dir.mkdir(parents=True, exist_ok=True)

        # When the cascade was told to skip writing scaffolded BEDs
        # (write_scaffolded_beds=False), the binning helper needs the
        # scaffold map so it can apply the rename + flip at bin time
        # rather than reading from a non-existent scaffolded BED.
        map_path = out_dir / f"{stem}.{db_id_resolved}.scaffold_map.tsv"
        map_rows = read_map(map_path) if map_path.is_file() else None

        coarse_path = _ensure_binned_scaffolded(
            out_dir=out_dir,
            stem=stem,
            db_id=db_id_resolved,
            fs=centromere_fs,
            bin_size=coarse_bin_size,
            leaf_set=centromere_leaves,
            auto=auto,
            input_name=spec.path.name,
            threads=threads,
            map_rows=map_rows,
        )
        coarse_bins = _load_binned_bed(coarse_path)

        fine_bins: OrderedDict[str, list[Interval]] | None = None
        if fine_bin_size:
            fine_path = _ensure_binned_scaffolded(
                out_dir=out_dir,
                stem=stem,
                db_id=db_id_resolved,
                fs=centromere_fs,
                bin_size=fine_bin_size,
                leaf_set=centromere_leaves,
                auto=auto,
                input_name=spec.path.name,
                threads=threads,
                map_rows=map_rows,
            )
            fine_bins = _load_binned_bed(fine_path)

        logger.info(
            "finding centromere ranges for %s (contigs=%d)",
            spec.path.name,
            len(coarse_bins),
        )
        t_find = time.perf_counter()
        ranges = find_centromere_ranges(coarse_bins, fine_bins)
        logger.info(
            "found %d centromere range(s) for %s in %.1fs",
            len(ranges),
            spec.path.name,
            time.perf_counter() - t_find,
        )

        # Resolve hap label by reading the scaffold map (which scaffold
        # wrote during the cascade above) so it shows up in the result
        # summary. We don't strictly need it -- the centromere BED uses
        # the encoded contig name -- but it's nice for the CLI output.
        hap_label = _read_hap_label(out_dir, stem, db_id_resolved)

        plain = out_dir / f"{stem}.{db_id_resolved}.centromeres.bed"
        with plain.open("w") as h:
            for contig, (cstart, cend) in ranges.items():
                h.write(f"{contig}\t{cstart}\t{cend}\n")
        final = _bgzip_file(plain, threads=threads) if bgzip else plain

        results[spec.path.name] = CentromereResult(
            input_path=spec.path,
            hap_label=hap_label,
            centromeres_bed=final,
            num_contigs=len(ranges),
        )

    logger.info(
        "centromeres complete in %.1fs (%d input(s))",
        time.perf_counter() - t_cen_start,
        len(inputs),
    )
    return results


def _read_hap_label(out_dir: Path, stem: str, db_id: str) -> str:
    """Best-effort: read the hap label from the scaffold map file.

    Falls back to an empty string if the map file is unreadable for
    any reason. Only used for the CLI summary; not load-bearing.
    """
    map_path = out_dir / f"{stem}.{db_id}.scaffold_map.tsv"
    if not map_path.is_file():
        return ""
    try:
        with map_path.open("r") as h:
            header = h.readline()
            del header
            first = h.readline().rstrip("\n").split("\t")
            if len(first) >= 4:
                return first[3]
    except OSError:
        pass
    return ""
