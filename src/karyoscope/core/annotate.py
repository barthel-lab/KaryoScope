"""Orchestration for the ``annotate`` command.

The high-level flow:

1. Locate the database and validate its layout.
2. Invoke ``get_featureIDs`` (the vendored C++ helper) against the input
   FASTA. The result is a single combined BED whose 4th column holds
   integer feature ids.
3. Read the database's ``features.tsv`` to learn the mapping
   ``feature_id -> {feature_set: feature_name}``.
4. For each requested feature set, stream the combined BED, translate
   feature ids to feature names in that set, merge adjacent same-name
   intervals, and write the per-feature-set "presmoothed" BED.
5. Optionally bgzip the outputs (default on) and delete the combined
   intermediate (default on).

Smoothing — the hierarchy-aware "fill in less-specific labels using
flanking sequence" step — lives in :mod:`karyoscope.core.smooth`
(Stage 5c, not yet implemented). The BED written here is the
*presmoothed* track, identical to what smooth_features.py writes when
it processes the original (non-smoothed) variant.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from karyoscope import installed as _installed
from karyoscope.core.external import require_tool, run_tool
from karyoscope.core.io.features import Features, parse_features, render_feature
from karyoscope.core.io.kmc import run_get_featureids
from karyoscope.exceptions import (
    DatabaseNotFoundError,
    KaryoscopeError,
)
from karyoscope.manifest import validate_database_layout

logger = logging.getLogger(__name__)


#: FASTA extensions we recognise when deriving the output basename.
#: Order matters — longer extensions first so they win over shorter ones.
_FASTA_EXTENSIONS: tuple[str, ...] = (
    ".fasta.gz",
    ".fa.gz",
    ".fna.gz",
    ".fasta",
    ".fa",
    ".fna",
)


# --- result type ------------------------------------------------------


@dataclass
class AnnotateResult:
    """What :func:`annotate` produces.

    Attributes
    ----------
    output_paths
        ``{feature_set: path-to-output-BED}``. The path's extension is
        ``.bed.gz`` when bgzipped, ``.bed`` otherwise.
    combined_intermediate
        Path to the C++ combined BED, or ``None`` if it was deleted.
    """

    output_paths: dict[str, Path]
    combined_intermediate: Path | None


# --- helpers ----------------------------------------------------------


def _derive_input_basename(input_path: Path) -> str:
    """Strip recognised FASTA extensions from a path's filename.

    ``my_assembly.fa.gz`` -> ``my_assembly``. Falls back to the raw
    stem if no known extension matches.
    """
    name = input_path.name
    name_lower = name.lower()
    for ext in _FASTA_EXTENSIONS:
        if name_lower.endswith(ext):
            return name[: -len(ext)]
    return input_path.stem


def _resolve_database(
    db_root: Path,
    db_id: str | None,
) -> tuple[str, Path]:
    """Locate and validate the database to use.

    Returns ``(database_id, database_directory)``. Raises a clean error
    if the user requested an unknown id or if the layout is invalid.
    """
    state = _installed.load(db_root)

    if db_id is not None:
        record = state.databases.get(db_id)
        if record is None:
            installed_ids = sorted(state.databases.keys())
            raise DatabaseNotFoundError(
                f"database {db_id!r} is not installed at {db_root}. "
                f"Installed databases: {installed_ids or '(none)'}. "
                "Run `karyoscope download --list` to see what's available."
            )
        db_dir = db_root / record.directory
    else:
        # No explicit id — use the only installed db if there's exactly
        # one, otherwise complain.
        if not state.databases:
            raise DatabaseNotFoundError(
                f"no databases installed at {db_root}. "
                "Install one with `karyoscope download`, or pass --db <ID>."
            )
        if len(state.databases) > 1:
            ids = sorted(state.databases.keys())
            raise DatabaseNotFoundError(
                f"multiple databases installed; pass --db to choose one. Installed: {ids}"
            )
        db_id, record = next(iter(state.databases.items()))
        db_dir = db_root / record.directory

    if not db_dir.is_dir():
        raise DatabaseNotFoundError(
            f"database {db_id!r} is recorded but its directory {db_dir} is missing on disk."
        )

    # Surface layout errors early rather than letting get_featureIDs fail
    # with a less informative message later.
    validate_database_layout(db_dir)

    return db_id, db_dir


def _iter_bed_records(handle: TextIO) -> Iterator[tuple[str, int, int, int]]:
    """Yield ``(seq_name, start, end, feature_id)`` from a 4-column BED stream.

    Skips blank lines and lines starting with ``#``. Raises ``ValueError``
    on malformed numeric fields — the caller should reframe as a clean
    user error.
    """
    for line_no, raw in enumerate(handle, start=1):
        line = raw.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            raise ValueError(
                f"line {line_no}: expected at least 4 tab-separated columns, "
                f"got {len(parts)}: {raw!r}"
            )
        try:
            start = int(parts[1])
            end = int(parts[2])
            feature_id = int(parts[3])
        except ValueError as e:
            raise ValueError(
                f"line {line_no}: could not parse integers in BED record: {raw!r}"
            ) from e
        yield parts[0], start, end, feature_id


def _split_combined_bed(
    combined_bed: Path,
    feature_sets: list[str],
    features: Features,
    output_paths: dict[str, Path],
) -> None:
    """Read ``combined_bed`` and write one per-feature-set BED per requested set.

    Adjacent records on the same sequence whose translated feature name
    matches are merged on the fly (the C++ already does this for
    same-feature-id runs, but translation across feature sets can collapse
    distinct ids into the same name — e.g., two adjacent rows mapping to
    chr1 via different region ids).
    """
    # Per-feature-set running record. Stored as 4-tuples
    # (seq_name, start, end, current_feature_name) or None.
    pending: dict[str, tuple[str, int, int, str] | None] = {fs: None for fs in feature_sets}
    handles: dict[str, TextIO] = {fs: output_paths[fs].open("w") for fs in feature_sets}

    def flush(fs: str) -> None:
        rec = pending[fs]
        if rec is None:
            return
        seq, start, end, name = rec
        handles[fs].write(f"{seq}\t{start}\t{end}\t{name}\n")
        pending[fs] = None

    try:
        with combined_bed.open() as f:
            for seq_name, start, end, fid in _iter_bed_records(f):
                for fs in feature_sets:
                    name = render_feature(fid, fs, features)
                    rec = pending[fs]
                    if (
                        rec is not None
                        and rec[0] == seq_name
                        and rec[3] == name
                        and rec[2] == start
                    ):
                        # Extend the running interval.
                        pending[fs] = (rec[0], rec[1], end, rec[3])
                    else:
                        flush(fs)
                        pending[fs] = (seq_name, start, end, name)

        for fs in feature_sets:
            flush(fs)
    finally:
        for h in handles.values():
            h.close()


def _bgzip_file(path: Path) -> Path:
    """Compress ``path`` in-place with ``bgzip``, returning the new path.

    ``bgzip`` removes the source file by default (matches gzip's behaviour).
    Returns ``Path(str(path) + ".gz")``.
    """
    bgzip = require_tool(
        "bgzip",
        install_hint="Install htslib (`conda install -c bioconda htslib`), "
        "or rerun with --no-bgzip to skip compression.",
    )
    run_tool([bgzip, "-f", str(path)])
    return Path(str(path) + ".gz")


# --- main entry point -------------------------------------------------


def annotate(
    *,
    input_path: Path,
    output_dir: Path,
    db_root: Path,
    db_id: str | None = None,
    feature_sets: list[str] | None = None,
    threads: int = 0,
    keep_intermediates: bool = False,
    bgzip: bool = True,
) -> AnnotateResult:
    """Run the full annotate pipeline for one input FASTA.

    Parameters
    ----------
    input_path
        FASTA (plain or gzipped) to annotate.
    output_dir
        Directory to write per-feature-set BEDs to. Created if absent.
    db_root
        KaryoScope database root (typically ``ensure_db_root(...)``).
    db_id
        Specific database id to use, or ``None`` to pick the unique
        installed database (errors if there are zero or many).
    feature_sets
        Restrict output to these feature sets. ``None`` means "all sets
        declared in the database's manifest".
    threads
        Threads to pass to ``get_featureIDs``. ``0`` means auto.
    keep_intermediates
        Keep the combined ``...combined.presmoothed.featureIDs.bed`` from
        the C++ step. Default: delete after splitting.
    bgzip
        bgzip the per-feature-set output BEDs. Default: yes.

    Raises
    ------
    DatabaseNotFoundError
        If the database can't be resolved or its layout is broken.
    ToolNotFoundError
        If ``get_featureIDs`` (or ``bgzip``, when requested) isn't found.
    ExternalToolError
        If a subprocess exits non-zero.
    KaryoscopeError
        If a requested feature set isn't in the database's manifest, or
        the input file doesn't exist.
    """
    if not input_path.is_file():
        raise KaryoscopeError(f"input file not found: {input_path}")

    db_id_resolved, db_dir = _resolve_database(db_root, db_id)
    manifest = validate_database_layout(db_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Pick feature sets, validating the user's choices against the manifest.
    available = list(manifest.feature_sets)
    if feature_sets is None:
        requested = available
    else:
        unknown = [fs for fs in feature_sets if fs not in available]
        if unknown:
            raise KaryoscopeError(
                f"feature set(s) {unknown!r} not declared in {db_id_resolved}'s "
                f"manifest. Available: {available!r}"
            )
        requested = list(feature_sets)
    logger.info("annotating %s against %s, sets=%s", input_path, db_id_resolved, requested)

    # Parse features.tsv up front so we fail fast if it's malformed.
    features = parse_features(db_dir / manifest.features)
    # Sanity-check: the manifest's feature_sets should be a subset of
    # what's actually in features.tsv. Don't enforce equality — extra
    # columns in features.tsv aren't fatal.
    missing_from_table = [fs for fs in requested if fs not in features.feature_sets]
    if missing_from_table:
        raise KaryoscopeError(
            f"feature set(s) {missing_from_table!r} declared in manifest but "
            f"missing from features.tsv columns ({features.feature_sets!r})"
        )

    # Run the C++ helper. We pass an explicit prefix so its output name
    # uses our database id, not the KMC index basename ("features").
    input_basename = _derive_input_basename(input_path)
    prefix = f"{input_basename}.{db_id_resolved}"
    kmc_db_basename = db_dir / manifest.index.basename
    combined_bed = run_get_featureids(
        db_path=kmc_db_basename,
        input_path=input_path,
        output_dir=output_dir,
        threads=threads,
        prefix=prefix,
        capture=True,
    )
    if not combined_bed.is_file():
        raise KaryoscopeError(f"get_featureIDs did not produce expected output at {combined_bed}")
    logger.debug("combined BED at %s", combined_bed)

    # Compute per-feature-set output paths (uncompressed names; bgzip
    # later if requested).
    output_paths: dict[str, Path] = {
        fs: output_dir / f"{prefix}.{fs}.presmoothed.bed" for fs in requested
    }
    logger.debug("writing %d per-feature-set BED(s)", len(output_paths))
    _split_combined_bed(combined_bed, requested, features, output_paths)

    # bgzip (or not).
    if bgzip:
        for fs in requested:
            output_paths[fs] = _bgzip_file(output_paths[fs])

    # Tidy up the combined intermediate unless asked to keep it.
    if not keep_intermediates:
        try:
            combined_bed.unlink()
            combined_kept: Path | None = None
            logger.debug("removed combined intermediate %s", combined_bed)
        except OSError as e:
            logger.warning("could not remove intermediate %s: %s", combined_bed, e)
            combined_kept = combined_bed
    else:
        combined_kept = combined_bed

    return AnnotateResult(
        output_paths=output_paths,
        combined_intermediate=combined_kept,
    )
