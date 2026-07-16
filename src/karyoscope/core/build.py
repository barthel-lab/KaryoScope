"""Orchestration for ``karyoscope build``: turn a genome + per-feature-set BED
annotations into a complete, registry-ready HKS database.

Pipeline (see :mod:`karyoscope.core.buildspec` for the input model):

1. Per feature set, turn the BED into per-leaf FASTAs — optionally flattening
   overlaps and always gap-filling with a named background leaf — then derive
   the hierarchy edges, colours, and (priority mode) a priorities file.
2. Build the shared base index once with ``hks build-base``.
3. Layer each feature set on with ``hks add-feature-set``: priority mode when a
   priority file was given, an opt-in variable-k mode, else a plain fixed-k
   labeling (the default, which is what ``annotate`` queries at k=s).
4. Emit ``hierarchy.tsv`` / ``colors.tsv`` / ``manifest.yaml``, validate them,
   and register the database so the data commands can find it.

The heavy per-feature FASTAs live under ``<db_dir>/_build_work`` and are removed
on success unless ``keep_intermediates`` is set.
"""

from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from karyoscope import installed as _installed
from karyoscope._version import __version__
from karyoscope.core.buildspec import BuildSpec, FeatureSetSpec
from karyoscope.core.external import require_tool, run_tool
from karyoscope.core.io import emit, partition
from karyoscope.core.io.hierarchy import (
    REQUIRED_ROOT,
    Hierarchy,
    HierarchyRow,
    validate_hierarchy,
)
from karyoscope.core.io.hks import (
    get_hks_binary,
    run_hks_add_feature_set,
    run_hks_build_base,
    validate_sibling_priorities,
)
from karyoscope.exceptions import BuildError
from karyoscope.installed import InstalledRecord, now_iso
from karyoscope.manifest import validate_database_layout

logger = logging.getLogger(__name__)


@dataclass
class FeatureSetResult:
    """What was produced for one feature set."""

    name: str
    mode: str  # "priority", "fixed", or "variable-k"
    leaves: list[str]
    background: str | None


@dataclass
class BuildResult:
    """Outcome of :func:`build_database`."""

    db_id: str
    db_dir: Path
    feature_sets: list[FeatureSetResult] = field(default_factory=list)
    registered: bool = False


# --- small input parsers ---------------------------------------------


def _parse_priority_file(path: Path) -> tuple[list[tuple[str, str]] | None, dict[str, int]]:
    """Parse a priority file.

    Accepts two whitespace-separated forms:

    * 3-column ``child priority parent`` (the archive format) — yields both the
      hierarchy edges and the child priorities.
    * 2-column ``name priority`` — priorities only (no edges).

    Returns ``(edges_or_None, {name: priority})``.
    """
    edges: list[tuple[str, str]] = []
    prios: dict[str, int] = {}
    saw_edges = False
    with path.open() as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 3:
                child, prio_s, parent = parts
                saw_edges = True
                edges.append((child, parent))
                prios[child] = _int_priority(prio_s, path, lineno)
            elif len(parts) == 2:
                name, prio_s = parts
                prios[name] = _int_priority(prio_s, path, lineno)
            else:
                raise BuildError(
                    f"{path}:{lineno}: priority line must have 2 (name priority) or "
                    f"3 (child priority parent) columns: {line!r}"
                )
    return (edges if saw_edges else None), prios


def _int_priority(value: str, path: Path, lineno: int) -> int:
    try:
        return int(value)
    except ValueError as e:
        raise BuildError(f"{path}:{lineno}: priority must be an integer, got {value!r}") from e


def _parse_edge_list(path: Path) -> list[tuple[str, str]]:
    """Parse a header-less ``child parent`` (whitespace-separated) edge list."""
    edges: list[tuple[str, str]] = []
    with path.open() as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise BuildError(f"{path}:{lineno}: hierarchy edge needs 'child parent': {line!r}")
            edges.append((parts[0], parts[1]))
    return edges


def _parse_set_colors(path: Path) -> dict[str, str]:
    """Parse a per-set colours file into ``{feature: hex}``.

    Accepts ``feature <tab> color`` or ``feature_set <tab> feature <tab> color``;
    a header row (last column literally ``color``) is skipped.
    """
    out: dict[str, str] = {}
    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if parts[-1].strip().lower() == "color":
                continue  # header
            if len(parts) >= 3:
                feature, color = parts[1].strip(), parts[2].strip()
            elif len(parts) == 2:
                feature, color = parts[0].strip(), parts[1].strip()
            else:
                continue
            out[feature] = color
    return out


# --- per-feature-set preparation -------------------------------------


@dataclass
class _PreparedSet:
    name: str
    fof_path: Path  # feature-file-list
    names_path: Path  # feature-names
    hierarchy_txt: Path  # index/features.<fs>.hierarchy.txt
    priorities_path: Path | None  # workdir priorities file (priority mode) or None
    edges: list[tuple[str, str]]  # child -> parent
    leaves: list[str]
    node_colors: dict[str, str]
    feature_fastas: list[Path]  # for base-index input in mode B


def _prepare_feature_set(
    fs: FeatureSetSpec,
    *,
    spec: BuildSpec,
    index_dir: Path,
    work_dir: Path,
) -> _PreparedSet:
    set_work = work_dir / fs.name
    set_work.mkdir(parents=True, exist_ok=True)

    # -- priorities / hierarchy inputs --
    prio_edges: list[tuple[str, str]] | None = None
    prios: dict[str, int] = {}
    if fs.priority is not None:
        prio_edges, prios = _parse_priority_file(fs.priority)

    explicit_edges: list[tuple[str, str]] | None = None
    if fs.hierarchy is not None:
        explicit_edges = _parse_edge_list(fs.hierarchy)
    elif prio_edges is not None:
        explicit_edges = prio_edges

    # -- obtain (label -> FASTA) inputs --
    if fs.mode == "bed":
        intervals = partition.parse_bed(fs.bed)
        if fs.flatten:
            order = (
                [n for n, _ in sorted(prios.items(), key=lambda kv: kv[1])]
                if prios
                else partition.labels_in(intervals)
            )
            intervals = partition.flatten_by_priority(intervals, order)
        slice_intervals = list(intervals)
        if fs.background is not None:
            fai = partition.read_fai(_ensure_fai(spec.sequence))
            slice_intervals += partition.compute_background_intervals(intervals, fai, fs.background)
        label_paths = partition.slice_features_to_fastas(
            spec.sequence, slice_intervals, k=spec.s, outdir=set_work / "fastas"
        )
        leaves = list(label_paths.keys())
        feature_fastas = list(label_paths.values())
    else:
        # Mode B: FASTAs are already per-feature.
        label_paths, feature_fastas = _mode_b_inputs(fs, set_work)
        leaves = list(label_paths.keys())

    if not leaves:
        raise BuildError(
            f"feature set {fs.name!r}: produced no feature sequences (check that the "
            "BED contig names match the genome, or that the FASTA inputs are non-empty)"
        )

    # -- hierarchy edges --
    edges = _resolve_edges(fs, leaves, explicit_edges)

    # Validate the tree shape via the shared validator.
    hierarchy = Hierarchy(rows=[HierarchyRow(fs.name, c, p) for c, p in edges])
    issues = validate_hierarchy(hierarchy)
    if issues:
        raise BuildError(
            f"feature set {fs.name!r}: hierarchy is not a well-formed tree:\n  "
            + "\n  ".join(issues)
        )

    nodes = hierarchy.nodes(fs.name)

    # -- colours --
    provided = _parse_set_colors(fs.colors) if fs.colors is not None else {}
    node_colors = emit.assign_colors(
        nodes=nodes, leaves=leaves, background=fs.background, provided=provided
    )

    # -- write per-set index artifacts --
    hierarchy_txt = index_dir / f"features.{fs.name}.hierarchy.txt"
    emit.write_feature_set_hierarchy_txt(hierarchy_txt, edges)

    fof_path = set_work / "feature_files.txt"
    names_path = set_work / "feature_names.txt"
    fof_path.write_text("".join(f"{p}\n" for p in feature_fastas))
    names_path.write_text("".join(f"{leaf}\n" for leaf in leaves))

    # -- priorities (priority mode only) --
    priorities_path: Path | None = None
    if fs.priority is not None:
        full_prio = {node: prios.get(node, 0) for node in nodes}
        parent_of = {c: p for c, p in edges}
        sib_issues = validate_sibling_priorities(parent_of, full_prio)
        if sib_issues:
            raise BuildError(
                f"feature set {fs.name!r}: invalid priorities:\n  " + "\n  ".join(sib_issues)
            )
        priorities_path = set_work / "priorities.txt"
        emit.write_priorities_file(priorities_path, full_prio)

    return _PreparedSet(
        name=fs.name,
        fof_path=fof_path,
        names_path=names_path,
        hierarchy_txt=hierarchy_txt,
        priorities_path=priorities_path,
        edges=edges,
        leaves=leaves,
        node_colors=node_colors,
        feature_fastas=feature_fastas,
    )


def _mode_b_inputs(fs: FeatureSetSpec, set_work: Path) -> tuple[dict[str, Path], list[Path]]:
    """Resolve mode-B (per-feature FASTA) inputs to ``({label: path}, [paths])``.

    For ``fastas`` the label is each file's stem. ``per_seq_file`` is passed
    through as a single file whose sequence *names* are the labels (HKS reads
    them directly), so there is no per-label path map to build.
    """
    if fs.per_seq_file is not None:
        # One file, one sequence per feature; labels come from sequence names.
        return {fs.per_seq_file.stem: fs.per_seq_file}, [fs.per_seq_file]
    label_paths: dict[str, Path] = {}
    for p in fs.fastas or []:
        label = p.name
        for suffix in (".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna"):
            if label.lower().endswith(suffix):
                label = label[: -len(suffix)]
                break
        if label in label_paths:
            raise BuildError(
                f"feature set {fs.name!r}: two input FASTAs map to the same feature "
                f"name {label!r}; give distinct filenames"
            )
        label_paths[label] = p
    return label_paths, list(label_paths.values())


def _resolve_edges(
    fs: FeatureSetSpec,
    leaves: list[str],
    explicit_edges: list[tuple[str, str]] | None,
) -> list[tuple[str, str]]:
    """Return the ``child -> parent`` edge list for one set.

    With no explicit hierarchy, produce a flat star: every leaf (including the
    background) is a child of ``categorized``. With an explicit hierarchy, every
    produced leaf must appear in it; the background leaf is auto-attached to the
    root if the author left it out.
    """
    if explicit_edges is None:
        return [(leaf, REQUIRED_ROOT) for leaf in leaves]

    edges = list(explicit_edges)
    known = {c for c, _ in edges} | {p for _, p in edges}
    # Auto-attach the background leaf if omitted.
    if fs.background is not None and fs.background in leaves and fs.background not in known:
        edges.append((fs.background, REQUIRED_ROOT))
        known.add(fs.background)

    missing = [leaf for leaf in leaves if leaf not in known]
    if missing:
        raise BuildError(
            f"feature set {fs.name!r}: {len(missing)} produced label(s) are absent from "
            f"the provided hierarchy: {missing if len(missing) <= 8 else [*missing[:8], '...']}"
        )
    return edges


def _ensure_fai(fasta: Path) -> Path:
    """Return the ``.fai`` for ``fasta``, creating it with ``samtools faidx`` if absent."""
    fai = Path(str(fasta) + ".fai")
    if fai.is_file():
        return fai
    samtools = require_tool(
        "samtools",
        install_hint="Install samtools (conda install -c bioconda samtools) or run "
        "`samtools faidx <genome>` yourself to create the .fai index used for gap-fill.",
    )
    logger.info("indexing %s with samtools faidx", fasta.name)
    run_tool([samtools, "faidx", str(fasta)])
    if not fai.is_file():
        raise BuildError(f"samtools faidx did not produce {fai}")
    return fai


# --- top-level orchestration -----------------------------------------


def build_database(
    spec: BuildSpec,
    db_root: Path,
    *,
    register: bool = True,
    force: bool = False,
    keep_intermediates: bool = False,
) -> BuildResult:
    """Build (and optionally register) an HKS database from ``spec``.

    Parameters
    ----------
    spec
        The validated build request.
    db_root
        Database root; the database is written to ``db_root/<spec.id>``.
    register
        Record the database in ``installed.json`` so the data commands find it.
    force
        Overwrite an existing database directory and/or install record.
    keep_intermediates
        Keep the per-feature FASTA working directory instead of deleting it.
    """
    spec.validate()
    get_hks_binary()  # fail fast if the builder isn't installed

    db_root = db_root.resolve()
    db_dir = (db_root / spec.id).resolve()
    if db_dir.exists():
        if not force:
            raise BuildError(
                f"database directory already exists: {db_dir}. Pass --force to overwrite."
            )
        shutil.rmtree(db_dir)
    index_dir = db_dir / "index"
    index_dir.mkdir(parents=True, exist_ok=True)
    work_dir = db_dir / "_build_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    basename = "index/features"
    base_path = db_dir / f"{basename}.hksb"

    try:
        # 1. Prepare every feature set (writes hierarchy.txt + FASTAs into work).
        prepared: list[_PreparedSet] = [
            _prepare_feature_set(fs, spec=spec, index_dir=index_dir, work_dir=work_dir)
            for fs in spec.feature_sets
        ]

        # 2. Base index (once). Inputs: the genome (mode A) plus any mode-B FASTAs.
        _build_base_index(spec, prepared, base_path, work_dir)

        # 3. Layer each feature set.
        results: list[FeatureSetResult] = []
        for fs, prep in zip(spec.feature_sets, prepared, strict=True):
            hksf = db_dir / f"{basename}.{fs.name}.hksf"
            use_priority = prep.priorities_path is not None
            if use_priority:
                mode = "priority"
            elif fs.variable_k:
                mode = "variable-k"
            else:
                mode = "fixed"
            logger.info(
                "adding feature set %r (%s mode, %d leaves)", fs.name, mode, len(prep.leaves)
            )
            run_hks_add_feature_set(
                base_path=base_path,
                output_path=hksf,
                feature_set_name=fs.name,
                feature_names=prep.names_path,
                feature_hierarchy=prep.hierarchy_txt,
                feature_file_list=None if fs.per_seq_file is not None else prep.fof_path,
                feature_per_seq_file=fs.per_seq_file if fs.per_seq_file is not None else None,
                feature_priorities=prep.priorities_path,
                variable_k_support=fs.variable_k,
                forward_only=spec.forward_only,
                threads=spec.threads,
            )
            results.append(
                FeatureSetResult(
                    name=fs.name,
                    mode=mode,
                    leaves=prep.leaves,
                    background=fs.background,
                )
            )

        # 4. Aggregate artifacts + manifest.
        _emit_aggregate_artifacts(spec, prepared, db_dir)

        # 5. Validate the finished layout, then register.
        manifest = validate_database_layout(db_dir)
        registered = False
        if register:
            registered = _register_built(db_root, db_dir, manifest.version, force=force)

        result = BuildResult(
            db_id=spec.id, db_dir=db_dir, feature_sets=results, registered=registered
        )
    finally:
        if not keep_intermediates:
            shutil.rmtree(work_dir, ignore_errors=True)

    return result


def _build_base_index(
    spec: BuildSpec,
    prepared: list[_PreparedSet],
    base_path: Path,
    work_dir: Path,
) -> None:
    # HKS variable-k needs a dummy node for each feature sequence's start, which
    # only exists when the base is built from the feature sequences themselves.
    # So when *any* set requests variable_k, build the base from every set's
    # per-feature FASTAs (mode-A genome slices included) instead of the genome —
    # the feature k-mer set is identical (each slice is k-1-extended, capturing
    # boundary k-mers), but now every feature run starts at a sequence start.
    # Otherwise keep the simpler whole-genome base for mode-A sets.
    any_variable_k = any(fs.variable_k for fs in spec.feature_sets)

    if any_variable_k:
        base_inputs: list[Path] = []
        seen: set[Path] = set()
        for prep in prepared:
            for p in prep.feature_fastas:
                if p not in seen:
                    seen.add(p)
                    base_inputs.append(p)
    else:
        mode_b_fastas: list[Path] = []
        for fs, prep in zip(spec.feature_sets, prepared, strict=True):
            if fs.mode != "bed":
                mode_b_fastas.extend(prep.feature_fastas)
        base_inputs = []
        if spec.sequence is not None:
            base_inputs.append(spec.sequence)
        base_inputs.extend(mode_b_fastas)

    logger.info(
        "building HKS base index (s=%d, %d input file(s), variable_k=%s)",
        spec.s,
        len(base_inputs),
        any_variable_k,
    )
    if len(base_inputs) == 1 and not any_variable_k and spec.sequence is not None:
        run_hks_build_base(
            output_path=base_path,
            s=spec.s,
            input_path=base_inputs[0],
            threads=spec.threads,
            mem_gigas=spec.mem_gigas,
            external_memory=spec.external_memory,
            forward_only=spec.forward_only,
        )
    else:
        fof = work_dir / "base_inputs.txt"
        fof.write_text("".join(f"{p}\n" for p in base_inputs))
        run_hks_build_base(
            output_path=base_path,
            s=spec.s,
            input_file_list=fof,
            threads=spec.threads,
            mem_gigas=spec.mem_gigas,
            external_memory=spec.external_memory,
            forward_only=spec.forward_only,
        )


def _emit_aggregate_artifacts(spec: BuildSpec, prepared: list[_PreparedSet], db_dir: Path) -> None:
    hierarchy_rows: list[tuple[str, str, str]] = []
    color_rows: list[tuple[str, str, str]] = []
    for prep in prepared:
        for child, parent in prep.edges:
            hierarchy_rows.append((prep.name, child, parent))
        for node, color in prep.node_colors.items():
            color_rows.append((prep.name, node, color))

    emit.write_hierarchy_tsv(db_dir / "hierarchy.tsv", hierarchy_rows)
    emit.write_colors_tsv(db_dir / "colors.tsv", color_rows)

    # Advertise variable-k only when every feature set supports it: a fixed-mode
    # .hksf can only be queried at k=s, so a mixed database is fixed from a
    # querying standpoint.
    kmer_type = (
        "variable"
        if spec.feature_sets and all(fs.variable_k for fs in spec.feature_sets)
        else "fixed"
    )
    manifest = emit.build_manifest_dict(
        db_id=spec.id,
        version=spec.version,
        karyoscope_min_version=__version__,
        basename="index/features",
        s=spec.s,
        feature_sets=[fs.name for fs in spec.feature_sets],
        kmer_type=kmer_type,
        roles=spec.roles,
        smoothing=spec.smoothing,
    )
    emit.write_manifest(db_dir / "manifest.yaml", manifest)


def _register_built(db_root: Path, db_dir: Path, version: str, *, force: bool) -> bool:
    """Record the built database in ``installed.json`` (mirrors ``register``)."""
    rel = db_dir.relative_to(db_root)
    state = _installed.load(db_root)
    if db_dir.name in state.databases and not force:
        raise BuildError(
            f"database {db_dir.name!r} is already registered; pass --force to overwrite "
            "the install record."
        )
    _installed.record_install(
        db_root,
        db_dir.name,
        InstalledRecord(
            version=version,
            installed_at=now_iso(),
            source_url="local",
            source_sha256="",
            directory=str(rel),
            registry_doi=None,
        ),
    )
    return True
