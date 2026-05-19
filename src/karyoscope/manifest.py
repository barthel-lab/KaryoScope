"""Per-database manifest parsing and on-disk layout validation.

A KaryoScope database is a directory containing a ``manifest.yaml`` plus
hierarchy, features, colors, and index files. See ``DATABASE_LAYOUT.md`` in
the KaryoScope-registry repository for the canonical specification.

This module reads a database directory and verifies it conforms to that
specification. It does *not* perform deep validity checks (e.g., opening
the KMC index and confirming it queries correctly) — that's the job of the
commands that actually use the database. The goal here is to catch
malformed archives at install time, before they end up in a state where a
later command fails with a confusing error.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

from karyoscope.exceptions import DatabaseLayoutError, ManifestError

#: Supported index types. New backends (e.g. ``hks``) will be added here as
#: they land. See DATABASE_LAYOUT.md for the manifest schema details.
_SUPPORTED_INDEX_TYPES = frozenset({"kmc"})


@dataclass
class IndexSpec:
    """The ``index:`` block of a database manifest.

    For ``type == "kmc"``, ``basename`` is a path relative to the database
    directory; KMC opens ``<basename>.kmc_pre`` and ``<basename>.kmc_suf``.
    """

    type: str
    basename: str | None = None  # for type=kmc


@dataclass
class KmerSpec:
    """The ``kmer:`` block of a database manifest."""

    size: int
    type: str  # "fixed" or "variable"
    max_size: int


@dataclass
class Manifest:
    """A parsed and validated database manifest.

    Attribute names mirror the YAML keys. ``directory`` is set to the
    directory the manifest was loaded from, for convenience when callers
    need to resolve relative paths like ``hierarchy``.
    """

    id: str
    version: str
    karyoscope_min_version: str
    index: IndexSpec
    hierarchy: str
    features: str
    colors: str
    kmer: KmerSpec
    feature_sets: list[str]
    roles: dict[str, str] = field(default_factory=dict)
    smoothing: dict[str, object] = field(default_factory=dict)
    directory: Path | None = None


def _require(d: dict, key: str, where: str) -> object:
    if key not in d:
        raise ManifestError(f"missing required field '{key}' in {where}")
    return d[key]


def parse_manifest(manifest_path: Path) -> Manifest:
    """Parse and validate a ``manifest.yaml`` file.

    Performs schema-level checks (required fields, types, allowed values).
    Does *not* verify that referenced files exist on disk; use
    :func:`validate_database_layout` for that.

    Raises
    ------
    ManifestError
        If the manifest is missing, malformed, or references unsupported features.
    """
    if not manifest_path.is_file():
        raise ManifestError(f"manifest not found: {manifest_path}")

    try:
        data = yaml.safe_load(manifest_path.read_text())
    except yaml.YAMLError as e:
        raise ManifestError(f"could not parse YAML in {manifest_path}: {e}") from e

    if not isinstance(data, dict):
        raise ManifestError(f"top-level YAML in {manifest_path} must be a mapping")

    where = manifest_path.name

    db_id = _require(data, "id", where)
    if not isinstance(db_id, str) or not db_id:
        raise ManifestError(f"'id' in {where} must be a non-empty string")

    version = _require(data, "version", where)
    if not isinstance(version, str):
        raise ManifestError(f"'version' in {where} must be a string")

    min_ver = _require(data, "karyoscope_min_version", where)
    if not isinstance(min_ver, str):
        raise ManifestError(f"'karyoscope_min_version' in {where} must be a string")

    # index block
    index_raw = _require(data, "index", where)
    if not isinstance(index_raw, dict):
        raise ManifestError(f"'index' in {where} must be a mapping")
    index_type = _require(index_raw, "type", f"{where}:index")
    if index_type not in _SUPPORTED_INDEX_TYPES:
        raise ManifestError(
            f"unsupported index type '{index_type}' in {where} "
            f"(supported: {sorted(_SUPPORTED_INDEX_TYPES)})"
        )
    if index_type == "kmc":
        basename = _require(index_raw, "basename", f"{where}:index")
        if not isinstance(basename, str):
            raise ManifestError(f"'index.basename' in {where} must be a string")
        index = IndexSpec(type=index_type, basename=basename)
    else:  # pragma: no cover — guarded by the type check above
        index = IndexSpec(type=index_type)

    # kmer block
    kmer_raw = _require(data, "kmer", where)
    if not isinstance(kmer_raw, dict):
        raise ManifestError(f"'kmer' in {where} must be a mapping")
    kmer_size = _require(kmer_raw, "size", f"{where}:kmer")
    kmer_type = _require(kmer_raw, "type", f"{where}:kmer")
    kmer_max = _require(kmer_raw, "max_size", f"{where}:kmer")
    if not isinstance(kmer_size, int) or kmer_size < 1:
        raise ManifestError(f"'kmer.size' in {where} must be a positive integer")
    if kmer_type not in ("fixed", "variable"):
        raise ManifestError(
            f"'kmer.type' in {where} must be 'fixed' or 'variable', got {kmer_type!r}"
        )
    if not isinstance(kmer_max, int) or kmer_max < 1:
        raise ManifestError(f"'kmer.max_size' in {where} must be a positive integer")
    kmer = KmerSpec(size=kmer_size, type=kmer_type, max_size=kmer_max)

    feature_sets = _require(data, "feature_sets", where)
    if not isinstance(feature_sets, list) or not all(isinstance(f, str) for f in feature_sets):
        raise ManifestError(f"'feature_sets' in {where} must be a list of strings")
    if not feature_sets:
        raise ManifestError(f"'feature_sets' in {where} must not be empty")

    hierarchy = _require(data, "hierarchy", where)
    features = _require(data, "features", where)
    colors = _require(data, "colors", where)
    for k, v in (("hierarchy", hierarchy), ("features", features), ("colors", colors)):
        if not isinstance(v, str):
            raise ManifestError(f"'{k}' in {where} must be a string path")

    roles = data.get("roles", {})
    if not isinstance(roles, dict):
        raise ManifestError(f"'roles' in {where} must be a mapping if present")

    smoothing = data.get("smoothing", {})
    if not isinstance(smoothing, dict):
        raise ManifestError(f"'smoothing' in {where} must be a mapping if present")

    return Manifest(
        id=db_id,
        version=version,
        karyoscope_min_version=min_ver,
        index=index,
        hierarchy=hierarchy,
        features=features,
        colors=colors,
        kmer=kmer,
        feature_sets=feature_sets,
        roles=dict(roles),
        smoothing=dict(smoothing),
        directory=manifest_path.parent,
    )


def validate_database_layout(db_dir: Path) -> Manifest:
    """Validate that ``db_dir`` is a well-formed KaryoScope database.

    Checks:

    * ``manifest.yaml`` is present and parses cleanly.
    * Files referenced from the manifest (hierarchy, features, colors) exist.
    * For ``index.type == "kmc"``: ``<basename>.kmc_pre`` and ``<basename>.kmc_suf``
      both exist.

    Returns the parsed :class:`Manifest` on success. The manifest's
    ``directory`` attribute is set to ``db_dir.resolve()``.

    Raises
    ------
    ManifestError
        For manifest-level problems.
    DatabaseLayoutError
        If files are missing or paths escape ``db_dir``.
    """
    db_dir = db_dir.resolve()
    manifest = parse_manifest(db_dir / "manifest.yaml")
    manifest.directory = db_dir

    def _check_exists(rel: str, kind: str) -> None:
        # Defend against path-traversal manifests that point outside db_dir.
        full = (db_dir / rel).resolve()
        try:
            full.relative_to(db_dir)
        except ValueError as e:
            raise DatabaseLayoutError(
                f"manifest '{kind}' path escapes database directory: {rel}"
            ) from e
        if not full.is_file():
            raise DatabaseLayoutError(f"missing {kind} file: {rel} (looked at {full})")

    _check_exists(manifest.hierarchy, "hierarchy")
    _check_exists(manifest.features, "features")
    _check_exists(manifest.colors, "colors")

    if manifest.index.type == "kmc":
        assert manifest.index.basename is not None  # guaranteed by parse_manifest
        _check_exists(manifest.index.basename + ".kmc_pre", "KMC index (.kmc_pre)")
        _check_exists(manifest.index.basename + ".kmc_suf", "KMC index (.kmc_suf)")

    return manifest
