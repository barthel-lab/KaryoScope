"""Fetch, cache, and parse the KaryoScope database registry.

The registry is a YAML file maintained at:

    https://raw.githubusercontent.com/barthel-lab/KaryoScope-registry/main/registry.yaml

It lists all databases available for download. KaryoScope fetches it at
runtime and caches the result for 24 hours (the default TTL) to avoid
hitting the network on every command. Users can force a refresh by passing
``--refresh-registry`` to ``karyoscope download``, or by deleting
``<db_root>/registry_cache.yaml``.

The :func:`load_registry` function is the high-level entry point used by
commands. Tests pass a custom ``registry_url`` (typically a ``file://`` URL
pointing at a fixture) to exercise the parsing logic offline.
"""

from __future__ import annotations

import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from karyoscope._fetch import fetch
from karyoscope.exceptions import RegistryError

#: Default canonical registry URL.
DEFAULT_REGISTRY_URL = (
    "https://raw.githubusercontent.com/barthel-lab/KaryoScope-registry/main/registry.yaml"
)

#: Default cache lifetime, in seconds. 24 hours.
DEFAULT_CACHE_TTL_SECONDS = 24 * 60 * 60

#: Filename inside the database root where the cached registry lives.
CACHE_FILENAME = "registry_cache.yaml"

#: Schema versions of the registry format that this code understands.
_SUPPORTED_SCHEMA_VERSIONS = frozenset({1})


@dataclass
class TaxonomyEntry:
    """A single taxonomic identifier for a database."""

    genus: str
    species: str
    common_name: str | None = None
    ncbi_taxid: int | None = None


@dataclass
class DatabaseEntry:
    """A parsed registry entry describing one downloadable database."""

    id: str
    title: str
    version: str
    karyoscope_min_version: str
    taxonomy: list[TaxonomyEntry]
    feature_sets: list[str]
    size_gb: float
    source: str
    url: str
    sha256: str
    kmer_size: int
    kmer_type: str
    kmer_max_size: int
    is_default: bool = False
    description: str | None = None
    release_date: str | None = None
    doi: str | None = None
    citation: str | None = None
    maintained_by: str | None = None
    tags: list[str] = field(default_factory=list)
    roles: dict[str, str] = field(default_factory=dict)
    recommended_smoothing_window_bp: int | None = None
    # 'official' if from the top-level 'databases' list, 'community' otherwise.
    source_category: str = "official"


@dataclass
class Registry:
    """The full parsed registry."""

    schema_version: int
    databases: list[DatabaseEntry]

    def find(self, db_id: str) -> DatabaseEntry | None:
        """Look up a database by ID. Returns None if not found."""
        for db in self.databases:
            if db.id == db_id:
                return db
        return None

    def default_database(self) -> DatabaseEntry | None:
        """Return the entry marked ``default: true``, if any."""
        for db in self.databases:
            if db.is_default:
                return db
        return None


def _require(d: dict, key: str, where: str) -> object:
    if key not in d:
        raise RegistryError(f"missing required field '{key}' in {where}")
    return d[key]


def _parse_taxonomy(raw: object, where: str) -> list[TaxonomyEntry]:
    if not isinstance(raw, list) or not raw:
        raise RegistryError(f"'taxonomy' in {where} must be a non-empty list")
    out: list[TaxonomyEntry] = []
    for i, entry in enumerate(raw):
        if not isinstance(entry, dict):
            raise RegistryError(f"taxonomy[{i}] in {where} must be a mapping")
        genus = _require(entry, "genus", f"{where}:taxonomy[{i}]")
        species = _require(entry, "species", f"{where}:taxonomy[{i}]")
        if not isinstance(genus, str) or not isinstance(species, str):
            raise RegistryError(f"taxonomy[{i}] in {where}: genus and species must be strings")
        common_name = entry.get("common_name")
        ncbi_taxid = entry.get("ncbi_taxid")
        if common_name is not None and not isinstance(common_name, str):
            raise RegistryError(f"taxonomy[{i}].common_name in {where} must be a string")
        if ncbi_taxid is not None and not isinstance(ncbi_taxid, int):
            raise RegistryError(f"taxonomy[{i}].ncbi_taxid in {where} must be an integer")
        out.append(
            TaxonomyEntry(
                genus=genus, species=species, common_name=common_name, ncbi_taxid=ncbi_taxid
            )
        )
    return out


def _parse_database_entry(raw: dict, source_category: str) -> DatabaseEntry:
    db_id = _require(raw, "id", "<database>")
    if not isinstance(db_id, str):
        raise RegistryError("'id' must be a string")

    where = f"database[{db_id}]"

    def _str(key: str, *, optional: bool = False) -> str | None:
        if key not in raw:
            if optional:
                return None
            raise RegistryError(f"missing required field '{key}' in {where}")
        v = raw[key]
        if v is None and optional:
            return None
        if not isinstance(v, str):
            raise RegistryError(f"'{key}' in {where} must be a string")
        return v

    title = _str("title")
    version = _str("version")
    min_ver = _str("karyoscope_min_version")
    source = _str("source")
    url = _str("url")
    sha256 = _str("sha256")

    size_gb_raw = _require(raw, "size_gb", where)
    if not isinstance(size_gb_raw, (int, float)):
        raise RegistryError(f"'size_gb' in {where} must be a number")
    size_gb = float(size_gb_raw)

    taxonomy = _parse_taxonomy(_require(raw, "taxonomy", where), where)

    fs_raw = _require(raw, "feature_sets", where)
    if not isinstance(fs_raw, list) or not all(isinstance(f, str) for f in fs_raw):
        raise RegistryError(f"'feature_sets' in {where} must be a list of strings")
    if not fs_raw:
        raise RegistryError(f"'feature_sets' in {where} must not be empty")

    kmer_raw = _require(raw, "kmer", where)
    if not isinstance(kmer_raw, dict):
        raise RegistryError(f"'kmer' in {where} must be a mapping")
    kmer_size = _require(kmer_raw, "size", f"{where}:kmer")
    kmer_type = _require(kmer_raw, "type", f"{where}:kmer")
    kmer_max = _require(kmer_raw, "max_size", f"{where}:kmer")
    if not isinstance(kmer_size, int) or kmer_size < 1:
        raise RegistryError(f"'kmer.size' in {where} must be a positive integer")
    if kmer_type not in ("fixed", "variable"):
        raise RegistryError(f"'kmer.type' in {where} must be 'fixed' or 'variable'")
    if not isinstance(kmer_max, int) or kmer_max < 1:
        raise RegistryError(f"'kmer.max_size' in {where} must be a positive integer")

    roles = raw.get("roles", {})
    if not isinstance(roles, dict):
        raise RegistryError(f"'roles' in {where} must be a mapping")
    for k, v in roles.items():
        if not isinstance(k, str) or not isinstance(v, str):
            raise RegistryError(f"'roles' entries in {where} must be string->string")

    tags = raw.get("tags", [])
    if not isinstance(tags, list) or not all(isinstance(t, str) for t in tags):
        raise RegistryError(f"'tags' in {where} must be a list of strings")

    smoothing = raw.get("smoothing", {})
    rec_window = None
    if smoothing:
        if not isinstance(smoothing, dict):
            raise RegistryError(f"'smoothing' in {where} must be a mapping")
        rec_window = smoothing.get("recommended_window_bp")
        if rec_window is not None and not isinstance(rec_window, int):
            raise RegistryError(f"'smoothing.recommended_window_bp' in {where} must be an integer")

    is_default = bool(raw.get("default", False))

    return DatabaseEntry(
        id=db_id,
        title=title or "",
        version=version or "",
        karyoscope_min_version=min_ver or "",
        taxonomy=taxonomy,
        feature_sets=list(fs_raw),
        size_gb=size_gb,
        source=source or "",
        url=url or "",
        sha256=sha256 or "",
        kmer_size=kmer_size,
        kmer_type=kmer_type,
        kmer_max_size=kmer_max,
        is_default=is_default,
        description=_str("description", optional=True),
        release_date=_str("release_date", optional=True),
        doi=_str("doi", optional=True),
        citation=_str("citation", optional=True),
        maintained_by=_str("maintained_by", optional=True),
        tags=list(tags),
        roles=dict(roles),
        recommended_smoothing_window_bp=rec_window,
        source_category=source_category,
    )


def parse_registry(yaml_text: str) -> Registry:
    """Parse the contents of a registry YAML file into a :class:`Registry`."""
    try:
        data = yaml.safe_load(yaml_text)
    except yaml.YAMLError as e:
        raise RegistryError(f"could not parse registry YAML: {e}") from e

    if not isinstance(data, dict):
        raise RegistryError("registry top-level must be a YAML mapping")

    schema = data.get("schema_version")
    if not isinstance(schema, int):
        raise RegistryError("'schema_version' must be an integer")
    if schema not in _SUPPORTED_SCHEMA_VERSIONS:
        raise RegistryError(
            f"unsupported registry schema_version: {schema} "
            f"(supported: {sorted(_SUPPORTED_SCHEMA_VERSIONS)}). "
            "Upgrade KaryoScope to use this registry."
        )

    official = data.get("databases", []) or []
    community = data.get("community_databases", []) or []
    if not isinstance(official, list) or not isinstance(community, list):
        raise RegistryError("'databases' and 'community_databases' must be lists")

    entries: list[DatabaseEntry] = []
    seen_ids: set[str] = set()
    for raw in official:
        if not isinstance(raw, dict):
            raise RegistryError("each entry in 'databases' must be a mapping")
        entry = _parse_database_entry(raw, source_category="official")
        if entry.id in seen_ids:
            raise RegistryError(f"duplicate database id in registry: {entry.id}")
        seen_ids.add(entry.id)
        entries.append(entry)
    for raw in community:
        if not isinstance(raw, dict):
            raise RegistryError("each entry in 'community_databases' must be a mapping")
        entry = _parse_database_entry(raw, source_category="community")
        if entry.id in seen_ids:
            raise RegistryError(f"duplicate database id in registry: {entry.id}")
        seen_ids.add(entry.id)
        entries.append(entry)

    # Enforce at-most-one-default-true.
    defaults = [e for e in entries if e.is_default]
    if len(defaults) > 1:
        ids = ", ".join(d.id for d in defaults)
        raise RegistryError(f"multiple databases marked default=true: {ids}")

    return Registry(schema_version=schema, databases=entries)


def _cache_path(db_root: Path) -> Path:
    return db_root / CACHE_FILENAME


def _cache_is_fresh(cache_file: Path, ttl_seconds: float) -> bool:
    if not cache_file.is_file():
        return False
    age = time.time() - cache_file.stat().st_mtime
    return age < ttl_seconds


def load_registry(
    db_root: Path,
    registry_url: str = DEFAULT_REGISTRY_URL,
    *,
    refresh: bool = False,
    cache_ttl_seconds: float = DEFAULT_CACHE_TTL_SECONDS,
    show_progress: bool = False,
) -> Registry:
    """Load the registry, fetching it from the network if the cache is stale.

    Parameters
    ----------
    db_root
        The KaryoScope database root (e.g., ``~/.karyoscope/db/``). The
        registry cache is written here as ``registry_cache.yaml``.
    registry_url
        Where to fetch the registry from. Defaults to the canonical GitHub
        raw URL. Tests typically pass a ``file://`` URL.
    refresh
        If True, ignore any cached copy and always fetch fresh.
    cache_ttl_seconds
        Maximum cache age in seconds. Defaults to 24 hours.
    show_progress
        Whether to show a tqdm progress bar while fetching the registry.
        Usually False — the registry is small.

    Raises
    ------
    RegistryError
        If the registry cannot be fetched or is malformed.
    """
    db_root.mkdir(parents=True, exist_ok=True)
    cache_file = _cache_path(db_root)

    if not refresh and _cache_is_fresh(cache_file, cache_ttl_seconds):
        return parse_registry(cache_file.read_text())

    # Fetch fresh.
    with tempfile.NamedTemporaryFile(
        mode="wb",
        delete=False,
        dir=db_root,
        prefix=".registry_",
        suffix=".tmp",
    ) as tmp:
        tmp_path = Path(tmp.name)
    try:
        fetch(registry_url, tmp_path, show_progress=show_progress, resume=False)
        # Validate before atomically replacing the cache.
        registry = parse_registry(tmp_path.read_text())
        tmp_path.replace(cache_file)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        # If we have a stale cache, fall back to it on network errors.
        if not refresh and cache_file.is_file():
            return parse_registry(cache_file.read_text())
        raise

    return registry
