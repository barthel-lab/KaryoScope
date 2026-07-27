"""Tests for :mod:`karyoscope.registry`."""

from __future__ import annotations

import time
from pathlib import Path

import pytest

from karyoscope.exceptions import RegistryError
from karyoscope.registry import (
    CACHE_FILENAME,
    load_registry,
    parse_registry,
)

_MINIMAL_REGISTRY = """\
schema_version: 1
databases:
  - id: KS_test_v1
    title: Test database
    version: "1.0.0"
    karyoscope_min_version: "0.1.0"
    default: true
    taxonomy:
      - genus: Homo
        species: sapiens
        common_name: human
    size_gb: 17.0
    feature_sets: [chromosome, centromere]
    kmer:
      size: 21
      type: fixed
      max_size: 21
    source: synthetic
    url: https://example.com/x.tar.gz
    sha256: "0000000000000000000000000000000000000000000000000000000000000000"
    tags: [reference]
community_databases: []
"""


def test_parse_registry_minimal() -> None:
    reg = parse_registry(_MINIMAL_REGISTRY)
    assert reg.schema_version == 1
    assert len(reg.databases) == 1
    db = reg.databases[0]
    assert db.id == "KS_test_v1"
    assert db.is_default is True
    assert db.size_gb == pytest.approx(17.0)
    assert db.taxonomy[0].common_name == "human"
    assert db.kmer_size == 21
    assert db.tags == ["reference"]
    assert db.source_category == "official"


def test_parse_registry_finds_default() -> None:
    reg = parse_registry(_MINIMAL_REGISTRY)
    default = reg.default_database()
    assert default is not None
    assert default.id == "KS_test_v1"


def test_parse_registry_find_by_id() -> None:
    reg = parse_registry(_MINIMAL_REGISTRY)
    assert reg.find("KS_test_v1") is not None
    assert reg.find("nope") is None


def test_parse_registry_rejects_unknown_schema() -> None:
    bad = _MINIMAL_REGISTRY.replace("schema_version: 1", "schema_version: 999")
    with pytest.raises(RegistryError, match="schema_version"):
        parse_registry(bad)


def test_parse_registry_rejects_duplicate_ids() -> None:
    duped = """\
schema_version: 1
databases:
  - id: KS_test_v1
    title: First
    version: "1.0.0"
    karyoscope_min_version: "0.1.0"
    taxonomy:
      - genus: Homo
        species: sapiens
    size_gb: 1.0
    feature_sets: [chromosome]
    kmer: {size: 21, type: fixed, max_size: 21}
    source: synthetic
    url: https://example.com/x.tar.gz
    sha256: "0"
  - id: KS_test_v1
    title: Duplicate
    version: "1.0.0"
    karyoscope_min_version: "0.1.0"
    taxonomy:
      - genus: Homo
        species: sapiens
    size_gb: 1.0
    feature_sets: [chromosome]
    kmer: {size: 21, type: fixed, max_size: 21}
    source: synthetic
    url: https://example.com/y.tar.gz
    sha256: "0"
community_databases: []
"""
    with pytest.raises(RegistryError, match="duplicate"):
        parse_registry(duped)


def test_parse_registry_rejects_multiple_defaults() -> None:
    text = """\
schema_version: 1
databases:
  - id: KS_test_v1
    title: First
    version: "1.0.0"
    karyoscope_min_version: "0.1.0"
    default: true
    taxonomy:
      - genus: Homo
        species: sapiens
    size_gb: 1.0
    feature_sets: [chromosome]
    kmer: {size: 21, type: fixed, max_size: 21}
    source: synthetic
    url: https://example.com/x.tar.gz
    sha256: "0"
  - id: KS_test_v2
    title: Second default
    version: "1.0.0"
    karyoscope_min_version: "0.1.0"
    default: true
    taxonomy:
      - genus: Homo
        species: sapiens
    size_gb: 1.0
    feature_sets: [chromosome]
    kmer: {size: 21, type: fixed, max_size: 21}
    source: synthetic
    url: https://example.com/z.tar.gz
    sha256: "0"
community_databases: []
"""
    with pytest.raises(RegistryError, match="default=true"):
        parse_registry(text)


def test_parse_registry_marks_community_entries() -> None:
    text = """\
schema_version: 1
databases:
  - id: KS_test_v1
    title: Official
    version: "1.0.0"
    karyoscope_min_version: "0.1.0"
    taxonomy:
      - genus: Homo
        species: sapiens
    size_gb: 1.0
    feature_sets: [chromosome]
    kmer: {size: 21, type: fixed, max_size: 21}
    source: synthetic
    url: https://example.com/x.tar.gz
    sha256: "0"
community_databases:
  - id: KS_community_v1
    title: Community database
    version: "1.0.0"
    karyoscope_min_version: "0.1.0"
    taxonomy:
      - genus: Mus
        species: musculus
    size_gb: 2.0
    feature_sets: [chromosome]
    kmer: {size: 21, type: fixed, max_size: 21}
    source: synthetic
    url: https://example.com/c.tar.gz
    sha256: "0"
"""
    reg = parse_registry(text)
    community = [d for d in reg.databases if d.source_category == "community"]
    assert len(community) == 1
    assert community[0].id == "KS_community_v1"


# --- load_registry (cache behavior) ------------------------------------


def test_load_registry_fetches_when_cache_missing(tmp_path: Path) -> None:
    registry_file = tmp_path / "registry.yaml"
    registry_file.write_text(_MINIMAL_REGISTRY)
    db_root = tmp_path / "db_root"

    reg = load_registry(db_root, registry_file.absolute().as_uri())
    assert reg.find("KS_test_v1") is not None

    # Cache should now exist.
    assert (db_root / CACHE_FILENAME).is_file()


def test_load_registry_uses_cache_when_fresh(tmp_path: Path) -> None:
    db_root = tmp_path / "db_root"
    db_root.mkdir()
    cache = db_root / CACHE_FILENAME
    cache.write_text(_MINIMAL_REGISTRY)

    # Point at a registry file that does NOT exist; we should never touch it
    # because the cache is fresh.
    missing = tmp_path / "definitely_missing.yaml"
    reg = load_registry(db_root, missing.absolute().as_uri())
    assert reg.find("KS_test_v1") is not None


def test_load_registry_refreshes_when_cache_stale(tmp_path: Path) -> None:
    db_root = tmp_path / "db_root"
    db_root.mkdir()
    cache = db_root / CACHE_FILENAME
    cache.write_text(_MINIMAL_REGISTRY)
    # Age the cache by 48 hours.
    old = time.time() - (48 * 3600)
    import os

    os.utime(cache, (old, old))

    # Provide a new registry on disk.
    new_registry_text = _MINIMAL_REGISTRY.replace("KS_test_v1", "KS_test_refreshed")
    new_registry_file = tmp_path / "new_registry.yaml"
    new_registry_file.write_text(new_registry_text)

    reg = load_registry(db_root, new_registry_file.absolute().as_uri())
    assert reg.find("KS_test_refreshed") is not None
    assert reg.find("KS_test_v1") is None


def test_load_registry_refresh_flag_bypasses_fresh_cache(tmp_path: Path) -> None:
    db_root = tmp_path / "db_root"
    db_root.mkdir()
    cache = db_root / CACHE_FILENAME
    cache.write_text(_MINIMAL_REGISTRY)

    new_registry_text = _MINIMAL_REGISTRY.replace("KS_test_v1", "KS_forced_refresh")
    new_registry_file = tmp_path / "new_registry.yaml"
    new_registry_file.write_text(new_registry_text)

    reg = load_registry(db_root, new_registry_file.absolute().as_uri(), refresh=True)
    assert reg.find("KS_forced_refresh") is not None


def test_load_registry_falls_back_to_stale_on_network_failure(tmp_path: Path) -> None:
    db_root = tmp_path / "db_root"
    db_root.mkdir()
    cache = db_root / CACHE_FILENAME
    cache.write_text(_MINIMAL_REGISTRY)
    old = time.time() - (48 * 3600)
    import os

    os.utime(cache, (old, old))

    # Point at a non-existent file. load_registry should fall back to the
    # stale cache rather than raising.
    missing = tmp_path / "missing.yaml"
    reg = load_registry(db_root, missing.absolute().as_uri())
    assert reg.find("KS_test_v1") is not None


def test_load_registry_raises_when_no_cache_and_fetch_fails(tmp_path: Path) -> None:
    db_root = tmp_path / "db_root"  # no cache will exist
    missing = tmp_path / "missing.yaml"
    # Either FetchError (the source URL isn't reachable) or RegistryError
    # (after a partial cache write got cleaned up) is an acceptable outcome.
    from karyoscope.exceptions import FetchError

    with pytest.raises((FetchError, RegistryError)):
        load_registry(db_root, missing.absolute().as_uri())


def test_load_registry_with_real_test_registry(dummy_registry_url: str, tmp_path: Path) -> None:
    """End-to-end: the actual test-registry fixture parses cleanly."""
    db_root = tmp_path / "db_root"
    reg = load_registry(db_root, dummy_registry_url)
    assert reg.find("KS_dummy_test_v1") is not None
    assert reg.default_database().id == "KS_dummy_test_v1"


# --- size accounting --------------------------------------------------


def test_size_gb_is_the_installed_size_and_download_falls_back_to_it() -> None:
    """An entry without ``download_size_gb`` still yields usable numbers.

    Older registry entries (and any cached copy fetched before the field
    existed) only carry ``size_gb``. Rather than fail, both accessors
    collapse onto it -- and say so, so callers can flag the peak figure as
    a guess rather than a measurement.
    """
    entry = parse_registry(_MINIMAL_REGISTRY).find("KS_test_v1")
    assert entry.installed_size_bytes == 17_000_000_000
    assert entry.download_size_bytes == 17_000_000_000
    assert entry.download_size_is_declared is False


def test_explicit_download_size_is_used_when_declared() -> None:
    """The real HKS entry: a 13.3 GB archive that unpacks to 22.7 GB.

    Reporting only the download size is what let a user with 12 GB free
    start a 25-minute install that could never have finished.
    """
    text = _MINIMAL_REGISTRY.replace(
        "    size_gb: 17.0\n", "    size_gb: 22.7\n    download_size_gb: 13.3\n"
    )
    entry = parse_registry(text).find("KS_test_v1")
    assert entry.download_size_is_declared is True
    assert entry.download_size_bytes == 13_300_000_000
    assert entry.installed_size_bytes == 22_700_000_000
    # Both coexist on disk until extraction succeeds.
    assert entry.peak_install_bytes == 36_000_000_000


def test_download_size_gb_must_be_a_non_negative_number() -> None:
    for bad in ('"lots"', "-1"):
        text = _MINIMAL_REGISTRY.replace(
            "    size_gb: 17.0\n", f"    size_gb: 17.0\n    download_size_gb: {bad}\n"
        )
        with pytest.raises(RegistryError, match="download_size_gb"):
            parse_registry(text)
