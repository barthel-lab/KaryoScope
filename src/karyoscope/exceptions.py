"""Exception hierarchy for KaryoScope.

All KaryoScope errors derive from :class:`KaryoscopeError`. This makes it
straightforward for the CLI layer to catch known-bad-but-non-fatal conditions
and present them to the user as clean ``click.ClickException`` messages
rather than tracebacks.
"""

from __future__ import annotations


class KaryoscopeError(Exception):
    """Base class for all KaryoScope-specific errors."""


class RegistryError(KaryoscopeError):
    """Problems fetching or parsing the database registry."""


class ManifestError(KaryoscopeError):
    """Problems parsing or validating a per-database ``manifest.yaml``."""


class DatabaseLayoutError(KaryoscopeError):
    """A database directory or tarball does not conform to the expected layout."""


class DatabaseNotFoundError(KaryoscopeError):
    """A requested database is not present in the registry or not installed locally."""


class ChecksumError(KaryoscopeError):
    """An archive's SHA-256 did not match the value declared in the registry."""


class FetchError(KaryoscopeError):
    """Network or filesystem error while fetching a resource."""


class UnsupportedSchemeError(FetchError):
    """A URL scheme other than ``http(s)://`` or ``file://`` was passed."""


class IncompatibleVersionError(KaryoscopeError):
    """A database requires a newer KaryoScope version than the one installed."""


class BinError(KaryoscopeError):
    """Problems aggregating annotation BED records into fixed-size bins."""


class ScaffoldError(KaryoscopeError):
    """Problems classifying, orienting, or rewriting per-input scaffolded outputs."""


class CentromereError(KaryoscopeError):
    """Problems extracting centromere coordinates from binned scaffolded BEDs."""
