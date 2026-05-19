# Changelog

All notable changes to KaryoScope will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial project scaffold with command-line interface skeleton.
- `karyoscope` console entry point with subcommand dispatch.
- Foundational documentation: `README.md`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`,
  and issue/PR templates.
- GPL-3.0 license.
- Python packaging via `hatchling`, requiring Python ≥3.10.
- GitHub Actions CI workflow running lint, format, and tests on every push.
- Integration with the `KaryoScope-registry` for pre-built database discovery.
- Full implementation of the `karyoscope download` command, supporting installing
  databases by id or as the registry default, listing available databases
  (`--list`) with `--organism` / `--tag` / `--community` filters, inspecting
  individual entries (`--info ID`), showing locally installed databases
  (`--status`), and uninstalling them (`--remove ID`).
- Registry fetching with a 24-hour TTL cache at `<db_root>/registry_cache.yaml`
  and graceful fall-back to a stale cache on transient network errors. The
  default registry URL points at `barthel-lab/KaryoScope-registry`; private or
  test registries can be supplied via `--registry-url`.
- SHA-256 verification of downloaded tarballs, with `--no-checksum` for
  debugging. HTTP downloads stage to a `.part` file and resume via `Range`
  requests when interrupted.
- Safe tarball extraction that refuses entries containing `..`, symlinks or
  other special types, or files outside the database's expected top-level
  directory.
- Per-database manifest schema (`manifest.yaml`) and on-disk layout validator.
- `installed.json` for tracking which databases are installed locally, with
  atomic writes and recovery from corrupted files.
- Test fixtures: a 22 KB dummy database tarball containing a real (tiny) KMC
  index, plus a build script for regenerating it.

### Changed
- Lint and format tooling is now managed by [`pre-commit`](https://pre-commit.com).
  Tool versions are pinned in `.pre-commit-config.yaml` (currently ruff 0.15.13
  and pre-commit-hooks v6.0.0); CI runs `pre-commit run --all-files`, so
  developers and CI use identical tool versions. See `CONTRIBUTING.md` for
  setup instructions.
- `ruff` removed from `pyproject.toml`'s `[dev]` extras since pre-commit
  manages its own isolated installations.
- `tests/conftest.py` now surfaces a clear actionable error when pytest is
  launched from a Python interpreter where `karyoscope` is not installed
  (a common pitfall on macOS where multiple Python installations coexist).
  `CONTRIBUTING.md` updated to recommend `python -m pytest` over bare
  `pytest` to sidestep the issue.

### Fixed
- `paths.default_db_root()` now resolves `~` at call time rather than module
  import time, so `$HOME` changes during a process lifetime are honored.
- `CITATION.cff` and `CODE_OF_CONDUCT.md` normalized to a single trailing
  newline (auto-fixed by the new `end-of-file-fixer` pre-commit hook).

<!--
Use the following sections when adding entries:

### Added
### Changed
### Deprecated
### Removed
### Fixed
### Security

For releases, copy the [Unreleased] section to a new heading like:

## [1.0.0] - 2026-MM-DD
-->

[Unreleased]: https://github.com/barthel-lab/KaryoScope/compare/HEAD...HEAD
