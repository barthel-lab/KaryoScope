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
- Top-level verbosity flags ``-v`` / ``-vv`` (raise log level to INFO / DEBUG)
  and ``-q`` (lower to ERROR-only). Logging is reserved for diagnostics
  (URL fetches, cache hits, SHA-256 verifications, subprocess invocations);
  program output remains visible regardless of verbosity. Default is WARNING.
- ``karyoscope.core.external``: small wrapper around :mod:`subprocess`
  with consistent error reporting (:class:`ExternalToolError` includes the
  command, exit code, and a tail of stderr) and a :func:`require_tool`
  helper for "binary not found on PATH" cases with actionable error messages.
- ``karyoscope.core.io.hierarchy``: parser for the ``hierarchy.tsv``
  format used inside KaryoScope databases.
- The ``karyoscope info`` command now provides a real implementation
  rather than a stub. With no argument it lists installed databases
  with version, install date, size, and source. Given a database id, it
  prints the parsed manifest plus feature-set counts derived from the
  database's ``hierarchy.tsv``. Given a filesystem path, it probes
  whether the path is a KaryoScope database directory.
- `native/` directory holding the C++ helper code that KaryoScope uses
  for the performance-critical k-mer querying step. Contains a vendored
  copy of the [KMC](https://github.com/refresh-bio/KMC) 3.2.4 API
  (`native/external/kmc_api/`, GPL-3) with a one-line patch wrapping the
  libstdc++-only `<ext/algorithm>` include in `#ifdef __GLIBCXX__` so
  the source compiles under macOS libc++. Also vendors
  [cxxopts](https://github.com/jarro2783/cxxopts) 3.2.0 as a single
  header (MIT-licensed).
- `native/get_featureIDs/` — the first C++ binary: a multi-threaded
  FASTA/FASTQ reader that queries every k-mer in each input sequence
  against a KMC database and emits run-length-encoded BED records of
  per-position feature ids. Supports plain FASTA/FASTQ, gzipped input,
  and stdin. Adapted from the version in the KaryoScope archive repo;
  the only changes are include paths driven by the new Makefile's `-I`
  flags rather than the previous `../external/...` relative paths.
- `native/get_featureIDs/Makefile` — cross-platform build with
  incremental compilation. Defaults `CXX` to `c++` (the system's
  default C++ compiler — g++ on Linux, Apple Clang on macOS). Locates
  zlib via `pkg-config` with a fall-back to plain `-lz`, and picks up
  `$CONDA_PREFIX` paths when set. Requires C++20.
- `karyoscope.core.io.kmc` — Python wrapper around `get_featureIDs`.
  Resolves the binary in precedence order:
  `$KARYOSCOPE_GET_FEATUREIDS` → `shutil.which("get_featureIDs")` →
  walking up from the package source to find
  `<repo>/native/get_featureIDs/build/get_featureIDs` (works in
  editable installs). Raises `ToolNotFoundError` with actionable
  instructions if none resolve.
- Integration tests (`tests/test_kmc.py`, marked
  `@pytest.mark.integration`) verifying the full pipeline:
  Python wrapper → subprocess → BED output → schema validation. CI
  builds the binary on Linux and macOS in a separate
  `cpp_build_and_integration` job and runs `pytest -m integration`.
- `native/README.md` documenting the layout, vendored code, build
  prerequisites, and the binary-lookup logic.

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
- Default `pytest` invocations now deselect tests marked
  `@pytest.mark.integration`. CI runs them explicitly via
  `pytest -m integration` after building the C++ helper. This keeps the
  unit-test workflow fast and friendly for Python-only contributors,
  while still ensuring the full pipeline is exercised in CI.

### Fixed
- `paths.default_db_root()` now resolves `~` at call time rather than module
  import time, so `$HOME` changes during a process lifetime are honored.
- `CITATION.cff` and `CODE_OF_CONDUCT.md` normalized to a single trailing
  newline (auto-fixed by the new `end-of-file-fixer` pre-commit hook).
- Test fixtures (`unpacked_dummy_db`, `populated_db_root`) no longer crash
  on Python versions older than the `tarfile.TarFile.extractall(filter=...)`
  backport (3.10.12 / 3.11.4 / 3.12+). A shared `_extractall_compat`
  helper in `conftest.py` wraps the call in a `try`/`except TypeError`
  fall-through, mirroring the same pattern that production code in
  `download.py` already uses.
- `native/get_featureIDs/Makefile` now adds `-Wl,-rpath,$(CONDA_PREFIX)/lib`
  to the link flags when `CONDA_PREFIX` is set. Without this, on
  macOS + conda the binary built fine but failed at runtime with
  "dyld: Library not loaded: @rpath/libz.1.dylib ... no LC_RPATH's
  found" — conda-forge dylibs use `@rpath/libfoo.dylib` install-names
  and require the consuming binary to embed an `LC_RPATH`. Conda's own
  compilers do this automatically; the system `c++` doesn't.
  `native/README.md` documents the gotcha for future reference.

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
