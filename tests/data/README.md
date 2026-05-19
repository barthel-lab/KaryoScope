# Test fixtures

This directory holds fixtures used by the KaryoScope test suite.

## Files committed to the repo

- **`build_dummy_db.py`** — script that builds the dummy database from
  scratch. Requires `kmc` to be on `$PATH` (or set via the `KMC` env var).
  Run this only when fixtures need to be regenerated (e.g., manifest
  schema changes). CI does **not** run it.

- **`dummy_db.tar.gz`** — gzipped tarball of the dummy database, produced
  by `build_dummy_db.py`. This is the canonical form used by tests.

- **`dummy_db.sha256`** — SHA-256 of `dummy_db.tar.gz`. Embedded in
  `dummy_registry.yaml` and verified during install tests.

- **`dummy_registry.yaml`** — a minimal `registry.yaml` that lists the
  dummy database. The `url` field contains a placeholder string; tests
  rewrite it at runtime to a `file://` URL pointing at the local
  `dummy_db.tar.gz` (see `tests/conftest.py`).

## Files NOT committed

- **`dummy_db/`** — the unpacked dummy database. Produced by
  `build_dummy_db.py` for local inspection and ignored by `.gitignore`.
  Tests that need an unpacked database extract `dummy_db.tar.gz` to a
  temp directory via the `unpacked_dummy_db` pytest fixture.

## What's inside the dummy database

The dummy database is a structurally complete, but biologically
meaningless, KaryoScope database. It contains a real KMC index built from
a few short synthetic DNA sequences, three TSV/text metadata files
(`hierarchy.tsv`, `features.tsv`, `colors.txt`), and a `manifest.yaml`
describing the layout. The whole thing weighs ~5 MB unpacked but
compresses to ~22 KB as a tarball.

It exercises every code path that doesn't require biologically meaningful
k-mers: download, checksum verification, archive extraction, layout
validation, install-state tracking, and uninstallation.

## Regenerating

```bash
cd <repo root>
python tests/data/build_dummy_db.py
```

The script produces deterministic output (sorted entries, zeroed
metadata) so re-running on the same machine should leave `dummy_db.tar.gz`
byte-identical, and the SHA-256 in `dummy_db.sha256` should match. If
anything about the database contents or manifest schema changes, also
update the `sha256` field in `dummy_registry.yaml` to match the new
`dummy_db.sha256`.
