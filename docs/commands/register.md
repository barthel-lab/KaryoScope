# karyoscope register

Register a database already present under the database root so the data commands can use it.

## Synopsis

```
karyoscope register [OPTIONS] DATABASE_ID_OR_PATH
```

## Description

`download` is the normal way to install a database: it fetches the tarball, extracts it, validates it, and records it in `installed.json`. But databases are sometimes produced locally — built by hand or copied from another machine — and unpacked into the database root directly. Such a database is valid on disk but invisible to the data commands (`annotate`, `bin`, `scaffold`, `centromeres`, `karyotype`), which resolve databases through `installed.json` only. `register` closes that gap: point it at an already-present database directory (a database id under the root, or a path) and it writes the `installed.json` entry, deriving the id and version from `manifest.yaml` and recording the source as `local`. It validates the on-disk layout (the same check `download` runs) and refuses to clobber an existing entry without `--force`. The database must live inside the database root.

## Options

| Option | Description |
| --- | --- |
| `--db-root DIRECTORY` | Override the database root directory (default: `$KARYOSCOPE_DB` or `~/.karyoscope/db/`). |
| `--force` | Overwrite an existing installed.json entry for this database. |
| `-h`, `--help` | Show this message and exit. |

## Examples

```bash
# Register a database already extracted under the database root, by id
karyoscope register KS_human_CHM13_cytoband

# Register by path, with a custom database root
karyoscope register ./KS_human_CHM13_cytoband --db-root /data/ksdb

# Inspect it afterwards
karyoscope info KS_human_CHM13_cytoband
```

## See also

- [`karyoscope download`](download.md) — fetch and install a pre-built database
- [`karyoscope info`](info.md) — inspect installed databases
