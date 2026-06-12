# karyoscope info

Inspect installed databases, or probe a path to see whether it is a KaryoScope database.

## Synopsis

```
karyoscope info [OPTIONS] [TARGET]
```

## Description

When run with no argument, `karyoscope info` lists installed databases along with their version, install date, size, and source. When given a database id, it prints the parsed manifest plus feature-set counts derived from the database's `hierarchy.tsv`, and runs hierarchy and colors validation, printing any issues as warnings without exiting non-zero. When given a filesystem path, it probes whether the path is a KaryoScope database directory.

## Options

| Option | Description |
| --- | --- |
| `--db-root DIRECTORY` | Override the database root directory (default: `$KARYOSCOPE_DB` or `~/.karyoscope/db/`). |
| `-h`, `--help` | Show this message and exit. |

## Examples

```bash
# List installed databases
karyoscope info

# Inspect a specific installed database (manifest, feature sets, validation)
karyoscope info KS_human_CHM13_v2

# Probe whether a directory is a KaryoScope database
karyoscope info ./some/path
```

## Example output

```
KS_dummy_test_v1
  Path:                   /.../db/KS_dummy_test_v1
  Installed at:           2026-06-12T05:25:08Z
  Source URL:             local
  Version:                1.0.0
  KaryoScope min version: 0.1.0
  k-mer:                  size=21, type=fixed, max=21
  Index type:             kmc
  Index basename:         index/features
  Roles:                  chromosome_assignment=chromosome
  Size on disk:           5.0 MB
  Feature sets:
    chromosome: 3 edges
    region: 6 edges
```

## See also

- [`karyoscope download`](download.md) — install databases
- [`karyoscope register`](register.md) — register an on-disk database
