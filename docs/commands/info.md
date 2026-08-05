# karyoscope info

Inspect installed databases, or probe a path to see whether it is a KaryoScope database.

## Synopsis

```
karyoscope info [OPTIONS] [TARGET]
```

## Description

When run with no argument, `karyoscope info` lists installed databases along with their version, install date, size, and source. When given a database id, it prints the parsed manifest plus feature-set counts derived from the database's `hierarchy.tsv`, and runs hierarchy and colors validation, printing any issues as warnings without exiting non-zero. When given a filesystem path, it probes whether the path is a KaryoScope database directory or a database archive.

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

# Validate a database archive before publishing it
karyoscope info ./HKS_arabidopsis_ColCEN.tar.gz
```

## Validating an archive

Given a `.tar.gz`, `.tgz`, or `.tar` file, `info` checks the database layout
*inside* the archive. This is what the
[registry](https://github.com/barthel-lab/KaryoScope-registry) asks a
contributor to run before opening a pull request.

```
Path: HKS_arabidopsis_ColCEN.tar.gz
  Type: database archive (435.7 MB)
  Top-level directory: HKS_arabidopsis_ColCEN
  Contents: 12 files, 646.1 MB extracted
  Layout valid: yes
  Database id:  HKS_arabidopsis_ColCEN
  Version:                1.2.0
  ...
```

It reports `Layout valid: NO` when the archive holds more or fewer than one
top-level directory, has no `manifest.yaml`, or is missing a file the manifest
declares — including any of the per-feature-set index files an HKS database
needs. It also warns when the manifest's `id` differs from the top-level
directory name, since `download` installs the directory and then looks it up
by id.

The archive is streamed once and never unpacked: only the manifest and the
TSV/TXT sidecars are read into a temporary directory, and the index files are
checked by name and size alone. Reading a multi-gigabyte archive still takes as
long as decompressing it, but costs no disk space.

Layout problems are reported, not raised — `info` exits 0 either way. Read the
output rather than the exit status.

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

`k-mer` reports the query length, the index type, and the largest queryable
length. On a `fixed` index — the default, and what a priority-resolved database
is necessarily built as — `max` equals `size` and that one length is the only
one `annotate` accepts. `max` is a range you may query within only when `type`
is `variable`. See [build → Fixed-k and variable-k](build.md#fixed-k-and-variable-k).

## See also

- [`karyoscope download`](download.md) — install databases
- [`karyoscope register`](register.md) — register an on-disk database
