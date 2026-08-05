# karyoscope download

Download and manage pre-built KaryoScope databases.

## Synopsis

```
karyoscope download [OPTIONS] [DATABASE_ID ...]
```

## Description

Databases are resolved through an `installed.json` file kept under the database root, and running `karyoscope download` is the normal way to install a pre-built database. Several "action" flags (`--list`, `--info`, `--status`, `--remove`) also make this command the entry point for discovering and managing databases. The default database root is `$KARYOSCOPE_DB` or `~/.karyoscope/db/`.

Downloads are SHA-256 verified; HTTP downloads stage to a `.part` file and resume via HTTP Range requests if interrupted. If the download completes but the *install* then fails — extraction runs out of space, the run is interrupted — the verified archive is kept at `<db_root>/.<DATABASE_ID>.tar.gz` and any partially-extracted directory is removed. Re-running `karyoscope download <ID>` re-verifies that archive and extracts it directly, so a failed install costs one extraction rather than a second full transfer. Delete the staged file to reclaim its space instead. The registry is cached for 24 hours at `<db_root>/registry_cache.yaml`, and a stale cache is used on transient network failure. The current default database is `KS_human_CHM13_v2` (17.2 GB installed).

### Disk space

A database has two sizes, and for the HKS databases they differ substantially:

| Database | Download (`.tar.gz`) | On disk after install | **Free space to install** |
|---|---|---|---|
| `KS_human_CHM13_v2` (default) | 16.3 GB | 17.2 GB | **~34 GB** |
| `HKS_human_CHM13_v2` | 13.3 GB | 22.7 GB | **~36 GB** |
| `HKS_arabidopsis_ColCEN` | 0.46 GB | 0.68 GB | **~1.1 GB** |

Installing needs the sum of the two columns, because the archive is only deleted once extraction succeeds. `download` verifies that up front and refuses with the exact shortfall rather than filling the disk part-way through extraction; `--no-space-check` overrides it. Reinstalling over an existing copy credits the space that copy will free, and a staged or partially-downloaded archive is credited too, so a retry doesn't require room for the archive twice.

Sizes come from the registry's `size_gb` (extracted) and `download_size_gb` (archive) fields, in decimal GB. `df -h` reports binary GiB, so 36 GB appears there as 34 GiB. An entry that predates `download_size_gb` falls back to `size_gb` for both, and the resulting figure is labelled as an estimate.

## Options

| Option | Description |
| --- | --- |
| `--list` | List databases available in the registry, then exit. |
| `--info DATABASE_ID` | Show detailed information about a single database, then exit. |
| `--status` | Show locally installed databases, then exit. |
| `--remove DATABASE_ID` | Uninstall a locally installed database, then exit. |
| `--organism NAME` | Filter `--list` by common name / genus / species (case-insensitive substring). |
| `--tag TAG` | Filter `--list` by tag (e.g. 'test', 'reference'). |
| `--community` | Include community-contributed databases in listings. |
| `--db-root DIRECTORY` | Override the database root directory (default: `$KARYOSCOPE_DB` or `~/.karyoscope/db/`). |
| `--registry-url URL` | Override the registry URL (advanced; for testing or private registries). |
| `--refresh-registry` | Force a fresh fetch of the registry, ignoring any cached copy. |
| `--force` | Re-download and re-install even if the database is already present. |
| `--no-checksum` | Skip SHA-256 verification (not recommended; useful for debugging). |
| `--no-space-check` | Install even if the database root looks too small to hold the archive and its extracted contents. Only useful when the registry's declared sizes are wrong for your copy of the database. |
| `-y, --yes` | Assume 'yes' to interactive prompts (e.g. `--remove`). |
| `-q, --quiet` | Suppress progress bars. |
| `-h, --help` | Show this message and exit. |

## Examples

```bash
# Install the registry's default database
karyoscope download

# List databases available in the registry
karyoscope download --list

# Show details for one database
karyoscope download --info KS_human_CHM13_v2

# Show locally installed databases
karyoscope download --status

# Uninstall a database
karyoscope download --remove KS_human_CHM13_v2
```

## See also

- [`karyoscope register`](register.md) — register a database already present on disk
- [`karyoscope info`](info.md) — inspect installed databases
