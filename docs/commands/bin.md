# karyoscope bin

Aggregate base-pair annotation BEDs into larger fixed-size bins.

## Synopsis

```
karyoscope bin -i INPUT -o OUTPUT -b BIN_SIZE [OPTIONS]
```

## Description

`karyoscope bin` aggregates a coordinate-sorted per-base (or run-length-encoded) BED into fixed-size windows, labelling each bin with the winning feature via a deterministic three-rule selection. First, the base pairs of overlap are summed per feature within the bin. Second, if a leaf-feature set is in play (via `--db` and `--feature-set`), leaf features compete first. Third, the `novel` sentinel only wins if it covers a strict majority of the bin. Adjacent bins with the same winning feature are coalesced. Input must be sorted by chrom then start; `-` means stdin/stdout, and `.gz` is handled transparently.

## Options

| Option | Description |
| --- | --- |
| `-i, --input FILE` | Input BED file (sorted by chrom then start). Use `-` for stdin. `.gz` supported. **[required]** |
| `-o, --output FILE` | Output BED file. Use `-` for stdout. Output is gzipped iff the path ends in `.gz`. **[required]** |
| `-b, --bin-size INTEGER` | Bin size in base pairs (e.g. `1000000` for 1 Mb). **[required]** |
| `--db TEXT` | Database id whose `hierarchy.tsv` defines the leaf-feature set. Default: the unique installed database if exactly one is installed and `--feature-set` is given. |
| `--db-root DIRECTORY` | Override the database root directory (default: `$KARYOSCOPE_DB` or `~/.karyoscope/db/`). |
| `--feature-set TEXT` | Feature set to use for leaf prioritisation. Required when `--db` is given (or implied). |
| `-t, --threads INTEGER` | Per-sequence-chunk parallelism. `0` = auto (`os.cpu_count()`). `1` = single-threaded. Stdin/stdout I/O is always single-threaded regardless of this flag. **[default: 1]** |
| `-h, --help` | Show this message and exit. |

Passing `--db` without `--feature-set` is an error: the leaf set is per-feature-set, so one without the other is ambiguous. A bare invocation with no `--db`/`--feature-set` skips leaf prioritisation entirely.

## Examples

```bash
# Bin a per-base BED into 1 Mb windows
karyoscope bin -i sample.region.smoothed.bed.gz -o sample.region.binned.bed.gz -b 1000000

# With leaf prioritisation for a specific feature set
karyoscope bin -i sample.region.smoothed.bed -o sample.region.binned.bed -b 1000000 --db KS_human_CHM13_v2 --feature-set region

# Stream via stdin/stdout
cat sample.bed | karyoscope bin -i - -o - -b 100000 > binned.bed
```

## See also

- [`karyoscope annotate`](annotate.md) — produces the per-base BEDs that `bin` aggregates
- [`karyoscope karyotype`](karyotype.md) — bins internally as part of rendering
