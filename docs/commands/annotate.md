# karyoscope annotate

Annotate sequences against a KaryoScope database, producing per-feature-set BED tracks.

## Synopsis

```
karyoscope annotate -i INPUT [OPTIONS]
```

## Description

`karyoscope annotate` assigns every k-mer in the input to a feature in a single alignment-free pass, by querying the database's KMC index via the bundled `get_featureIDs` helper. It produces one BED per feature set, and by default writes BOTH a "presmoothed" (raw) and a "smoothed" (hierarchy-smoothed) BED for each feature set. k-mers that are not present in the index render as `novel`. Input may be FASTA, FASTQ, or BAM. The expensive k-mer-query step is resumable across reruns, so an interrupted (for example, OOM-killed) run can resume straight into smoothing.

## Options

| Option | Description |
| --- | --- |
| `-i`, `--input FILE` | Input sequence file. Accepts FASTA (`.fasta`/`.fa`/`.fna`, plain or `.gz`), FASTQ (`.fastq`/`.fq`, plain or `.gz`), or BAM (`.bam`). BAM inputs are piped through `samtools fasta` (requires `samtools` on PATH); no intermediate file is written. **[required]** |
| `-o`, `--outdir DIRECTORY` | Directory to write output BEDs into. Default: same directory as `--input`. |
| `--db TEXT` | Database id to use (e.g., `KS_human_CHM13_v2`). Default: the unique installed database if there's exactly one. |
| `--db-root DIRECTORY` | Override the database root directory (default: `$KARYOSCOPE_DB` or `~/.karyoscope/db/`). |
| `--feature-set TEXT` | Restrict output to this feature set. Repeatable. Default: all feature sets declared in the database's manifest. |
| `-t`, `--threads INTEGER` | Threads for both k-mer querying and smoothing. `0` means auto-detect. [default: `0`] |
| `--smooth` / `--no-smooth` | Produce the hierarchy-smoothed BED in addition to the presmoothed BED. [default: `smooth`] |
| `--keep-presmoothed` / `--no-keep-presmoothed` | Keep the presmoothed BED. Pass `--no-keep-presmoothed` to write only the smoothed output. [default: `keep-presmoothed`] |
| `--keep-intermediates` | Keep the combined `.featureIDs.bed` from the C++ step (useful for debugging). |
| `--force` | Regenerate the combined intermediate even if a complete one already exists. By default a rerun reuses a verified combined BED left by a previous (e.g. OOM-killed) run and skips the `get_featureIDs` step, resuming straight into smoothing. A partial file from a killed run is never reused regardless of this flag. |
| `--bgzip` / `--no-bgzip` | bgzip the per-feature-set output BEDs. Pass `--no-bgzip` to write plain `.bed` files. [default: `bgzip`] |
| `--preserve-order` / `--no-preserve-order` | Write output BEDs with sequences in the same order as the input. Pass `--no-preserve-order` for the fastest path when order doesn't matter downstream (typically read data). [default: `preserve-order`] |
| `-h`, `--help` | Show this message and exit. |

## Examples

```bash
# Annotate an assembly with all feature sets (default), plain-text BEDs
karyoscope annotate --input hg002v1.1.fasta.gz --outdir results/ --threads 16 --no-bgzip

# Only the chromosome feature set
karyoscope annotate -i asm.fa --feature-set chromosome -o results/

# Read-level input (FASTQ): faster writes when output order doesn't matter
karyoscope annotate -i reads.fastq.gz -o results/ --no-preserve-order
```

## Output

For each feature set, `annotate` writes up to two BEDs:

- Presmoothed (raw): `<input>.<dbid>.<feature_set>.presmoothed.bed[.gz]`
- Smoothed (hierarchy-smoothed): `<input>.<dbid>.<feature_set>.smoothed.bed[.gz]`

Each BED's 4th column is the human-readable feature name. k-mers absent from the index render as `novel`. The `.gz` suffix is present unless `--no-bgzip` is passed; `--no-bgzip` keeps BEDs as plain text for easy inspection.

Smoothing promotes short noisy intervals (especially short `novel` runs flanked by specific features) to the lowest common ancestor of their flankers in the feature set's hierarchy.

For human-scale inputs, at least 16 threads and ~50 GB RAM are recommended; HG002 runs in ~30 min at `-t 16`.

## See also

- [`karyoscope karyotype`](karyotype.md) — render karyotype SVGs (runs annotate under the hood)
- [`karyoscope bin`](bin.md) — aggregate base-pair BEDs into larger bins
