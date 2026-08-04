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
| `--k INTEGER` | Query k-mer length. Defaults to the database's k. Only a variable-k HKS index (built with `karyoscope build --variable-k`) accepts a value other than its k — use it for a k-sweep; on a fixed-k index any other value is an error. A database built with priorities is always fixed-k. Outputs are tagged `.k<k>` so runs into one directory don't collide. |
| `--smooth` / `--no-smooth` | Produce the hierarchy-smoothed BED in addition to the presmoothed BED. [default: `smooth`] |
| `--keep-presmoothed` / `--no-keep-presmoothed` | Keep the presmoothed BED. Pass `--no-keep-presmoothed` to write only the smoothed output. [default: `keep-presmoothed`] |
| `--keep-intermediates` | Keep the combined `.featureIDs.bed` from the C++ step (useful for debugging). |
| `--force` | Regenerate the combined intermediate even if a complete one already exists. By default a rerun reuses a verified combined BED left by a previous (e.g. OOM-killed) run and skips the `get_featureIDs` step, resuming straight into smoothing. A partial file from a killed run is never reused regardless of this flag. |
| `--bgzip` / `--no-bgzip` | bgzip the per-feature-set output BEDs. Pass `--no-bgzip` to write plain `.bed` files. Note this shrinks the final output but not peak disk usage: compression runs after every BED has been written in full. [default: `bgzip`] |
| `--no-space-check` | Skip the up-front free-space check on `--outdir`. See [Disk space](#disk-space). |
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

For human-scale inputs, use at least 16 threads. Memory to request depends on the backend and input shape. With the **HKS backend**, the peak is set by the index rather than the input — the shared base index plus one feature set's labeling at a time — so it is ~10 GB whether you annotate a single haplotype or a combined diploid assembly such as HG002 v1.1 (request ≥ 16 GB). The **KMC backend** peaks at ~30–35 GB (request ≥ 50 GB). HG002 runs in ~20–30 min at `-t 16`.

## Progress output

`annotate` reports what it is doing as it goes, so a long run is distinguishable from a hung one:

```
Annotating hg002v1.1.fasta.gz against HKS_human_CHM13_v2
  6 feature set(s), 16 thread(s), ~34 GB estimated output
  [1/6] chromosome    4m05s
  [2/6] region        3m45s
  ...
  bgzip (12 file(s))  1m31s
Wrote:
  ...
```

The milestone shape follows the backend. HKS queries each feature set in sequence, so it reports `[i/N]` per set. KMC runs one combined query and then smooths every set in a single streaming pass, so it reports named phases instead (`k-mer query`, `smoothing 6 feature set(s)`).

Pass `-q` (before the subcommand: `karyoscope -q annotate ...`) to suppress the narration; the closing `Wrote:` block still prints, since that is the command's result. Detailed per-step timings remain on the logging channel — `karyoscope -v annotate ...` for those.

## Disk space

Output is large. Budget roughly **0.8 GB per feature set per Gbp of input**, plus one intermediate the size of the largest single feature set's presmoothed BED:

| Input | Feature sets | Output BEDs (uncompressed) | Peak including intermediate |
|---|---|---|---|
| Single human haplotype (~3.1 Gbp) | 6 | ~15 GB | ~18 GB |
| Diploid assembly, e.g. HG002 v1.1 (~6.0 Gbp) | 6 | ~29 GB | ~34 GB |

Measured on HG002 v1.1 against `HKS_human_CHM13_v2`: 21.7 GB of presmoothed BED, 7.3 GB of smoothed BED, and a 5.4 GB lookup TSV live at the peak. The BED content is the same for either backend, so the figures hold for both.

`--bgzip` (the default) does **not** lower the peak — every BED is written in full before the compression pass starts. What does lower it:

- `--feature-set` to annotate a subset at a time
- `--no-keep-presmoothed` (drops the larger of the two variants, ~0.6 of the 0.8 GB) or `--no-smooth`
- splitting a diploid assembly into per-haplotype runs

`annotate` estimates this footprint before starting and refuses if `--outdir` can't hold it, so a full disk is caught in a second rather than after the k-mer query. The estimate uses an exact base count from a sibling `.fai` when one exists, and otherwise scales the input file size. Pass `--no-space-check` if the estimate is wrong for your input.

## See also

- [`karyoscope karyotype`](karyotype.md) — render karyotype SVGs (runs annotate under the hood)
- [`karyoscope bin`](bin.md) — aggregate base-pair BEDs into larger bins
