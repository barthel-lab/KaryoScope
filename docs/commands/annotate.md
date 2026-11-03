# karyoscope annotate

Annotate sequences against a KaryoScope database, producing per-feature-set BED tracks.

## Synopsis

```
karyoscope annotate -i INPUT [OPTIONS]
```

## Description

`karyoscope annotate` assigns every k-mer in the input to a feature in a single alignment-free pass, by querying the database's index. The backend follows the database: an HKS-index database (e.g. `HKS_human_CHM13_v2`) is queried with `hks`, and a KMC-index database with the bundled `get_featureIDs` helper. It produces one BED per feature set, and by default writes BOTH a "presmoothed" (raw) and a "smoothed" (hierarchy-smoothed) BED for each feature set. k-mers that are not present in the index render as `novel`. Input may be FASTA, FASTQ, BAM or CRAM. The expensive k-mer-query step is resumable across reruns, so an interrupted (for example, OOM-killed) run can resume straight into smoothing.

## Options

| Option | Description |
| --- | --- |
| `-i`, `--input FILE` | Input sequence file. Accepts FASTA (`.fasta`/`.fa`/`.fna`, plain or `.gz`), FASTQ (`.fastq`/`.fq`, plain or `.gz`), BAM (`.bam`) or CRAM (`.cram`). BAM/CRAM inputs are converted with `samtools fasta` (requires `samtools` on PATH); CRAM also requires `--reference`. **[required]** |
| `--reference FILE` | Reference FASTA a CRAM input was aligned against. **Required for `.cram`**, ignored otherwise. See [CRAM input](#cram-input). |
| `--query-names` / `--no-query-names` | Identify output sequences by name rather than ordinal rank. Assemblies default to names; **refused for read-level input**, which always uses ranks. See [Paired-end reads](#paired-end-reads). |
| `--query-names-sidecar` | Also write `<outdir>/<input>.query_names.txt.gz` for each input: every record's name, one per line, in the order the query assigns ranks (line N+1 is rank N). This is how rank-identified read-level output is joined back to read names. For BAM/CRAM it is teed off the decode `annotate` performs anyway, on either backend, so it is nearly free; for FASTA/FASTQ input it is one extra streaming read of the file. See [Query names](#query-names). |
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
| `--format [bed\|kst]` | Output format. `bed` (the default) is unchanged. `kst` stores the presmoothed tracks as compressed `.kst` files instead, 12–44× smaller than the bgzipped BEDs, and stores no smoothed track: it is derived on read. Requires `karyozip` on PATH. See [Compact output: `--format kst`](#compact-output---format-kst). [default: `bed`] |
| `--kst-effort [low\|medium]` | karyozip compression effort for `--format kst`. `low` writes as fast as bgzip and keeps single-read queries at ~0.02 s; `medium` is about half the size, ~10× slower to write, and ~1.4 s per single-read query. [default: `low`] |
| `--smoothed-bed` | With `--format kst`, also write each set's smoothed BED (`<input>.<dbid>.<set>.smoothed.bed[.gz]`, bgzipped unless `--no-bgzip`), derived from the `.kst` with `karyozip decompress --smooth`. |
| `--no-resource-check` | Skip the up-front disk **and memory** checks. Disk: the output footprint on `--outdir`, see [Disk space](#disk-space). Memory: HKS databases must hold their index resident (~10 GB for a human database), so an under-allocated run is killed by the OS mid-query; the check reads the requirement from the index files themselves. |
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

# Compact output for a cohort: one merged .kst per input, no BEDs
karyoscope annotate -i sample.fa.gz --db HKS_human_CHM13_v2 -o results/ --format kst

# Read-level input, keeping the rank -> read-name mapping next to the output
karyoscope annotate -i tumor.cram --reference GRCh38.fasta -o results/ --query-names-sidecar
```

## CRAM input

CRAM stores bases as a diff against the reference used for alignment, so it cannot be
decoded without that same FASTA. Pass it with `--reference`:

```bash
karyoscope annotate -i tumor.cram --reference GRCh38.fasta \
                    --db HKS_human_CHM13_cytoband --feature-set cytoband -o results/
```

`--reference` is the **alignment** reference and has nothing to do with `--db`. `annotate`
is alignment-free — it only ever sees read sequences — so annotating GRCh38-aligned reads
against a CHM13 database is normal, not a mismatch.

Omitting `--reference` is refused up front rather than left to htslib, which would
otherwise resolve each contig's M5 checksum through `$REF_PATH`/`$REF_CACHE` and could
decode against a *different* build of the same genome, yielding plausible but incorrect
sequence.

**What gets converted.** `samtools fasta -F 0x900 -N` — primary records only, with the
`/1`,`/2` mate suffix forced on. That is exactly one full-length copy of every read:
each read has one primary record carrying full-length SEQ, whereas supplementary records
are hard-clipped slices of reads the primary already supplied. Unmapped reads and
duplicates are **kept** — they are genuine distinct reads. Minus-strand records are
reverse-complemented back to sequencing orientation.

**Scratch space.** The KMC backend streams the conversion with no temp file. The HKS
backend cannot — `hks lookup` needs a seekable path — so it materialises a temp FASTA
next to the output first (the system tempdir on cluster nodes is often small or
RAM-backed). **That file is full size**: a 64x human WGS CRAM measured 28.9 GB in and
254 GB out, so budget that space on the output filesystem on top of the BEDs.

## Paired-end reads

`annotate` has no concept of pairing: each mate is an independent query sequence.
What matters is whether the output can still be joined back into fragments afterwards.

Read-level output is identified by **ordinal rank**, and `--query-names` is *refused*
for read-level input. This is not a style preference — `hks` collects every sequence
name into memory before querying (`load_seq_names` builds a `Vec<String>` of the lot),
which is fine for an assembly's few thousand contigs and fatal for reads. Measured on a
64x human WGS CRAM: 1.315 G names at ~68 bytes each is ~90 GB on top of the ~12.3 GB
index, and the run is OOM-killed roughly 20 minutes in, having already spent 15 of them
decoding. Raising the memory limit buys nothing but longer identifiers.

Nothing is lost. Rank N is the Nth record of the query file, so the mapping back to read
names is a property of the input, and `--query-names-sidecar` records it next to the output
(see [Query names](#query-names)). Mates group by stripping the `/1`,`/2` suffix — which is
why `-N` is forced during the alignment decode. In a coordinate-sorted BAM/CRAM the two
mates of a fragment are nowhere near each other in rank, so pairing needs the names; in
split R1/R2 FASTQ files mate *i* of one file is mate *i* of the other, and the ranks pair
by themselves.

Note the coordinates are k-mer offsets **within each read**, not genomic positions: a
151 bp read queried at k=31 spans 0..121.

## Query names

`--query-names-sidecar` writes `<outdir>/<input>.query_names.txt.gz` for each input: one
record name per line, in the order the query assigns ranks, so line N+1 is rank N. The name
is derived exactly as `hks` derives it: drop the record's first byte, then take the first
whitespace-delimited token (`SRR123.7` from `@SRR123.7 1 length=151`; `read1` from `@ read1`;
a CRLF header loses its `\r`), with the `/1`,`/2` mate suffix for BAM/CRAM input. Like every `.gz` KaryoScope writes it is bgzipped, so any gzip reader opens it. It is
named from the input alone — not the database — because the mapping does not depend on
what the reads were annotated against, so one CRAM run against two databases yields one
sidecar.

What it costs depends on where the names come from:

| Input | Backend | How the names are produced |
|---|---|---|
| BAM / CRAM | HKS | Teed off the `samtools fasta` decode `annotate` runs anyway, on the way to the temp FASTA. Nearly free; regenerating the mapping afterwards would be a second full decode (~25 min on a 56 GB CRAM). |
| BAM / CRAM | KMC | Teed off the same decode as it streams into `get_featureIDs`, via a FIFO. Nearly free. (If a rerun reuses the combined BED from an earlier run, no decode happens and the sidecar is produced by a decode of its own.) |
| FASTA / FASTQ | HKS | One streaming read of the file (`gzip -dc` in front of it if gzipped) selecting header lines with `sed`, then `bgzip -@ threads`. Measured at ~650 MB/s of FASTQ on a laptop. Small next to `hks`, which reads the input once per feature set. |
| FASTA / FASTQ | KMC | The same pass. `get_featureIDs` reads the input once for all feature sets, so this is one extra read of equal size. |

Only `hks` emits ranks; `get_featureIDs` (the KMC backend) always writes names into the BED, so there the sidecar is simply the list of record names in file order, written for consistency rather than out of necessity, and it lists every record whether or not `get_featureIDs` reported anything for it.

The scan mirrors how `hks` reads the file. Gzip is detected from the magic bytes and FASTA from FASTQ by the first byte, not the filename, so an extensionless FASTQ gets a correct sidecar. FASTQ headers are taken every fourth line without parsing, which is exactly the shape `hks` parses: a record spanning more lines is a parse error there, not a rank. Headers are picked by position rather than by a leading `@`, since a quality line may begin with `@`. sed rather than awk because on the BSD awk of macOS the awk form was ten times slower than the whole hks lookup on the same file.

Without the flag the mapping is reproducible by hand, since the query file is a
deterministic function of the input:

```bash
# BAM / CRAM: the same decode annotate performs
samtools fasta -F 0x900 -N --reference GRCh38.fasta tumor.cram | grep '^>'
# FASTA
sed -n '/^>/p' asm.fa | sed 's/^.//;s/^[[:space:]]*//;s/[[:space:]].*//'
# FASTQ (four lines per record, as hks reads it)
gzip -dc reads.fastq.gz | sed -n 'p;n;n;n' | sed 's/^.//;s/^[[:space:]]*//;s/[[:space:]].*//'
```

## Output

With `--format kst`, see [Compact output](#compact-output---format-kst). Otherwise, for each feature set, `annotate` writes up to two BEDs:

- Presmoothed (raw): `<input>.<dbid>.<feature_set>.presmoothed.bed[.gz]`
- Smoothed (hierarchy-smoothed): `<input>.<dbid>.<feature_set>.smoothed.bed[.gz]`

With `--query-names-sidecar`, one more file per input, independent of the database: `<input>.query_names.txt.gz` (see [Query names](#query-names)).

Each BED's 4th column is the human-readable feature name. k-mers absent from the index render as `novel`. The `.gz` suffix is present unless `--no-bgzip` is passed; `--no-bgzip` keeps BEDs as plain text.

Smoothing promotes short noisy intervals (especially short `novel` runs flanked by specific features) to the lowest common ancestor of their flankers in the feature set's hierarchy.

For human-scale inputs, use at least 16 threads. Memory to request depends on the backend and input shape. With the **HKS backend**, the peak is set by the index rather than the input — the shared base index plus one feature set's labeling at a time — so it is ~10 GB whether you annotate a single haplotype or a combined diploid assembly such as HG002 v1.1 (request ≥ 16 GB). The **KMC backend**'s largest process peaks at ~36 GB on a haplotype and ~46 GB on a combined diploid, and its smoothing workers can collectively need as much again (request ≥ 96 GB). At `-t 16`, a haplotype takes ~7–9 min (HKS) or ~13–14 min (KMC); HG002 v1.1 combined takes ~14–16 min (HKS) or ~23–27 min (KMC).

## Compact output: `--format kst`

BED stays the default, because most downstream processing wants it. `--format kst` is for groups annotating many samples who need to minimise disk space. It writes the presmoothed tracks as [karyozip](https://github.com/AndreaGuarracino/karyozip) `.kst` files and nothing else:

- **One merged `<input>.<dbid>.kst`** when every feature set of the database is annotated and the database has a `features.tsv` (e.g. `HKS_human_CHM13_v2`, and every KMC database). It holds all sets as `features.tsv` featureIDs.
- **One `<input>.<dbid>.<set>.kst` per set** otherwise: a `--feature-set` subset, or a database without `features.tsv` (e.g. `HKS_hpv_PaVE`).
- Each `.kst` is indexed (`<file>.kzi`), so single-read queries stay fast.
- `<input>.<dbid>.annotate_manifest.json` records the karyozip version, the effort, and each output's size and SHA-256.

No smoothing pass and no bgzip pass run: the smoothed track is not stored, because karyozip rebuilds it on read, row for row as `annotate` would have written it. `--smoothed-bed` writes the smoothed BEDs as well, for tools that want them on disk. `--no-keep-presmoothed` does not apply.

Sizes for all six sets of a 2,356-read HG002 sample (64 Mbp) against `HKS_human_CHM13_v2`:

| Storage | Size | Smaller than BED by |
|---|---|---|
| 6 presmoothed + 6 smoothed `.bed.gz` (`--format bed`) | 70.2 MB | 1× |
| 6 per-set `.kst`, `-e low` | 5.95 MB | 12× |
| 1 merged `.kst`, `--kst-effort low` (2.0 s to write) | 3.35 MB | 21× |
| 1 merged `.kst`, `--kst-effort medium` (17.7 s to write) | 1.61 MB | 44× |

**Peak disk** depends on the installed `hks`:

- With an `hks` whose `lookup` has `--format kst` (and karyozip ≥ 0.3.0), on an HKS database, `hks` writes each set's `.kst` itself and **no presmoothed BED is written at all**. Per-set tracks land under their final names, indexed by `hks`; for a merged track `hks` writes per-set `.kst` into `<input>.<dbid>.kst-sets.tmp/`, which `karyozip compress -d` merges (to the same bytes as from BEDs) before the directory is removed. The peak is then a few percent of a BED run's.
- Otherwise (an older `hks`, or the KMC backend) the backend writes each set's presmoothed BED and karyozip compresses it, so the peak is that of the BEDs (a merged track needs every set's at once; per-set tracks are compressed as each set finishes).

`annotate` picks the first when it can, falls back to the second with a log line, and sizes its disk check for whichever it runs. The manifest's `kst.written_by` says which one wrote the track.

**Reading a `.kst` back.** Every KaryoScope command that reads annotate output (`bin`, `remap-bed`, `scaffold`, `centromeres`, `karyotype`) accepts the `.kst` in place of the BEDs: the pipeline commands find `<input>.<dbid>.kst` next to the input as they would the BEDs, and `bin`/`remap-bed` take it with `--kst-track smoothed|presmoothed` and `--feature-set`/`--db`. Outside KaryoScope, use `karyozip` directly:

```bash
# A per-set track: the presmoothed BED, or the smoothed one
karyozip decompress -i sample.HKS_hpv_PaVE.type.kst
karyozip decompress -i sample.HKS_hpv_PaVE.type.kst --smooth --hier $DB/hierarchy.tsv

# A merged track: name the set, and give the database's features.tsv and hierarchy.tsv
DB=~/.karyoscope/db/HKS_human_CHM13_v2
karyozip decompress -i sample.HKS_human_CHM13_v2.kst -d $DB/features.tsv --hier $DB/hierarchy.tsv --set chromosome
karyozip decompress -i sample.HKS_human_CHM13_v2.kst -d $DB/features.tsv --hier $DB/hierarchy.tsv --set chromosome --smooth

# Every set as a BED under a directory
karyozip decompress -i sample.HKS_human_CHM13_v2.kst -d $DB/features.tsv --hier $DB/hierarchy.tsv --out-dir beds/

# One read, or a range of one; with --smooth the whole read is smoothed first,
# so a range sees exactly what a full run would
karyozip query -i sample.HKS_human_CHM13_v2.kst -d $DB/features.tsv --hier $DB/hierarchy.tsv \
    --set region --smooth -r 'read_123:1000-5000'
```

Read-level input is identified by rank (see [Query names](#query-names)); karyozip stores consecutive ranks without storing names at all, and names (`--query-names`, assemblies) are stored tokenised. Both come back exactly.

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
- `--format kst` on an HKS database with an `hks` that has `lookup --format kst`: no presmoothed BED is written at all (see [Compact output](#compact-output---format-kst))
- `--format kst` otherwise, with a `--feature-set` subset or a database without `features.tsv`: each set's BED is compressed as soon as it is written, so the peak is one set's presmoothed BED (~0.6 GB per Gbp) rather than all of them
- splitting a diploid assembly into per-haplotype runs

`annotate` estimates this footprint before starting and refuses if `--outdir` can't hold it, so a full disk is caught in a second rather than after the k-mer query. The estimate uses an exact base count from a sibling `.fai` when one exists, and otherwise scales the input file size. Pass `--no-resource-check` if the estimate is wrong for your input.

## See also

- [`karyoscope karyotype`](karyotype.md) — render karyotype SVGs (runs annotate under the hood)
- [`karyoscope bin`](bin.md) — aggregate base-pair BEDs into larger bins
