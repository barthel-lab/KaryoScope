# karyoscope build

Build an HKS index database from a genome and per-feature-set BED annotations, then register it so the data commands can use it.

## Synopsis

```
karyoscope build --id ID --sequence GENOME --feature-set NAME=BED [NAME=BED ...] [OPTIONS]
karyoscope build --spec build.yaml [OPTIONS]
```

## Description

`build` produces a complete, registry-ready [HKS](https://github.com/jnalanko/HKS) database — the same layout `download`/`register` expect (`manifest.yaml`, `hierarchy.tsv`, `colors.tsv`, and an `index/` directory) — from a reference genome and one or more feature-set BED files. It runs the HKS `build-base` and `add-feature-set` construction steps for you, gap-fills unannotated regions with a named background feature, derives the label hierarchy and colours, writes the manifest, and records the result in `installed.json`.

The command's contract begins at a **final per-feature-set BED** whose 4th column is the **feature label**: a leaf of that feature set's hierarchy. Turning raw annotation sources (RepeatMasker, censat, GENCODE, …) into that BED is up to you — see the worked example under [Preparing a feature-set BED](#preparing-a-feature-set-bed).

## Required inputs

Only two inputs are required — a genome and at least one annotation BED:

| Input | Required | What it is |
| --- | --- | --- |
| **Genome FASTA** (`--sequence`) | **Yes**, for BED (mode A) feature sets | The assembly your annotations are in coordinates of; plain or bgzipped. If a `.fai` is missing, `build` creates one with `samtools faidx`. |
| **Feature-set BED** (`--feature-set NAME=BED`) | **Yes**, at least one | A BED whose **4th column is the feature label** — a leaf of this feature set's hierarchy. With no hierarchy file, each distinct label simply becomes its own leaf. Repeat the flag for more sets. Overlaps are allowed and gaps are filled automatically. |
| Hierarchy (`--hierarchy NAME=PATH`) | No | `child<tab>parent` edge list. Without one, every leaf becomes a child of the root `categorized` (a flat star). |
| Priorities (`--priority NAME=PATH`) | No | Resolves which label wins a k-mer claimed by several. In its 3-column `child priority parent` form it **also supplies the hierarchy**, so one file covers both. |
| Colours (`--colors NAME=PATH`) | No | `feature<tab>color`, or the full `feature_set<tab>feature<tab>color`. Without one, a distinct palette is generated per leaf. |
| Background label (`--background NAME=LABEL`) | No | Names the gap-fill leaf (default `background`). |

Nothing else is needed. In particular there is **no** alignment step, no pre-existing index, and no per-feature FASTA files — `build` slices those out of the genome itself.

The only external requirements are the [`hks`](https://github.com/jnalanko/HKS) binary on `PATH` (or `$KARYOSCOPE_HKS`) and `samtools` for the `.fai`. So the minimum viable build is:

```bash
karyoscope build --id HKS_mygenome --sequence genome.fa.gz --feature-set repeat=repeat.bed
```

For what this costs in time, memory, and disk, see [Resource requirements](#resource-requirements).

### Input modes

Each feature set is one of:

- **BED (mode A):** `--feature-set NAME=annot.bed` together with a shared `--sequence` genome. The genome is sliced into one FASTA per feature label (each region extended by `k-1` bp so every k-mer that starts in the region is captured). Gaps are filled with a background leaf (see below).
- **FASTAs (mode B, spec file only):** a feature set entry with `fastas:` (one file per feature) or `per_seq_file:` (one sequence per feature). No genome slicing and no gap-fill.

### Overlaps, priorities, and variable-k

HKS resolves a k-mer that belongs to several labels via its hierarchy, so **overlapping BED regions are allowed** — you do not need a non-overlapping partition. To control which label wins an overlap, supply a **priority file** for the set (`--priority NAME=file`); HKS then keeps the higher-priority (lower integer) label per k-mer, climbing to the common ancestor only on ties. This is the per-k-mer equivalent of pre-flattening the BED, so it usually replaces it. If you *do* want a hard flattened partition, `--flatten` collapses overlaps to one label per base before indexing.

The priority file accepts either the 3-column `child priority parent` form (which also supplies the hierarchy) or the 2-column `name priority` form. Within any group of siblings, priorities must be either all equal or all distinct (an HKS requirement, checked up front).

By default — no priority file — a set is built as a plain **fixed-k** labeling, which is exactly what `annotate` queries at (k = s). This is the recommended default for a production database.

`--variable-k` (or per set, `variable_k: true` in the spec) instead builds an index that can be queried at **any k ≤ s from a single build** — useful for a k-sweep, e.g. checking that results are stable across k without rebuilding. It works from BED (mode A) input: because HKS variable-k needs a dummy node for each feature sequence's start, the base index is built from the per-feature FASTAs `build` already generates (the genome slices), not from the whole genome. `--variable-k` is mutually exclusive with `--priority`, and it produces a larger index, so use it when you actually need multi-k queries. The manifest's `kmer.type` is set to `variable` when every feature set is variable-k.

### Background (gap-fill) feature

Bases that no feature annotates are, by default, filled with a **background** leaf so they render as a real feature rather than as `novel` (HKS's `none`, meaning "k-mer not in the reference at all"). The background label defaults to `background`; name it per set with `--background NAME=LABEL` (e.g. `--background repeat=nonrepeat`). It is attached under the hierarchy root and coloured grey (`#808080`). Set `background: null` in the spec file to disable gap-fill for a set. Where a BED already tiles the whole genome, no background records are produced.

### Defaults

- No hierarchy for a set → a flat star: every leaf is a child of the root, `categorized`.
- No colours → an auto-generated distinct palette per leaf; grouping/root nodes are `#B0C4DE`, background is `#808080`. Provide your own with `--colors NAME=file` (`feature<tab>color`, or the full `feature_set<tab>feature<tab>color`).
- `features.tsv` is **not** produced: the HKS backend reads label names from the index and never uses it.

## Options

| Option | Description |
| --- | --- |
| `--id TEXT` | Database id (also the directory name). Required unless `--spec`. |
| `--sequence FILE` | Genome FASTA (plain or bgzipped). Required for BED (mode A) feature sets. |
| `--feature-set NAME=BED` | A feature set as `NAME=annotation.bed` (4th column = the feature label, a leaf of its hierarchy). Repeatable. |
| `--background NAME=LABEL` | Gap-fill label for a feature set (default `background`). Repeatable. |
| `--hierarchy NAME=PATH` | Edge-list (`child parent`) hierarchy for a feature set. Repeatable. |
| `--priority NAME=PATH` | Priority file for a feature set; enables priority mode. Repeatable. |
| `--colors NAME=PATH` | Colours file for a feature set: `feature<tab>color`, or the full `feature_set<tab>feature<tab>color`. Repeatable. |
| `--flatten` | Pre-flatten overlapping BED regions to one label per base. |
| `--variable-k` | Build a variable-k index (queryable at any k ≤ s, e.g. a k-sweep). Not combinable with `--priority`. |
| `--spec FILE` | Build-spec YAML (alternative to `--id`/`--sequence`/`--feature-set`). |
| `--db-version TEXT` | Database version, semver (default `1.0.0`). |
| `-s`, `--s INTEGER` | Maximum query length / k-mer size (default `31`). |
| `-t`, `--threads INTEGER` | Threads for HKS construction (default `4`). |
| `--mem-gigas INTEGER` | RAM budget (GB) for base-index construction (default `8`). |
| `--external-memory DIRECTORY` | Build the base index in external-memory mode using this scratch dir (lower RAM, slower). |
| `--forward-only` | Do not add reverse-complemented k-mers. |
| `--exclude TEXT` | Sequence name to exclude from the whole build (e.g. an organelle `ChrM`). Repeatable / comma-separated. See [Excluding sequences](#excluding-sequences). |
| `--db-root DIRECTORY` | Override the database root (default: `$KARYOSCOPE_DB` or `~/.karyoscope/db/`). |
| `--no-register` | Build only; do not record in `installed.json`. |
| `--force` | Overwrite an existing database directory / install record. |
| `--keep-intermediates` | Keep the per-feature FASTA working directory (`<db>/_build_work`). |
| `-h`, `--help` | Show this message and exit. |

## Build-spec file

For databases with several feature sets, a YAML spec is clearer than flags and is checked in alongside the data:

```yaml
id: HKS_mygenome
version: "1.0.0"
sequence: /path/genome.fa.gz
kmer: { s: 31 }
build: { threads: 16, mem_gigas: 8, external_memory: /scratch/tmp }  # last two optional
feature_sets:
  - name: repeat
    bed: /path/repeat.bed            # 4th col = feature label (a hierarchy leaf)
    background: nonrepeat            # default "background"; null disables gap-fill
    priority: /path/repeat.priority.txt   # optional; enables priority mode
    hierarchy: /path/repeat.hierarchy.txt # optional; else a flat star
    colors: /path/repeat.colors.tsv       # optional; else an auto palette
  - name: gene
    bed: /path/gene.bed
exclude: [ChrM, ChrC]                          # optional; sequences to leave out (see below)
roles: { chromosome_assignment: chromosome }   # optional
smoothing: { recommended_window_bp: 1000 }     # optional
```

Relative paths in the spec are resolved against the spec file's directory.

## Excluding sequences

The **`chromosome` feature set declares the karyotype chromosomes** — its leaves are what `karyotype` lays out, and (unless a chromosome is present in the sample) an empty column is drawn for each. So keep **non-karyotype sequences out of the `chromosome` set**: organelles (`ChrM`, `ChrC`, …), unplaced contigs, decoys, and the like.

The simplest way is `--exclude` (or `exclude:` in the spec). Excluded sequences are dropped from **every** feature BED and from the gap-fill index, so no feature set covers them and they read as `none` everywhere — uniform across sets, and absent from the karyotype. They're still real sequence; they're just not karyotype chromosomes. (This is why, for example, human CHM13 databases leave out `chrM`.) Note that `--exclude` filters BEDs and the gap-fill index, not `hierarchy`/`colors`/`priority` files — so also drop the excluded names from those if you list them there.

```bash
karyoscope build --spec build.yaml --exclude ChrM,ChrC
```

## Examples

```bash
# One feature set from a genome + BED, with a named background feature
karyoscope build --id HKS_mygenome --sequence genome.fa.gz \
    --feature-set repeat=repeat.bed --background repeat=nonrepeat

# Priority-resolved overlaps for the repeat set
karyoscope build --id HKS_mygenome --sequence genome.fa.gz \
    --feature-set repeat=repeat.bed --priority repeat=repeat.priority.txt

# Multi-feature-set build from a spec file
karyoscope build --spec build.yaml

# Inspect and use the result
karyoscope info HKS_mygenome
karyoscope annotate -i query.fa.gz --db HKS_mygenome

# Variable-k build + k-sweep from one index (query at several k)
karyoscope build --id HKS_sweep --sequence genome.fa.gz \
    --feature-set gene=gene.bed --variable-k
for k in 21 25 31; do
    karyoscope annotate -i query.fa.gz --db HKS_sweep --k "$k" -o sweep/
done
```

## Preparing a feature-set BED

`build` starts from a final BED; producing one is dataset-specific. A typical recipe is: take a reference annotation source (RepeatMasker, censat, GENCODE, EDTA, a satellite-monomer catalog, …), reduce it to a BED whose 4th column is the feature label, optionally priority-merge overlaps, and hand that BED to `build`. Overlaps and gaps are fine — `build` gap-fills automatically and HKS resolves overlaps per k-mer.

> **Planned: `karyoscope prep-bed`.** Turning those raw annotation sources into a final labelled BED — GFF3/GTF gene models into exon/intron/intergenic, RepeatMasker/EDTA tables into labelled repeat BEDs with a hierarchy, satellite monomer files into merged array bands — is the main friction in building a database today, and each source currently needs its own conversion. A dedicated `karyoscope prep-bed` helper for these common conversions is planned. It will be a **separate subcommand**, not folded into `build`: `build`'s contract stays "a final labelled BED", so it never has to sniff and guess at raw file formats.

## Resource requirements

Measured runs, all with `-t 16` on one cluster node. Wall time is the whole command (genome slicing, base index, and every feature set); peak RSS is the maximum resident set across the process tree.

| Genome | Feature sets | k | Wall | Peak RSS | Database on disk |
| --- | --- | --- | --- | --- | --- |
| *A. thaliana* Col-CEN (135 Mb) | 4 | 31 | 50 s | ~6 GB | 673 MB |
| Human hg38, cytoband only (3.1 Gb) | 1 | 31 | 7 m | ~81 GB | 12.2 GB |
| Human CHM13 v2 (3.1 Gb) | 6 | 11 | 8 m | ~76 GB | 19 MB |
| Human CHM13 v2 | 6 | 21 | 18 m | ~73 GB | 19 GB |
| Human CHM13 v2 | 6 | 31 | 21 m | ~82 GB | 21 GB |
| Human CHM13 v2 | 6 | 41 | 29 m | ~147 GB | 23 GB |
| Human CHM13 v2 | 6 | 51 | 31 m | ~170 GB | 23 GB |
| Human CHM13 v2 | 6 | 61 | 30 m | ~151 GB | 24 GB |

> **Treat the memory figures as approximate and size jobs with headroom.** They
> are the larger of two measurements that disagree by up to 15% on the same run
> (`/usr/bin/time -v`, which reports each process's own peak, versus Slurm's
> sampled `MaxRSS`); neither is a guaranteed upper bound. The variation between
> k=41, 51 and 61 is measurement noise, not a real effect of k — see below.
> **Request ~1.5x the figure here**, e.g. 128 GB for a human build at k≤31.

Three things follow, and the first two are easy to get wrong:

### `--mem-gigas` is not a memory cap

It is the budget handed to SBWT construction, not a ceiling on the process. Every human-scale row above used `mem_gigas: 32` and still peaked far higher. **Do not size a job from `--mem-gigas`** — size it from this table, or use `--external-memory` (below).

### Peak memory doubles above k=32, and is otherwise flat

Peak RSS tracks *number of distinct k-mers × bytes per k-mer*, and has almost nothing to do with the size of the database produced — the k=11 build peaked at 72 GB while emitting a 19 MB index. The bytes-per-k-mer term is a step function: SBWT packs a k-mer into `⌈k/32⌉` 64-bit words, **padded**, so cost jumps at k=32, 64, 96, 128 and is flat between. That is why k=11/21/31 all sit near 70 GB and k=41/51/61 all sit near 148 GB — a k=41 k-mer occupies the same 16 bytes as a k=64 one.

Measured across the six builds above: roughly 70–82 GB for every k ≤ 31, and roughly 147–170 GB for every k ≥ 41 — flat within each band, with the spread inside a band being measurement noise rather than an effect of k.

Practical consequence: **k ≤ 31 costs about half the RAM of k > 32**, and within a band, raising k is nearly free in memory.

### Low on RAM? Use `--external-memory`, not fewer threads

Lowering `--threads` does **not** reduce peak memory — it only makes the build slower. `--external-memory DIR` is the lever, switching to the disk-based construction algorithm. Both measured on the CHM13 v2 six-set build at k=31:

| | Peak RSS | Wall |
| --- | --- | --- |
| `-t 16` (default algorithm) | 78.0 GB | 23 m 44 s |
| `-t 4` | 78.0 GB | 39 m 39 s |
| `-t 1` | 74.8 GB | 1 h 40 m |
| `-t 16 --external-memory` | **13.8 GB** | 32 m 27 s |

Dropping from 16 threads to 1 costs **4.2x the wall time to save 4% of the
memory** — so threads are not the lever, even at the extreme.

`--external-memory` gives a **5.7× smaller memory footprint for ~1.4× the wall time**, and produces the same index. This is what makes a human-scale build feasible on a workstation rather than a cluster node. Point `DIR` at a filesystem with room for the intermediates.

### Disk

Budget roughly **2× the final database size** during the build. `build` slices the genome into one FASTA per feature label before indexing; with feature sets that tile the genome, that working directory is about one genome copy per feature set (~18 GB for CHM13 v2's six sets). It is removed at the end unless you pass `--keep-intermediates`.

## Notes

- Requires the `hks` binary on `PATH` (or `$KARYOSCOPE_HKS`); `--background` gap-fill also uses `samtools faidx` to index the genome if a `.fai` is not already present.
- Human-scale base-index construction is memory-intensive; see [Resource requirements](#resource-requirements) before sizing a job, and run on a cluster node rather than a login node (or use `--external-memory`).

## See also

- [`karyoscope register`](register.md) — register an already-built database
- [`karyoscope info`](info.md) — inspect the built database
- [`karyoscope annotate`](annotate.md) — annotate sequences against it
