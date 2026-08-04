# karyoscope build

Build an HKS index database from a genome and per-feature-set BED annotations, then register it so the data commands can use it.

## Synopsis

```
karyoscope build --id ID --sequence GENOME --feature-set NAME=BED [NAME=BED ...] [OPTIONS]
karyoscope build --spec build.yaml [OPTIONS]
```

## Description

`build` produces a complete, registry-ready [HKS](https://github.com/jnalanko/HKS) database — the same layout `download`/`register` expect (`manifest.yaml`, `hierarchy.tsv`, `colors.tsv`, and an `index/` directory) — from a reference genome and one or more feature-set BED files. It runs the HKS `build-base` and `add-feature-set` construction steps for you, gap-fills unannotated regions with a named background feature, derives the label hierarchy and colours, writes the manifest, and records the result in `installed.json`.

The command's contract begins at a **final per-feature-set BED** whose 4th column is the **feature label**: a leaf of that feature set's hierarchy. Turning raw annotation sources (RepeatMasker, censat, GENCODE, …) into that BED is [`karyoscope prep-bed`](prep-bed.md)'s job — see [Preparing a feature-set BED](#preparing-a-feature-set-bed).

## Required inputs

Only two inputs are required — a genome and at least one annotation BED:

| Input | Required | What it is |
| --- | --- | --- |
| **Genome FASTA** (`--sequence`) | **Yes** | The assembly your annotations are in coordinates of; plain or bgzipped. If a `.fai` is missing, `build` creates one with `samtools faidx`. |
| **Feature-set BED** (`--feature-set NAME=BED`) | **Yes**, at least one | A BED whose **4th column is the feature label** — a leaf of this feature set's hierarchy. With no hierarchy file, each distinct label simply becomes its own leaf. Repeat the flag for more sets. Overlaps are allowed and gaps are filled automatically. |
| Hierarchy (`--hierarchy NAME=PATH`) | No | `child<tab>parent` edge list. Without one, every leaf becomes a child of the root `categorized` (a flat star). See [Hierarchy](#hierarchy). |
| Priorities (`--priority NAME=PATH`) | No | Resolves which label wins a k-mer claimed by several. In its 3-column `child priority parent` form it **also supplies the hierarchy**, so one file covers both. See [Priorities](#priorities). |
| Colours (`--colors NAME=PATH`) | No | `feature<tab>color`, or the full `feature_set<tab>feature<tab>color`, with an optional 4th `legend_group` column. Without one, a distinct palette is generated per leaf. See [Colours](#colours). |
| Background label (`--background NAME=LABEL`) | No | Names the gap-fill leaf (default `background`). |

Nothing else is needed. In particular there is **no** alignment step, no pre-existing index, and no per-feature FASTA files — `build` slices those out of the genome itself.

The only external requirements are the [`hks`](https://github.com/jnalanko/HKS) binary on `PATH` (or `$KARYOSCOPE_HKS`) and `samtools` for the `.fai`. So the minimum viable build is:

```bash
karyoscope build --id HKS_mygenome --sequence genome.fa.gz --feature-set repeat=repeat.bed
```

For what this costs in time, memory, and disk, see [Resource requirements](#resource-requirements).

### Hierarchy

An edge list, `child<tab>parent`, one edge per line, rooted at `categorized`.
Every label appearing in the BED must be a node in it. Without a hierarchy file
each distinct BED label becomes a child of the root — a flat star, under which
every shared k-mer resolves to `categorized`.

The hierarchy also drives smoothing and rendering: `bin` and `karyotype` prefer
leaf labels, so anything you want drawn must be a leaf rather than an interior
node.

### How labels are resolved

HKS indexes **k-mers**, not coordinates. A feature set's BED is sliced into one
FASTA per label, and a k-mer is claimed by every label whose sequence contains
it. Two consequences follow:

- Overlapping BED regions are fine; any set of labelled regions works.
- The same k-mer sequence occurring in two different regions is claimed by both
  **whether or not those regions overlap in coordinates** — including regions on
  different chromosomes.

Where a k-mer has more than one claimant, the label is the claimants' **lowest
common ancestor** in the hierarchy, unless priorities say otherwise.

Take this hierarchy and BED:

```
# repeat.hierarchy.txt          # repeat.bed
LINE        Transposon          chr1  100  200  LINE
SINE        Transposon          chr1  150  250  SINE
Transposon  categorized
```

A k-mer found only in 100–150 is `LINE`; one found only in 200–250 is `SINE`.
One found in both labels' sequence — anywhere — is `Transposon`, the lowest
common ancestor. Under a flat hierarchy that ancestor is the root,
`categorized`; a hierarchy is what lets a shared k-mer resolve to something more
specific.

### Priorities

A priority file overrides the common-ancestor rule: the claimant with the
**lowest priority number wins**, and the ancestor is used only when the top
claimants tie.

```
# repeat.priority.txt
LINE        1  Transposon
SINE        2  Transposon
Transposon  1  categorized
```

With that file, the shared k-mer above is `LINE` rather than `Transposon`.

Both priorities and `--flatten` resolve competing claims when the index is
built, but at different granularity. `--flatten` picks a winner per *base*
before the k-mers are extracted, so the losing label is gone. Priorities keep
every claim and resolve per *k-mer*, which is the finer of the two and does not
depend on the regions overlapping in coordinates.

The file takes either form:

| Form | Columns | Also supplies the hierarchy? |
| --- | --- | --- |
| 3-column | `child priority parent` | Yes — one file covers both |
| 2-column | `name priority` | No — pass `--hierarchy` separately |

**Within any group of siblings, priorities must be either all equal or all
distinct.** This is an HKS requirement (a mix makes its priority-aware ancestor
operation non-associative) and is checked before the build starts:

```
LINE  1  Transposon      LINE  1  Transposon      LINE  1  Transposon
SINE  1  Transposon      SINE  2  Transposon      SINE  1  Transposon
LTR   1  Transposon      LTR   3  Transposon      LTR   2  Transposon
      valid                    valid                    rejected
```

Interior nodes need priorities only to satisfy that rule; they are never
assigned to a k-mer directly.

### Colours

`feature<tab>color`, or the full `feature_set<tab>feature<tab>color`, with an
optional 4th `legend_group` column. Without a colours file, `build` generates a
distinct palette per leaf; interior nodes are `#B0C4DE` and the background leaf
is `#808080`.

Rows are written to the database's `colors.tsv` in hierarchy order — each root,
then every child in edge order.

#### Grouping the legend

A feature set with hundreds of leaves in a handful of colours produces a legend
that dwarfs the figure. Features sharing a `legend_group` collapse to **one
legend row**, labelled by the group:

```
feature_set	feature	color	legend_group
cytoband	p11.1	#000000	gpos100
cytoband	q21.3	#000000	gpos100
cytoband	p13.2	#ffffff	gneg
```

A blank cell means "ungrouped" — that feature keeps its own row. For the CHM13
cytoband database this turns 833 legend entries into 9, labelled by Giemsa
stain. Groups appear in the legend in order of first appearance in `colors.tsv`.

`build` carries the column through to the database's `colors.tsv`, so the file
you supply and the file the database ships have the same shape. The column is
written only if at least one feature declares a group, so a build that doesn't
ask for grouping produces the same 3-column file as before.

A legend row is one swatch and one label, so **every feature in a group must
share a colour**; `build` fails if a group spans two. Two groups may share a
colour, which is legal and merely redundant.

### Background (gap-fill) feature

Bases that no feature annotates are, by default, filled with a **background**
leaf so they render as a real feature rather than as `novel` (HKS's `none`,
meaning "k-mer not in the reference at all"). The background label defaults to
`background`; name it per set with `--background NAME=LABEL` (e.g.
`--background repeat=nonrepeat`). It is attached under the hierarchy root and
coloured grey (`#808080`). Set `background: null` in the spec file to disable
gap-fill for a set. Where a BED already tiles the genome, no background records
are produced.

### Fixed-k and variable-k

By default a database is **fixed-k**: the index answers exactly one k-mer
length, the `s` it was built with, and that is what `annotate` queries at. This
is the recommended default for a production database.

`--variable-k` (per set, `variable_k: true`) instead builds an index that
answers **any k ≤ s from a single build**, for a k-sweep. Three things to know
before using it:

- **It cannot be combined with priorities.** HKS rejects the combination, so a
  priority-resolved database is necessarily fixed-k. `build` checks this up
  front.
- **The base index is built differently.** Variable-k needs a dummy node per
  feature sequence start, so the base is built from the per-feature FASTAs
  rather than from the genome. With feature sets that tile the genome, its input
  is roughly one genome copy *per set* — six tiling sets over a 3.1 Gb genome
  means ~19 Gb of input rather than 3.1 Gb.
- **The index is larger.**

In the manifest, `kmer.type` is `variable` only when every feature set is
variable-k. `kmer.max_size` is the largest queryable k, which on a fixed-k
index equals `kmer.size` — it is **not** a range you may query within. On a
fixed-k database `annotate --k` accepts only `kmer.size` and rejects anything
else.

### Defaults

- No hierarchy for a set → a flat star: every leaf is a child of the root, `categorized`.
- No colours → an auto-generated distinct palette per leaf.

## Options

| Option | Description |
| --- | --- |
| `--id TEXT` | Database id (also the directory name). Required unless `--spec`. |
| `--sequence FILE` | Genome FASTA (plain or bgzipped). Required for BED feature sets. |
| `--feature-set NAME=BED` | A feature set as `NAME=annotation.bed` (4th column = the feature label, a leaf of its hierarchy). Repeatable. |
| `--background NAME=LABEL` | Gap-fill label for a feature set (default `background`). Repeatable. |
| `--hierarchy NAME=PATH` | Edge-list (`child parent`) hierarchy for a feature set ([Hierarchy](#hierarchy)). Repeatable. |
| `--priority NAME=PATH` | Priority file for a feature set ([Priorities](#priorities)); enables priority mode. Repeatable. |
| `--colors NAME=PATH` | Colours file for a feature set: `feature<tab>color`, or the full `feature_set<tab>feature<tab>color`, optionally with a 4th `legend_group` column ([Colours](#colours)). Repeatable. |
| `--flatten` | Pre-flatten overlapping BED regions to one label per base ([Priorities](#priorities)). |
| `--variable-k` | Build a variable-k index, queryable at any k ≤ s ([Fixed-k and variable-k](#fixed-k-and-variable-k)). Not combinable with `--priority`. |
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

### Feature sets from FASTAs instead of a BED

A feature set can name its sequences directly rather than being sliced out of a
genome. This form is spec-only, and replaces `bed:` with one of:

```yaml
feature_sets:
  - name: marker
    fastas: [/path/markerA.fa, /path/markerB.fa]   # one file per feature
  - name: probe
    per_seq_file: /path/probes.fa                  # one sequence per feature
```

There is no coordinate system here, so `background:` and `flatten:` do not
apply and are rejected. `--sequence` is needed only if some *other*
feature set uses `bed:`.

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

`build` starts from a final BED whose 4th column is the feature label. Overlaps and gaps are fine — `build` gap-fills automatically and HKS resolves overlaps per k-mer.

Use [`karyoscope prep-bed`](prep-bed.md) to produce that BED from the formats annotation usually arrives in: RepeatMasker and EDTA tables, GFF3/GTF gene models, UCSC cytoband tables, satellite monomer catalogs, and a plain `.fai`. It writes the BED and hierarchy and prints the `feature_sets:` stanza to paste (or append) here.

For complete worked examples, see the [database recipes](../recipes/).

`prep-bed` is a **separate subcommand** rather than a mode of `build`. That keeps `build`'s input contract to one thing — a final labelled BED — so it reads a single known format instead of identifying whichever the source annotation happened to use. It also keeps gap-filling, flattening and sequence exclusion in one place: `background:`, `flatten:` and `exclude:` below.

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

Three things follow:

### `--mem-gigas` sets the SBWT construction budget

It is the allowance handed to the SBWT construction step, and total process memory runs well above it: every human-scale row above used `mem_gigas: 32` and peaked far higher. **Size a job from the table above**, or from `--external-memory` (below).

### Peak memory doubles above k=32, and is otherwise flat

Peak RSS tracks *number of distinct k-mers × bytes per k-mer*; the size of the database produced is unrelated — the k=11 build peaked at 72 GB while emitting a 19 MB index. The bytes-per-k-mer term is a step function: SBWT packs a k-mer into `⌈k/32⌉` 64-bit words, **padded**, so cost jumps at k=32, 64, 96, 128 and is flat between. That is why k=11/21/31 all sit near 70 GB and k=41/51/61 all sit near 148 GB — a k=41 k-mer occupies the same 16 bytes as a k=64 one.

Measured across the six builds above: roughly 70–82 GB for every k ≤ 31, and roughly 147–170 GB for every k ≥ 41 — flat within each band, with the spread inside a band being measurement noise rather than an effect of k.

Practical consequence: **k ≤ 31 costs about half the RAM of k > 32**, and within a band, raising k is nearly free in memory.

### `--external-memory` is the lever on a memory-limited machine

`--external-memory DIR` switches base-index construction to a disk-based algorithm: a **5.7× smaller memory footprint for ~1.4× the wall time**, producing the same index. It is what makes a human-scale build feasible on a workstation rather than a cluster node. Point `DIR` at a filesystem with room for the intermediates.

`--threads` controls speed rather than memory. Both measured on the CHM13 v2 six-set build at k=31:

| | Peak RSS | Wall |
| --- | --- | --- |
| `-t 16` (default algorithm) | 78.0 GB | 23 m 44 s |
| `-t 4` | 78.0 GB | 39 m 39 s |
| `-t 1` | 74.8 GB | 1 h 40 m |
| `-t 16 --external-memory` | **13.8 GB** | 32 m 27 s |

Dropping from 16 threads to 1 costs **4.2× the wall time and saves 4% of the memory**, so the thread count is worth setting for speed alone.

### Disk

Budget roughly **2× the final database size** during the build. `build` slices the genome into one FASTA per feature label before indexing; with feature sets that tile the genome, that working directory is about one genome copy per feature set (~18 GB for CHM13 v2's six sets). It is removed at the end unless you pass `--keep-intermediates`.

## Notes

- Requires the `hks` binary on `PATH` (or `$KARYOSCOPE_HKS`); `--background` gap-fill also uses `samtools faidx` to index the genome if a `.fai` is not already present.
- Human-scale base-index construction is memory-intensive; see [Resource requirements](#resource-requirements) before sizing a job, and run on a cluster node rather than a login node (or use `--external-memory`).

## See also

- [`karyoscope register`](register.md) — register an already-built database
- [`karyoscope info`](info.md) — inspect the built database
- [`karyoscope annotate`](annotate.md) — annotate sequences against it
