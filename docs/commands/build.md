# karyoscope build

Build an HKS index database from a genome and per-feature-set BED annotations, then register it so the data commands can use it.

## Synopsis

```
karyoscope build --id ID --sequence GENOME --feature-set NAME=BED [NAME=BED ...] [OPTIONS]
karyoscope build --spec build.yaml [OPTIONS]
```

## Description

`build` produces a complete, registry-ready [HKS](https://github.com/jnalanko/HKS) database — the same layout `download`/`register` expect (`manifest.yaml`, `hierarchy.tsv`, `colors.tsv`, and an `index/` directory) — from a reference genome and one or more feature-set BED files. It runs the HKS `build-base` and `add-feature-set` construction steps for you, gap-fills unannotated regions with a named background feature, derives the label hierarchy and colours, writes the manifest, and records the result in `installed.json`.

The command's contract begins at a **final per-feature-set BED** whose 4th column is the leaf label. Turning raw annotation sources (RepeatMasker, censat, GENCODE, …) into that BED is up to you — see the worked example under [Preparing a feature-set BED](#preparing-a-feature-set-bed).

### Input modes

Each feature set is one of:

- **BED (mode A):** `--feature-set NAME=annot.bed` together with a shared `--sequence` genome. The genome is sliced into one FASTA per leaf label (each region extended by `k-1` bp so every k-mer that starts in the region is captured). Gaps are filled with a background leaf (see below).
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
| `--feature-set NAME=BED` | A feature set as `NAME=annotation.bed` (4th column = leaf label). Repeatable. |
| `--background NAME=LABEL` | Gap-fill label for a feature set (default `background`). Repeatable. |
| `--hierarchy NAME=PATH` | Edge-list (`child parent`) hierarchy for a feature set. Repeatable. |
| `--priority NAME=PATH` | Priority file for a feature set; enables priority mode. Repeatable. |
| `--colors NAME=PATH` | Colours file (`feature<tab>color`) for a feature set. Repeatable. |
| `--flatten` | Pre-flatten overlapping BED regions to one label per base. |
| `--variable-k` | Build a variable-k index (queryable at any k ≤ s, e.g. a k-sweep). Not combinable with `--priority`. |
| `--spec FILE` | Build-spec YAML (alternative to `--id`/`--sequence`/`--feature-set`). |
| `--db-version TEXT` | Database version, semver (default `1.0.0`). |
| `-s`, `--s INTEGER` | Maximum query length / k-mer size (default `31`). |
| `-t`, `--threads INTEGER` | Threads for HKS construction (default `4`). |
| `--mem-gigas INTEGER` | RAM budget (GB) for base-index construction (default `8`). |
| `--external-memory DIRECTORY` | Build the base index in external-memory mode using this scratch dir (lower RAM, slower). |
| `--forward-only` | Do not add reverse-complemented k-mers. |
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
    bed: /path/repeat.bed            # 4th col = leaf label
    background: nonrepeat            # default "background"; null disables gap-fill
    priority: /path/repeat.priority.txt   # optional; enables priority mode
    hierarchy: /path/repeat.hierarchy.txt # optional; else a flat star
    colors: /path/repeat.colors.tsv       # optional; else an auto palette
  - name: gene
    bed: /path/gene.bed
roles: { chromosome_assignment: chromosome }   # optional
smoothing: { recommended_window_bp: 1000 }     # optional
```

Relative paths in the spec are resolved against the spec file's directory.

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

`build` starts from a final BED; producing one is dataset-specific. A worked example for a RepeatMasker-derived repeat feature set (download reference, simplify the RepeatMasker output, optionally priority-merge overlaps, then hand the BED to `build`) is kept with the archived workflow scripts in the KaryoScope-archive repository. The key point is that overlaps and gaps are fine — `build` gap-fills automatically and HKS resolves overlaps per k-mer.

## Notes

- Requires the `hks` binary on `PATH` (or `$KARYOSCOPE_HKS`); `--background` gap-fill also uses `samtools faidx` to index the genome if a `.fai` is not already present.
- Human-scale base-index construction is memory-intensive; use `--mem-gigas` and, if needed, `--external-memory` on a scratch filesystem, and run on a cluster node rather than a login node.

## See also

- [`karyoscope register`](register.md) — register an already-built database
- [`karyoscope info`](info.md) — inspect the built database
- [`karyoscope annotate`](annotate.md) — annotate sequences against it
