# karyoscope centromeres

Extract per-contig centromere coordinates from a genome assembly.

## Synopsis

```
karyoscope centromeres -i [NAME=]PATH [OPTIONS]
```

## Description

For each scaffolded contig, this command identifies the centromere's start and end from the binned scaffolded region BED. A coarse pass (1 Mb bins by default) finds the minimum start and maximum stop of bins classified as centromeric, and an optional fine pass (100 kb bins by default) tightens the call within a window around the coarse range. Like `scaffold`, it auto-derives its prerequisites (the scaffold cascade and bin); passing `--no-auto` turns missing inputs into hard errors. Which feature set drives detection comes from the manifest role `centromere_detection`, with a fallback chain to `region_assignment` and then the literal `region`. Coordinates are reported in the scaffolded (post-flip) coordinate system, so `karyoscope karyotype --mode centromere` can consume them directly. Read-level inputs (FASTQ/BAM) are rejected.

## Options

| Option | Description |
| --- | --- |
| `-i, --input TEXT` | FASTA-format genome assembly. Repeat per haplotype. Form: `NAME=PATH` or bare `PATH` to auto-infer the label from the filename stem. Read-level inputs (FASTQ / BAM) are rejected. [required] |
| `--telo TEXT` | Optional precomputed seqtk telo output. Form: `NAME=PATH`. Without this, the cascade runs seqtk telo internally. |
| `--split-haps TEXT` | Optional regex applied per contig name; capture group 1 is the hap label. |
| `--db TEXT` | Database id. Default: the unique installed database if exactly one is installed. |
| `--db-root DIRECTORY` | Override the database root directory (default: `$KARYOSCOPE_DB` or `~/.karyoscope/db/`). |
| `--centromere-feature-set TEXT` | Override which feature set drives centromere detection. Default: `manifest.roles.centromere_detection`, with chain fallback to `roles.region_assignment` and then the literal `region`. |
| `--coarse-bin-size INTEGER` | Bin size (bp) for the coarse pass. 1 Mb matches the manuscript benchmarks. [default: 1000000] |
| `--fine-bin-size INTEGER` | Bin size (bp) for the optional fine-refinement pass. Pass 0 to disable refinement. [default: 100000] |
| `--min-scaffold-length INTEGER` | Drop contigs shorter than this (no telomere) during the scaffold step. [default: 5000000] |
| `--acrocentric TEXT` | Chromosome name to treat as acrocentric during scaffold's flip decision. Repeatable; accepts comma-separated lists. Default: human acrocentrics with a warning. |
| `-t, --threads INTEGER` | Threads for auto-run annotate invocations. [default: 0] |
| `--bgzip / --no-bgzip` | bgzip the output centromere BED. [default: bgzip] |
| `--auto / --no-auto` | Auto-derive missing inputs (scaffold cascade, bin). Disable to require everything upfront. [default: auto] |
| `-o, --outdir DIRECTORY` | Where to write the centromere BEDs. Default: same directory as each `--input`. |
| `-h, --help` | Show this message and exit. |

## Examples

```bash
# Extract centromere coordinates for a diploid assembly
karyoscope centromeres -i hap1=hap1.fa.gz -i hap2=hap2.fa.gz --db KS_human_CHM13_v2 -o results/

# Coarse pass only (disable fine refinement)
karyoscope centromeres -i asm.fa --fine-bin-size 0 -o results/
```

## Output

Per input, the command writes `<input_stem>.<dbid>.centromeres.bed[.gz]`. This is a 3-column BED (contig, start, end) — it is coordinate-only with no feature column, which is why the filename carries no feature-set segment. Contigs with no centromeric content are omitted. The file is `.gz` unless `--no-bgzip` is passed.

## See also

- [`karyoscope karyotype`](karyotype.md) — consumes centromere coordinates for `--mode centromere`
- [`karyoscope scaffold`](scaffold.md) — the scaffolding step centromere detection builds on
