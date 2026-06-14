# karyoscope scaffold

Order, orient, and rename assembly contigs into canonical chromosome scaffolds, alignment-free.

## Synopsis

```
karyoscope scaffold -i [NAME=]PATH [OPTIONS]
```

## Description

`karyoscope scaffold` takes one or more FASTA inputs (typically one per haplotype) and uses each contig's k-mer feature profile to assign it to a chromosome and decide its orientation, then emits scaffolded outputs. Each `-i` input produces its own outputs, preserving topology. The authoritative artifact is `scaffold_map.tsv`, which maps every encoded contig name back to its source. Missing prerequisites (annotate, seqtk telo, bin) are auto-derived; `--no-auto` turns any missing input into a hard error. Which feature sets drive chromosome assignment and orientation comes from the manifest roles (`chromosome_assignment`, `region_assignment`), falling back to the literal `chromosome` and `region` feature sets. Read-level inputs (FASTQ / BAM) are rejected.

## Options

| Option | Description |
| --- | --- |
| `-i, --input TEXT` | FASTA-format genome assembly. Repeat per haplotype. Form: `NAME=PATH` (e.g. `hap1=hap1.fa.gz`) or bare `PATH` to auto-infer the label from the filename stem. Read-level inputs (FASTQ / BAM) are rejected. [required] |
| `--telo TEXT` | Optional precomputed seqtk telo output. Form: `NAME=PATH`. Without this, scaffold runs seqtk telo on each input automatically. |
| `--split-haps TEXT` | Optional regex applied per contig name; capture group 1 is the hap label. Overrides the built-in patterns (hifiasm h[12]tg, hap1/hap2, maternal/paternal). |
| `--db TEXT` | Database id to use. Default: the unique installed database if exactly one is installed. |
| `--db-root DIRECTORY` | Override the database root directory (default: `$KARYOSCOPE_DB` or `~/.karyoscope/db/`). |
| `--feature-set TEXT` | Restrict the scaffolded-BED outputs to this feature set. Repeatable. Default: every feature set declared in the database's manifest. |
| `--bin-size INTEGER` | Bin size (bp) for the orientation BEDs. The manuscript uses 1 Mb. [default: 1000000] |
| `--min-scaffold-length INTEGER` | Drop contigs shorter than this that have no telomere. [default: 5000000] |
| `--acrocentric TEXT` | Chromosome name to treat as acrocentric in the flip decision. Repeatable; accepts comma-separated lists. Default: human acrocentrics (chr13/14/15/21/22) with a warning to set it explicitly for non-human assemblies. |
| `-t, --threads INTEGER` | Threads for auto-run annotate invocations. [default: 0] |
| `--mode [fasta\|bed\|both]` | What to write per input. `fasta` writes a scaffolded FASTA only; `bed` writes per-feature-set scaffolded BEDs only (used by karyotype); `both` writes both. The map and legacy stats files are always written. [default: fasta] |
| `--keep-unscaffolded / --drop-unscaffolded` | In FASTA mode, append contigs that did not get scaffolded (too short, no leaf chromosome) at the end of the output under their original names. Disable to drop them entirely. [default: keep-unscaffolded] |
| `--combine-chromosomes / --no-combine-chromosomes` | Concatenate all contigs of each chromosome+haplotype into a single `<chrom>_<hap>` sequence separated by N gaps. FASTA modes only. Adds a `combined_chromosomes` tag to output filenames and writes an AGP file. In `--mode both`, the per-feature-set BEDs are written in the combined coordinate system. [default: no-combine-chromosomes] |
| `--scaffold-gap-size INTEGER` | Number of N bases to insert between concatenated contigs when `--combine-chromosomes` is set. [default: 100000] |
| `--combine-acrocentrics / --no-combine-acrocentrics` | Also combine acrocentric chromosomes. Off by default. Only has effect with `--combine-chromosomes`. [default: no-combine-acrocentrics] |
| `--bgzip / --no-bgzip` | bgzip the scaffolded output BEDs and FASTA. [default: bgzip] |
| `--auto / --no-auto` | Auto-derive missing inputs (annotate, seqtk telo, bin). Disable to require everything upfront. [default: auto] |
| `-o, --outdir DIRECTORY` | Where to write scaffolded outputs. Default: same directory as each `--input`. |
| `-h, --help` | Show this message and exit. |

## Examples

```bash
# Scaffold a diploid assembly (one FASTA per haplotype) into a scaffolded FASTA
karyoscope scaffold -i hap1=hap1.fa.gz -i hap2=hap2.fa.gz --db KS_human_CHM13_v2 -o results/

# Write per-feature-set scaffolded BEDs as well as the FASTA
karyoscope scaffold -i asm.fa --mode both -o results/

# Concatenate each chromosome's contigs into a single sequence (less fragmented)
karyoscope scaffold -i hap1.fa -i hap2.fa --combine-chromosomes -o results/
```

## Output

- `scaffold_map.tsv` — the authoritative 8-column map (`new_name`, `original_name`, `input_file`, `hap`, `chromosome`, `flipped`, `length`, `stats`).
- `scaffold_stats.tsv` — legacy 2-column stats, kept for back-compat.
- `--mode bed` / `both`: per-feature-set `<input_stem>.<dbid>.<feature_set>.smoothed.scaffolded.bed[.gz]`.
- `--mode fasta` (default) / `both`: `<input_stem>.<dbid>.scaffolded.fa[.gz]`.
- Encoded contig name format: `<chrom>_<hap>_<contig>[_rc]`, where the `_rc` suffix marks a reverse-complemented contig.
- With `--combine-chromosomes`: filenames carry a `combined_chromosomes` tag, and an AGP 2.1 file `<stem>.<dbid>.scaffolded.combined_chromosomes.agp` documents every component placement and gap. Each `(chrom, hap)` group becomes one `<chrom>_<hap>` record; acrocentric groups left uncombined (the default) instead become one record per contig, named in canonical order as `<chrom>_<hap>_<A|B|C...>` (the original contig name and `_rc` suffix are dropped).

`scaffold_map.tsv` is the contract that downstream stages parse; the encoded contig name can change between releases.

## See also

- [`karyoscope karyotype`](karyotype.md) — render karyotypes (drives scaffold under the hood)
- [`karyoscope centromeres`](centromeres.md) — extract centromere coordinates from scaffolded contigs
- [`karyoscope annotate`](annotate.md) — the per-base annotation that feeds scaffolding
