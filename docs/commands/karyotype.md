# karyoscope karyotype

Render karyotype visualizations (SVG, optionally PDF/PNG) from a genome assembly.

## Synopsis

```
karyoscope karyotype -i [NAME=]PATH [OPTIONS]
```

## Description

`karyotype` sits at the top of the pipeline cascade. Given one or more FASTA inputs, it runs the full pipeline through to a karyotype image (annotate -> seqtk telo -> bin -> scaffold -> centromeres for centromere mode -> render). Existing intermediates are reused; passing `--no-auto` turns missing inputs into hard errors. It supports three render modes (genome, subtelomere, centromere) at different bin sizes; `--mode` is repeatable and the default renders every mode. One image is written per (mode, feature_set) combination, where `--feature-set` is repeatable and defaults to every feature set in the manifest. Sex handling via `--sex` and `--sex-determination-system` decides which (chrom, hap) cells render. Read-level inputs (FASTQ/BAM) are rejected.

### Rendering scale and non-human genomes

By default the render **scales to the data**: the longest chromosome fills a fixed height in genome view (the longest centromere in centromere view), the genome-view bin size is derived from the longest sequence (≈ longest / 250), and the scale bar is a "nice" round length for that zoom. So a small genome (e.g. Arabidopsis, ~32 Mb) fills the plot the same way a human one does, and the finer bins preserve feature diversity that a fixed coarse bin would wash out — human output is unchanged from the old fixed scale. Pin a fixed scale with `--pixels-per-mb` (to compare assemblies) or a fixed bin with `--bin-size`.

For non-human assemblies, two things matter: pass `--telo-motif` for the organism's telomere (the seqtk default is the human `TTAGGG`, which finds nothing on e.g. plants), and make sure the database's `chromosome` feature set declares the right karyotype chromosomes with organelles excluded (`build --exclude`) — the layout is seeded from that set, so no `--no-human-chroms` flag is needed.

## Options

| Option | Description |
| --- | --- |
| `-i, --input TEXT` | FASTA-format genome assembly. Repeat per haplotype. Form: `NAME=PATH` or bare `PATH`. Read-level inputs (FASTQ / BAM) are rejected. **[required]** |
| `--telo TEXT` | Optional precomputed seqtk telo output. Form: `NAME=PATH`. |
| `--telo-motif TEXT` | Telomere repeat motif for the auto-run `seqtk telo` (its `-m`). Default: seqtk's `CCCTAA` (the vertebrate `TTAGGG` telomere). Non-vertebrate genomes need their own — e.g. Arabidopsis / plants use `CCCTAAA` (`TTTAGGG`); the human default finds no telomeres on them. |
| `--split-haps TEXT` | Optional regex applied per contig name; capture group 1 is the hap label. |
| `--db TEXT` | Database id. Default: the unique installed database if exactly one is installed. |
| `--db-root DIRECTORY` | Override the database root directory. |
| `--scaffold-db TEXT` | Layout database id. When set, chromosome ordering, region orientation, and centromere detection come from THIS database (which must declare the chromosome/region roles), while the `--feature-set`(s) are plotted and coloured from `--db`. Use to render a plot-only database (e.g. a cytoband database with no chromosome/region sets) by borrowing the layout from a roles-bearing database. Default: unset (`--db` supplies both layout and plotting). |
| `--scaffold-db-root DIRECTORY` | Override the database root directory for `--scaffold-db`. Default: same root as `--db-root`. |
| `--feature-set TEXT` | Feature set to render. Repeatable. Default: every feature set in the manifest (one SVG per set). |
| `--colors FILE` | Custom colour file (same format as a database `colors.tsv`) to use instead of the database default. The colour-file stem is appended to the output filename so it doesn't clash with the default-colour render. Default: the database's `colors.tsv`. |
| `--mode [genome\|subtelomere\|centromere]` | Which view(s) to render. Repeatable. Default: render every mode. |
| `--sex [male\|female\|reference\|unknown]` | Sample sex. `unknown` draws sex-chromosome haps only where data is present. **[default: unknown]** |
| `--sex-determination-system [xy\|x0\|zw\|zo]` | Sex-determination system. **[default: XY]** |
| `--background-color [white\|black]` | Background colour. `white` draws sequence outlines; `black` uses light text. **[default: white]** |
| `--bin-size INTEGER` | Bin size (bp) for the SVG. Default: **data-driven** for genome view (≈ longest sequence / 250, so ~1Mb for human, ~130kb for Arabidopsis — finer bins keep the feature diversity that coarse bins wash out via the plurality-per-bin rule); 100Kb for centromere, 100bp for subtelomere. Only valid with exactly one `--mode`. |
| `--pixels-per-mb FLOAT` | Fix the vertical zoom at this many pixels per Mb (e.g. to compare plots across assemblies at the same scale). Default: **data-driven** so the longest chromosome (genome view) / longest centromere (centromere view) fills a fixed height regardless of genome size. |
| `--subtelomere-boundary INTEGER` | Subtelomere window size (bp). Only used in `--mode subtelomere`. **[default: 250000]** |
| `--min-scaffold-length INTEGER` | Drop contigs shorter than this (no telomere) during the scaffold step. **[default: 5000000]** |
| `--acrocentric TEXT` | Chromosome name to treat as acrocentric during scaffold's flip decision. Repeatable; accepts comma-separated lists. Default: human acrocentrics with a warning. |
| `--combine-chromosomes / --no-combine-chromosomes` | Render combined chromosomes: cascade scaffold with `--combine-chromosomes` so each chromosome+haplotype's contigs are concatenated into a single `<chrom>_<hap>` sequence (separated by N gaps), then lay out the karyotype from those combined BEDs. Tags the output filenames with `combined_chromosomes`. **[default: no-combine-chromosomes]** |
| `--scaffold-gap-size INTEGER` | Number of N bases inserted between concatenated contigs when `--combine-chromosomes` is set. **[default: 100000]** |
| `--combine-acrocentrics / --no-combine-acrocentrics` | Also combine acrocentric chromosomes when `--combine-chromosomes` is set. Off by default: their contigs stay as separate records, each renamed in canonical order to `<chrom>_<hap>_<A\|B\|C...>` (e.g. `chr14_hap2_A`, `chr14_hap2_B`). **[default: no-combine-acrocentrics]** |
| `--no-human-chroms` | Don't seed the layout with the database's declared chromosome set (the `chromosome` feature-set leaves); draw only the chromosomes present in the data. By default a chromosome missing from the sample still gets an empty column. The database's chromosome set defines the karyotype chromosomes, so keep non-karyotype sequences (organelles) out of it — or `build --exclude` them. |
| `--format [svg\|pdf\|png]` | Output format(s). Repeatable. Default: svg only. PDF and PNG are generated by converting the SVG via cairosvg. |
| `--sample-label TEXT` | Sample label rendered in the SVG title band. Default: joined stems of the input FASTAs. |
| `--no-title` | Don't draw the title band at the top of the SVG. |
| `--no-legend` | Don't draw the color legend in the right margin of the SVG. |
| `-t, --threads INTEGER` | Threads for auto-run annotate invocations. **[default: 0]** |
| `--auto / --no-auto` | Auto-derive missing inputs. Disable to require everything upfront. **[default: auto]** |
| `--bgzip / --no-bgzip` | bgzip the intermediate scaffolded BEDs (and the centromeres BED, when `--mode` includes centromere). The SVG/PDF/PNG outputs themselves are unaffected. **[default: bgzip]** |
| `--scaffolding / --no-scaffolding` | Write full-resolution scaffolded BEDs to disk during the cascade. Pass `--no-scaffolding` to skip the expensive per-feature-set rewrite step (saves ~5-10 min on whole-genome HG002): the scaffold_map.tsv is still written and applied at bin time, so the binned-scaffolded BEDs and final SVGs are equivalent. **[default: scaffolding]** |
| `--smoothed / --presmoothed` | Which annotation variant to render. Default `--smoothed` uses the hierarchy-smoothed BEDs; `--presmoothed` uses the raw annotations. **[default: smoothed]** |
| `-o, --outdir DIRECTORY` | Where to write the SVGs. Default: same directory as the first `--input`. |
| `--output FILE` | Explicit output path base. The mode and feature_set will be appended; with `--output foo.svg` you get `foo.<dbid>.<mode>.<fs>.karyotype.svg`. Conflicts with `--outdir` when both are set. |
| `-h, --help` | Show this message and exit. |

## Examples

```bash
# Render all modes x all feature sets for a male diploid sample
karyoscope karyotype -i hap1=hap1.fa.gz -i hap2=hap2.fa.gz --sex male -o results/

# One mode + one feature set (genome view, chromosome set), faster cascade
karyoscope karyotype -i asm.fa --mode genome --feature-set chromosome \
    --sex male --no-scaffolding --no-bgzip -o results/

# Also emit PDF and PNG alongside the SVG
karyoscope karyotype -i asm.fa --mode genome --feature-set chromosome \
    --format svg --format pdf --format png -o results/
```

## Output

One image is written per (mode, feature_set) combination, named `<base>.<dbid>.<mode>.<feature_set>.smoothed.karyotype.svg`. `<base>` is the basename of `--output PATH` when given; otherwise it is derived from the input stems — a single input's stem, or for multi-input runs the longest common prefix of the stems (e.g. `GM04890.haplotype1` + `GM04890.haplotype2` → `GM04890`), falling back to the first stem when the stems share no separator-delimited prefix. Passing `--format pdf` / `--format png` (repeatable) additionally produces those formats by converting the SVG via cairosvg (which needs libcairo at runtime). Using `--presmoothed` renders from raw annotations, and the filename then reflects that variant. The cascade also writes the intermediate scaffolded / binned / centromere BEDs unless suppressed.

## Progress output

The cascade reports as it goes, so a run that takes tens of minutes is distinguishable from a hung one. Nested steps — the `annotate` the cascade runs to derive missing annotations, which is usually the bulk of the wall time — are indented one level, so their headline reads as part of this run rather than a separate command:

```
Rendering karyotypes for hg002v1.1.fasta.gz against HKS_human_CHM13_v2
  3 mode(s) x 6 feature set(s) = 18 render(s)
  Annotating hg002v1.1.fasta.gz against HKS_human_CHM13_v2
    6 feature set(s), 16 thread(s), ~34 GB estimated output
    [1/6] chromosome    4m05s
    ...
  [1/18] genome/chromosome      52.1s
  [2/18] genome/region          48.7s
  ...
Wrote:
  ...
```

Each render's time is measured from the top of its loop, so the first view of a mode carries its binning pass (and any cascade work its feature set needed) rather than reporting a misleadingly fast render.

Pass `-q` (before the subcommand: `karyoscope -q karyotype ...`) to suppress the narration; the closing `Wrote:` block still prints. Use `-v` for the full per-step logging on stderr.

## See also

- [`karyoscope annotate`](annotate.md), [`karyoscope scaffold`](scaffold.md), [`karyoscope centromeres`](centromeres.md) — the cascade stages karyotype drives
