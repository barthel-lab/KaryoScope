# karyoscope remap-bed

Apply an existing scaffold map to a separately-annotated BED, alignment-free.

## Synopsis

```
karyoscope remap-bed -b BED -m MAP -o OUTPUT [OPTIONS]
```

## Description

`karyoscope scaffold` builds a `scaffold_map.tsv` and, in the same run, rewrites the BEDs it annotated into scaffolded coordinates. `karyoscope remap-bed` is the *standalone* version: it takes a BED that was annotated separately — possibly against a **different database** than the one used to derive the map — and rewrites it into the scaffolded coordinate system using a previously-built map.

The motivating case is a two-database workflow. Lay out an assembly with a roles-bearing database (e.g. `KS_human_CHM13_v2`), annotate a feature set from a plot-only database (e.g. a cytoband database with no `chromosome`/`region` roles) on the *original* FASTA, then `remap-bed` that annotation into scaffold coordinates so it lines up with the region BEDs. It is the file-producing counterpart of the in-memory remap that `karyotype --scaffold-db` performs while rendering.

The BED carries no pointer back to its source FASTA, so before rewriting, `remap-bed` validates that the BED and map plausibly describe the same assembly using what the map records (`original_name` and per-contig `length`):

- **Hard error** when no contig name in the BED appears in the map — the two files describe different assemblies (or an argument was swapped).
- **Hard error** when a contig present in both has a BED interval whose end exceeds the contig length recorded in the map — the BED was annotated against a different (longer) sequence, which would corrupt the flipped-contig mirror math.
- **Warning** when the BED filename's stem does not match the map's source FASTA stem — advisory only, since files are routinely renamed. `--strict` promotes this (and any map contig with no records in the BED) to a hard error.

A BED contig that is *absent* from the map is not an error: the map only lists contigs that survived scaffolding (length filter + leaf chromosome), while the original-coordinate BED contains every contig. Those contigs are dropped from the output, and the count is reported.

## Options

| Option | Description |
| --- | --- |
| `-b, --bed FILE` | Annotation BED in original (unscaffolded) contig coordinates. `.gz` supported. [required] |
| `-m, --map FILE` | A `scaffold_map.tsv` written by `karyoscope scaffold` (or `karyotype`). [required] |
| `-o, --output FILE` | Output BED in scaffolded coordinates. Gzipped iff the path ends in `.gz`. [required] |
| `--strict` | Promote advisory checks (filename-stem mismatch; map contigs with no records in the BED) from warnings to hard errors. |

## Example

Two-database flow — layout from `KS_human_CHM13_v2`, plot a cytoband feature set:

```
karyoscope scaffold -i hap1=hap1.fa.gz --db KS_human_CHM13_v2 --mode bed
karyoscope annotate -i hap1.fa.gz --db KS_human_CHM13_cytoband
karyoscope remap-bed \
    -b hap1.KS_human_CHM13_cytoband.cytoband.smoothed.bed.gz \
    -m hap1.KS_human_CHM13_v2.scaffold_map.tsv \
    -o hap1.cytoband.scaffolded.bed.gz
```

## See also

- [`scaffold`](scaffold.md) — build the `scaffold_map.tsv` (and rewrite the layout database's own BEDs).
- [`annotate`](annotate.md) — produce the per-feature-set BED to remap.
- [`karyotype`](karyotype.md) — `--scaffold-db` applies the same two-database remap in memory while rendering SVGs.
