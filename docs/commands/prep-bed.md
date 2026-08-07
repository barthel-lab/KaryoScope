# karyoscope prep-bed

Convert a source annotation into the feature-set BED that [`karyoscope build`](build.md) consumes.

## Synopsis

```
karyoscope prep-bed repeatmasker --input FILE --output BED --hierarchy TSV [OPTIONS]
karyoscope prep-bed edta         --input FILE --output BED --hierarchy TSV [OPTIONS]
karyoscope prep-bed gff-gene     --input FILE --lengths FAI --output BED --hierarchy TSV [OPTIONS]
karyoscope prep-bed cytoband     --input FILE --lengths FAI --output BED --hierarchy TSV [OPTIONS]
karyoscope prep-bed fai          --lengths FAI --output BED [OPTIONS]
karyoscope prep-bed censat       --input FILE --lengths FAI --output BED --hierarchy TSV [OPTIONS]
karyoscope prep-bed asat         --input FILE --output BED --hierarchy TSV [OPTIONS]
karyoscope prep-bed satellite    --input FILE --lengths FAI --output BED --hierarchy TSV [OPTIONS]
karyoscope prep-bed pave         --input FILE --output BED --hierarchy TSV [OPTIONS]
```

## Description

`build` starts from a final labelled BED. `prep-bed` produces one from the annotation formats those BEDs usually come from, so that step is a documented command rather than a bespoke script per dataset.

There is one subcommand **per source format**, not per feature set. Unrelated formats produce the same kind of set — RepeatMasker output and an EDTA GFF3 both yield a `repeat` set but share no parsing at all — so keying on the format lets each subcommand carry exactly the options that apply to it, with no flags that are ignored depending on the input.

| Subcommand | Input | Produces |
| --- | --- | --- |
| `repeatmasker` | RepeatMasker `.out` or the UCSC BED repackaging | `repeat`: one leaf per RepeatMasker class |
| `edta` | EDTA TE GFF3 | `repeat`: one leaf per TE superfamily |
| `gff-gene` | GFF3 or GTF gene models | `gene`: `exon` / `intron` / `intergenic` |
| `cytoband` | UCSC `cytoBand.txt` or `cytoBandMapped` BED | `cytoband`: one leaf per band |
| `fai` | samtools `.fai` | `chromosome`: one leaf per sequence |
| `censat` | CenSat annotation BED | `region`: CenSat features plus `p_arm`/`q_arm` |
| `satellite` | centromeric satellite monomers (GFF or BED) | `region`: satellite bands plus `p_arm`/`q_arm` |
| `pave` | PaVE genome records (JSON) | `gene`: one leaf per papillomavirus ORF, plus the URR |

Each subcommand writes its BED (and, where the format implies one, a hierarchy) and prints the matching build-spec stanza. **The stanza goes to stdout and everything else to stderr**, so you can append it straight to a spec:

```bash
karyoscope prep-bed edta --input EDTA.TEanno.gff3.gz \
    --output repeat.bed --hierarchy repeat.tsv >> build.yaml
```

### What prep-bed does not do

It does not gap-fill, flatten overlaps, or drop sequences. `build` already owns all three — `background:`, `flatten:` and `exclude:` — and doing them twice would mean two places to get them wrong. So a converter emits only the bases its source annotates, and names the gap-fill label in the stanza rather than materialising it.

The exceptions are `censat` and `satellite`, which do tile. That tiling is *semantic*, not a gap-fill: `p_arm` and `q_arm` have to be told apart, and only the annotation knows where the centromere is. `build`'s gap-fill has exactly one label and could not make that distinction.

Sequences a set deliberately does not cover — `_alt`/`_random`/`_fix` scaffolds, `chrUn_*` contigs, organelles — are reported for the spec's top-level `exclude:` rather than given a placeholder label. Older hand-written sets used a literal `exclude` label for this, which predates `build`'s `exclude:` list. A label claims the sequence for the feature set; `exclude:` leaves it uncovered.

### Matching seqids to your assembly

Annotation seqids often differ from the assembly's. Two options, on every subcommand:

- `--seqid-map FILE` — a two-column `old new` table, for accession-style names (`NC_060925.1` → `chr1`).
- `--rename-prefix OLD:NEW` — a leading-prefix rewrite, for systematic differences (`Col-CEN_chr` → `Chr`). Names not starting with `OLD` are left alone.

`--seqid-map` takes precedence where both apply. Converters that take `--lengths` report any annotation seqid absent from the `.fai`, so a mismatch surfaces as a warning rather than an empty output.

### Hierarchies, colours and priorities

`--hierarchy` is written wherever the source format implies one. `--colors` and `--priority` are optional and produced only when you name an output path.

Colours are emitted only where an established KaryoScope palette exists — `repeatmasker`, `gff-gene`, `cytoband`, `asat` and `pave`. For `repeatmasker`, `gff-gene` and `cytoband`, the palette reproduces the shipped reference databases; `asat` colours by suprachromosomal family, and `pave` colours early genes cool, late genes warm and the URR grey. For `edta` there is no such convention, so no colours are written and `build` assigns them.

Where a set has many leaves in few colours, the colours file also fills in `legend_group` (the optional 4th column) so the legend collapses. `cytoband` groups by Giemsa stain, turning several hundred bands into nine legend rows; `asat` groups arrays by suprachromosomal family, also nine rows; `repeatmasker` groups the five RNA leaves, which share one colour; `pave` groups the E5 variants into one row. See [`build`'s legend section](build.md#grouping-the-legend).

## Options

Common to every subcommand:

| Option | Description |
| --- | --- |
| `--output PATH` | Output feature-set BED. Required. |
| `--hierarchy PATH` | Output `child<TAB>parent` hierarchy. Required except on `fai`, where it is optional. |
| `--name NAME` | Feature-set name used in the printed stanza. Defaults to the conventional name for the format. |
| `--seqid-map PATH` | Two-column `old new` seqid table. |
| `--rename-prefix OLD:NEW` | Leading-prefix seqid rewrite. |
| `--force` | Overwrite existing outputs. Without it, `prep-bed` refuses rather than clobbering. |

Format-specific options of note:

| Subcommand | Option | Description |
| --- | --- | --- |
| `repeatmasker`, `edta` | `--background LABEL` | Gap-fill label named in the stanza (default `nonrepeat`). |
| `repeatmasker`, `edta`, `cytoband`, `satellite` | `--priority PATH` | Also write the tree as a 3-column priority file. |
| `gff-gene` | `--feature-type TYPE` | Column-3 type read as an exon (default `exon`). |
| `cytoband` | `--primary-pattern REGEX` | Which seqids carry banding (default `^chr([0-9]+\|X\|Y)$`). |
| `censat` | `--priority PATH` | Priority file ranking centromeric over rDNA over arm. |
| `asat` | `--class CLASS` | α-satellite class to include: `active_hor`, `hor`, `dhor`, `mon`. Repeatable; default all four. |
| `asat` | `--priority PATH` | Priority file ranking `alpha_hor` over `dhor` over `mon`, with the background named explicitly. |
| `asat` | `--background LABEL` | Gap-fill label named in the stanza and the priority file (default `background`). |
| `satellite` | `--satellite LABEL` | Leaf label for the bands, e.g. `CEN180` or `aSat`. |
| `satellite` | `--merge-gap N` | Merge monomers within N bases into one band (default 10). |
| `satellite` | `--cluster-gap N` | Gap for clustering bands into the centromere core (default 500000). |
| `pave` | `--priority PATH` | Priority file ranking the ORFs in genome order. |
| `pave` | `--fasta PATH` | Also write the genome sequences the records carry. |
| `pave` | `--taxonomy PATH` | Also write the ICTV lineage as a hierarchy, for a `type` set. |
| `pave` | `--background LABEL` | Gap-fill label named in the stanza (default `intergenic`). |

Run `karyoscope prep-bed SUBCOMMAND --help` for the full list.

## Notes on individual converters

### `repeatmasker`

Reads either the native `.out` table or the UCSC BED repackaging of it. The dialect is detected from the file's own content and named in the output.

Leaves are the RepeatMasker classes, with the `?` uncertainty marker stripped (`DNA?` is still `DNA`).

A class the converter has no leaf for is labelled `other_repeat` and reported. It is not folded into `Unknown`, which is a real RepeatMasker class meaning the element could not be classified — a distinct case from a class this converter has no leaf for. It is also not dropped, which would leave those bases to the `nonrepeat` gap-fill. `other_repeat` joins the hierarchy and the palette only when something lands on it, so a file whose classes are all recognised produces the reference tree unchanged.

### `gff-gene`

Exons are read directly; introns are derived **per transcript** as the gaps between that transcript's consecutive exons; everything else is intergenic. Where transcripts disagree about a base the more specific label wins — `exon` over `intron` over `intergenic` — so alternative splicing never double-labels a base.

Deriving introns from a chromosome-wide exon list instead, rather than per transcript, makes every gap between neighbouring *genes* an intron. On a compact genome the error is large: it moves *A. thaliana* from ~21% intron to ~42%.

Only records of `--feature-type` (default `exon`) are read. If your annotation puts UTRs in separate `five_prime_UTR`/`three_prime_UTR` records with no matching `exon`, those bases are not exonic.

The result tiles every sequence in `--lengths`, so the stanza sets `background: null`.

### `cytoband`

Band labels keep their chromosome (`chr1` + `p36.33` → `1p36.33`) so no label is ambiguous across chromosomes, and the hierarchy nests three deep: chromosome → band group → band. A band with no sub-band (`1p33`) hangs directly off its chromosome.

UCSC ships the band column bare in both shapes — the 5-column golden-path table and the 6-column `cytoBandMapped` BED, where the qualified name sits in a *separate* column that is ignored and rebuilt, so one code path serves both. A file that puts the qualified name in column 4 directly is also handled: a band already beginning with its own chromosome is left alone, rather than becoming `11p36.33`. Bare cytogenetic bands always start with the arm, `p` or `q`, which distinguishes the two cases.

The hierarchy is checked before anything is written: every band must be a node with a clean path to the root, and no label may be both a band and a group.

### `censat`

For assemblies with a CenSat annotation (human), which already names its features. CenSat qualifies every label with the specific arrays it contains — `gSat(TAR1)`, `hor(S1C10H1-B)` — so the part before the parenthesis is taken as the leaf, collapsing several hundred distinct values onto the 14 features the hierarchy names. Abutting rows that end up sharing a label are merged.

The centromere is located from the `ct` (centromeric transition) features that bracket it: first `ct` start to last `ct` end. Where a sequence has no `ct`, the extent of all centromeric features is used instead. Everything outside the annotated features is then `p_arm` or `q_arm` by which side of that boundary it falls on — so pericentromeric remnants far out on an arm keep their own labels and only the gaps around them become arm.

`--priority` writes the CenSat v2.1 ranking (centromeric 1, rDNA 2, arm 3), which is what the shipped CHM13 v2 `region` set was built with.

Sequences with no CenSat annotation at all are left uncovered and reported for `exclude:`.

### `asat`

The same CenSat input as `censat`, read per array instead of per class: `hor_1_5(S1C1/5/19H1L)` becomes a record labelled `S1C1_5_19H1L`, where `censat` would label it `hor`. `/` becomes `_`, because `build` writes one FASTA per leaf named by the label.

The hierarchy mirrors the shipped `region` set, so the two describe the same arrays under the same names:

```
asat            → categorized
  alpha_hor     → asat
    active_hor  → alpha_hor
    hor         → alpha_hor
    dhor        → alpha_hor
  mon           → asat
```

Arrays sit flat under their class. A label that appears in both a `hor` and a `dhor` interval is placed under `dhor`.

CenSat names more than one array on intervals where two arrays' sequence is interleaved. Each named array gets its own record over the full interval, so a k-mer found in both resolves to their common ancestor. **Do not pass `build --flatten` for this set** — flattening assigns each base to a single label, which puts the shared sequence on whichever array sorts first.

Continuation names are expanded before labels are read. `hor_1_1(S3C1H2-A,B,C)` names S3C1H2-A, S3C1H2-B and S3C1H2-C; a bare capital continues the previous name with its variant token (`-A`, or the `L` of a live array) replaced. Splitting on the comma alone yields leaves called `A`, `B` and `C`.

Every class is included by default. Dropping one leaves that α-satellite to `build`'s gap-fill, and the gap-fill leaf sits at the hierarchy root, so every k-mer a named array shares with it resolves to the root. The shared sequence is the conserved α-satellite monomer core: on CHM13, excluding `mon` and `dhor` puts 36.9% of array bases on the root, against 4.6% with them included. See [background and the hierarchy root](build.md#background-and-the-hierarchy-root).

Arrays are left flat within their class because the structure among them is a phylogeny no annotation file carries. To add it, derive the tree separately — mashtree over the per-array sequences, for example — and replace the flat `alpha_hor` level. Nothing else changes.

### `satellite`

Monomers within `--merge-gap` coalesce into bands; the default bridges the 1–2 bp monomer-boundary artefacts tandem-repeat finders emit without swallowing real interior insertions, which run to kilobases. Bands then cluster within `--cluster-gap` and the **densest** cluster is taken as the centromere core — a raw min-to-max extent would be dragged across the arm by scattered pericentromeric remnants.

Arms are assigned by coordinate, so this assumes the assembly is oriented short-arm-first. Keep the labels `p_arm` and `q_arm`: `centromere_detection` treats anything else as the centromere catch-all.

Input coordinates follow the file suffix — `.gff`/`.gff3` are read 1-based inclusive, anything else 0-based half-open.

### `fai`

One whole-length record per sequence, labelled with its own name, in the `.fai`'s own order. No grouping hierarchy is derived: how sequences group (autosome vs sex vs organelle, metacentric vs acrocentric, haplotype) is organism-specific curation a `.fai` cannot supply. `--hierarchy` gives a flat list to edit by hand.

Keep non-karyotype sequences out with `build`'s `exclude:`, not by filtering here, so every feature set agrees about what exists.

### `pave`

Reads the JSON that [PaVE](https://pave.niaid.nih.gov/)'s `/api/genome/{id}` endpoint returns — a JSON array of them, a single one, or an object wrapping them under `data`. A summary from the `/api/genome` list endpoint carries no sequence or coordinates and is rejected by name.

Leaves are the ORFs — `E6`, `E7`, `E1`, `E2`, `E5` (as the genus-specific `E5_ALPHA`, `E5_BETA`, `E5_GAMMA`, `E5_DELTA` and `E5_OTHER`), `E10`, `L2`, `L1` — plus the `URR`, grouped `early` and `late`. A feature name with no leaf raises rather than being dropped, which would hand its bases to the gap-fill.

Two kinds of feature are excluded. The spliced transcripts `E1^E4`, `E8^E2` and `E6*` lie wholly inside the ORFs they are spliced from — across PaVE's 224 human reference genomes, 100% of their bases — so under a ranking that prefers the primary ORF they hold no sequence. The `E1BS` and `E2BS` binding motifs are 12–20 bp, shorter than any usable k.

Papillomavirus genomes are circular, and a feature spanning the origin arrives as two locations (`7157..7906` and `1..103`). Both become BED records with the same label. The k-mers that cross the junction are not present in a linear FASTA.

The reading frames overlap, so the set is built in priority mode rather than flattened: 10.2% of bases carry more than one annotated feature, and 1.2% once the spliced transcripts are dropped. `--priority` writes the ranking in genome order, which resolves each observed overlap toward the earlier ORF — `E6` over `E7`, `E7` over `E1`, `E1` over `E2`, `L2` over `L1`, and `E10` over `URR` through `early` outranking it.

A record also carries the genome sequence and the ICTV lineage, so `--fasta` writes the FASTA for `build`'s `sequence:` (bgzip and index it) and `--taxonomy` writes the family → genus → species → type hierarchy for a `type` feature set whose BED comes from `fai`. Taxon names have their spaces replaced with underscores, since hierarchy files are whitespace-separated. The constant `Papillomaviridae` level is dropped.

## Examples

```bash
# A repeat set from RepeatMasker, with the reference palette and legend groups
karyoscope prep-bed repeatmasker --input CHM13.RepeatMasker.bed \
    --output repeat.bed --hierarchy repeat.tsv --colors repeat_colors.tsv

# A gene set from a RefSeq GTF whose seqids are accessions
karyoscope prep-bed gff-gene --input CHM13.gtf.gz --lengths CHM13.fa.gz.fai \
    --seqid-map refseq_to_chr.tsv --output gene.bed --hierarchy gene.tsv

# A cytoband set, appending the stanza (and its exclude: list) to a spec
karyoscope prep-bed cytoband --input hg38.cytoBand.txt.gz --lengths hg38.fa.fai \
    --output cytoband.bed --hierarchy cytoband.tsv \
    --colors cytoband_colors.tsv >> build.yaml

# A centromeric region set for a non-human assembly
karyoscope prep-bed satellite --input ColCEN_CEN180.gff3 --lengths ColCEN.fa.gz.fai \
    --satellite CEN180 --rename-prefix Col-CEN_chr:Chr \
    --output region.bed --hierarchy region.tsv

# A per-array alpha-satellite set from the same CenSat file the region set uses
karyoscope prep-bed asat --input chm13v2.0.cenSatv2.1.bed \
    --output asat.bed --hierarchy asat.tsv

# A papillomavirus gene set, plus the FASTA and taxonomy from the same records
karyoscope prep-bed pave --input pave_human_ref.json \
    --output gene.bed --hierarchy gene.tsv --priority gene.priority.txt \
    --fasta HPV.fasta --taxonomy type.hierarchy.txt

# Then build from the assembled spec
karyoscope build --spec build.yaml
```

## See also

- [`karyoscope build`](build.md) — turn these BEDs into an HKS database.
- [Database recipes](../recipes/) — end-to-end worked examples, from download URL to built database.
