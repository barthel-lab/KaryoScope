# Recipe: `HKS_human_CHM13_v2`

Rebuild the human T2T-CHM13v2.0 database's `chromosome`, `region`, `repeat` and `gene` feature sets from published sources.

> The shipped database also carries `acrocentric` and `subtelomeric` feature sets, which are outside this recipe.

**Requirements:** `karyoscope`, `samtools`, `hks`, and UCSC's [`bigBedToBed`](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed).
**Cost:** roughly 25 GB of downloads, ~1 GB of intermediates, and a build that takes minutes at `-t 16` but needs substantial RAM — see [Resource requirements](../commands/build.md#resource-requirements).

Three files in [`files/`](files/) are used below. They are shipped because nothing derives them: two encode curation (which chromosomes are acrocentric; which repeat class outranks which) and one is an accession mapping.

## 1. Genome

```bash
curl -L -o chm13v2.0.fa.gz \
  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz

# Re-compress with bgzip, which build and samtools need for random access
zcat chm13v2.0.fa.gz | bgzip > CHM13.fasta.gz
samtools faidx CHM13.fasta.gz
```

Verify by **content**, not by file — bgzip changes the container:

```bash
zcat CHM13.fasta.gz | sha256sum
# 15a4ba1246f6021a89699bf5083da7f2bad3f79c86acd7bc1eb0ca3a13164e85
```

## 2. `chromosome`

One record per sequence, labelled with its own name. Derived from the `.fai`, so it needs no separate download.

```bash
karyoscope prep-bed fai \
    --lengths CHM13.fasta.gz.fai \
    --output chm13_chromosome.bed

sha256sum chm13_chromosome.bed
# da23139f348285a17085e264f8fac86644327965c647a12ced3ad810ea7e1827   (25 rows)
```

The shipped database groups these by centromere position — metacentric, submetacentric, acrocentric, sex. That is human cytogenetics and cannot be read off a `.fai`, so it is provided rather than derived, and passed to [the build](#6-build) below:

```bash
cp docs/recipes/files/chm13v2_chromosome.hierarchy.txt .
```

## 3. `region` — centromeric satellites

```bash
curl -L -o chm13v2.0.cenSatv2.1.bed \
  https://raw.githubusercontent.com/hloucks/CenSatData/refs/heads/main/CHM13/chm13v2.0.cenSatv2.1.bed

sha256sum chm13v2.0.cenSatv2.1.bed
# c6f35d83e4b28f33807c2f65df76b78444f246aa9fb82a3c262bde35e6a6a78c
```

```bash
karyoscope prep-bed censat \
    --input chm13v2.0.cenSatv2.1.bed \
    --lengths CHM13.fasta.gz.fai \
    --output chm13_region.bed \
    --hierarchy chm13_region.tsv \
    --priority chm13_region.priority.txt

sha256sum chm13_region.bed
# c3b0959f9221d89a854e680186589e428c101eb4de069e1c2dab01fb58540bd3   (3,452 rows)
```

`prep-bed` writes the hierarchy and priority files itself — there is nothing to supply.

CenSat labels every feature with the specific arrays it contains — `gSat(TAR1)`, `hor(S1C10H1-B)` — several hundred distinct values. `prep-bed` keeps the part before the parenthesis, giving the 14 leaves the hierarchy names, then labels everything else `p_arm` or `q_arm` by which side of the centromere it falls on. The centromere comes from the `ct` (centromeric transition) features that bracket it.

`chrM` has no CenSat annotation and is reported for `exclude:` rather than labelled.

## 4. `repeat` — RepeatMasker

The source is the UCSC bigBed, **not** the similarly-named BED in the T2T annotation bucket. Both are RepeatMasker 4.1.2p1 2022Apr14; the bucket's BED lists one row per fragment, while the bigBed joins fragments into blocked records, and only the latter matches this database.

```bash
curl -L -o chm13v2.0_rmsk.bb \
  https://hgdownload.soe.ucsc.edu/gbdb/hs1/t2tRepeatMasker/chm13v2.0_rmsk.bb

sha256sum chm13v2.0_rmsk.bb
# 92dfe2d85113752ebf140d9424d5979e56de64718a40b5c2374e1c1c454482f0
```

```bash
bigBedToBed chm13v2.0_rmsk.bb CHM13.RepeatMasker.bed

sha256sum CHM13.RepeatMasker.bed
# e061ca26b34aa9e9c0b6cbc70f0d6faefd123000af6cd2716b9452e4d717b0ed   (4,636,653 rows)
```

```bash
karyoscope prep-bed repeatmasker \
    --input CHM13.RepeatMasker.bed \
    --output chm13_repeat.bed \
    --hierarchy chm13_repeat.tsv \
    --colors chm13_repeat.colors.tsv

sha256sum chm13_repeat.bed
# 067d8020ec3480aafb30a9ff07280c1b2cb8ed58174d895212ebf9b2b96fc10e   (4,636,653 rows)
```

`prep-bed` writes the hierarchy and colours files itself — there is nothing to supply.

RepeatMasker names each element `name#class/family`, as in `L1MA#LINE/L1`. `prep-bed` keeps the class, giving the 15 leaves the hierarchy names. The colours file groups the five RNA leaves into a single legend row, since they share a colour.

Repeat annotations overlap heavily, and this database resolves that by assigning each base to the highest-priority class. The priority order is provided, and passed to [the build](#6-build) below:

```bash
cp docs/recipes/files/chm13v2_repeat.priority.txt .
```

## 5. `gene` — RefSeq

```bash
curl -L -o CHM13.gtf.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz

sha256sum CHM13.gtf.gz
# fd8a27c06da23b4defc140d1bb03f5f0eb8b5ba780eeee6730617fe709c6a16e
```

RefSeq names sequences by accession (`NC_060925.1`), so it needs a mapping to the assembly's `chr` names:

```bash
cp docs/recipes/files/chm13v2_refseq_to_chr.tsv .

karyoscope prep-bed gff-gene \
    --input CHM13.gtf.gz \
    --lengths CHM13.fasta.gz.fai \
    --seqid-map chm13v2_refseq_to_chr.tsv \
    --output chm13_gene.bed \
    --hierarchy chm13_gene.tsv \
    --colors chm13_gene.colors.tsv

sha256sum chm13_gene.bed
# 03ae30a1a9c90a574b96eb4dfaa34a0c834ef51f5b3eec4181f39d4e8c9063e9   (612,909 rows)
```

`prep-bed` writes the hierarchy and colours files itself — there is nothing to supply.

Introns are derived per transcript from that transcript's own consecutive exons; where transcripts disagree, `exon` beats `intron` beats `intergenic`. The result tiles every sequence, so the set needs no gap-fill.

## 6. Build

```yaml
# build.yaml
id: HKS_human_CHM13_v2
version: "2.0.0"
sequence: CHM13.fasta.gz
kmer:
  s: 31
build:
  threads: 16
  mem_gigas: 8
# chrM is not a karyotype chromosome and no feature set covers it, so it reads
# as "none" everywhere rather than being claimed by a placeholder label.
exclude: [chrM]
feature_sets:
  - name: chromosome
    bed: chm13_chromosome.bed
    hierarchy: chm13v2_chromosome.hierarchy.txt
    background: null            # one record per sequence already tiles everything
  - name: region
    bed: chm13_region.bed
    hierarchy: chm13_region.tsv
    priority: chm13_region.priority.txt
    background: null            # the arm split already tiles every base
  - name: repeat
    bed: chm13_repeat.bed
    hierarchy: chm13_repeat.tsv
    priority: chm13v2_repeat.priority.txt
    flatten: true               # one class per base, by the priority order
    colors: chm13_repeat.colors.tsv
    background: nonrepeat
  - name: gene
    bed: chm13_gene.bed
    hierarchy: chm13_gene.tsv
    colors: chm13_gene.colors.tsv
    background: null            # exon/intron/intergenic already tiles every base
roles:
  chromosome_assignment: chromosome
  region_assignment: region
  centromere_detection: region
smoothing:
  recommended_window_bp: 1000
```

```bash
karyoscope build --spec build.yaml
karyoscope info HKS_human_CHM13_v2
```

All four feature sets reproduce the shipped database's inputs exactly.

The one difference throughout is `chrM`. The original BEDs gave it a literal `exclude` *label*, a convention that predates `build`'s `exclude:` list. This recipe uses `exclude:`. A label claims the sequence for that feature set; `exclude:` leaves it uncovered, so it reads as `none` everywhere.

## See also

- [`karyoscope prep-bed`](../commands/prep-bed.md)
- [`karyoscope build`](../commands/build.md)
- [Recipe: Arabidopsis Col-CEN](arabidopsis-colcen.md)
