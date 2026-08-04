# Recipe: `HKS_arabidopsis_ColCEN`

Rebuild the *Arabidopsis thaliana* Col-CEN (T2T v1.2) database in full: `chromosome`, `region`, `repeat` and `gene`.

This is the whole database, and it is the cheap one — a 135 Mb genome, about 50 s to build at `-t 16` and ~6 GB peak RSS. It is the best recipe to try first, and a worked example of building KaryoScope for a non-human organism.

**Requirements:** `karyoscope`, `samtools`, `hks`.

One file in [`files/`](files/) is used below — the curated chromosome grouping, which nothing derives. Every other hierarchy here is written by `prep-bed`.

## 1. Genome

```bash
curl -L -o Col-CEN_v1.2.fasta.gz \
  https://raw.githubusercontent.com/schatzlab/Col-CEN/main/v1.2/Col-CEN_v1.2.fasta.gz

sha256sum Col-CEN_v1.2.fasta.gz
# b059bf9b589a7a6cd13c67179b80293b91809c3c61cb8b9393a518619d8b5fa8

samtools faidx Col-CEN_v1.2.fasta.gz
```

Seven sequences: `Chr1`–`Chr5`, plus the mitochondrion `ChrM` and chloroplast `ChrC`.

## 2. `chromosome`

```bash
karyoscope prep-bed fai \
    --lengths Col-CEN_v1.2.fasta.gz.fai \
    --output ath_chromosome.bed

sha256sum ath_chromosome.bed
# 28fb5d5a098ac105be5e2052ef53932e5c47d5817a97757f6e6c8c403b284b9b   (7 rows)
```

The organelles are excluded at build time rather than filtered here, so every feature set agrees about which sequences exist.

Grouping the five nuclear chromosomes is curation — nothing in a `.fai` says which sequences are nuclear — so it is provided rather than derived:

```bash
cp docs/recipes/files/colcen_chromosome.hierarchy.txt .
```

## 3. `region` — CEN180 centromeric satellite

The CEN180 monomer catalog comes from TAIR. Note the endpoint: `/download/file?path=…` is the *web page*, and returns an HTML shell. The file itself is served from `/api/download-files/download?filePath=…`.

```bash
curl -L -o ColCEN_CEN180.gff3 \
  "https://www.arabidopsis.org/api/download-files/download?filePath=Genes/Col-CEN_genome_assembly_release/ColCEN_CEN180.gff3"

sha256sum ColCEN_CEN180.gff3
# ccf54c224050a16dc65362623eee90621e803402bb95f8e7c191444028abedd1
```

```bash
karyoscope prep-bed satellite \
    --input ColCEN_CEN180.gff3 \
    --lengths Col-CEN_v1.2.fasta.gz.fai \
    --satellite CEN180 \
    --rename-prefix Col-CEN_chr:Chr \
    --output ath_region.bed \
    --hierarchy ath_region.tsv

sha256sum ath_region.bed
# 339e5e7e6a74bc6c41adb3b0972f3e02d0b695051909030e685f2e2abdcc7a7d   (2,041 rows)
```

`prep-bed` writes the hierarchy itself — there is nothing to supply.

Three things are worth understanding here, because each is a general problem rather than an Arabidopsis quirk:

- **`--rename-prefix` is required.** This file names sequences `Col-CEN_chr1`, the genome names them `Chr1`. Without the rewrite nothing matches and the feature set comes out empty.
- **It is a monomer catalog, not a feature annotation.** Unlike human CenSat it says only "CEN180 monomer here", 66,131 times. `prep-bed satellite` merges monomers within `--merge-gap` (default 10 bp, which bridges the 1–2 bp monomer-boundary artefacts tandem-repeat finders emit without swallowing real kilobase-scale insertions), then takes the densest cluster of bands as the centromere. A raw min-to-max extent would be dragged across the arm by scattered pericentromeric remnants.
- **The file has an uncommented header row.** `prep-bed` skips rows whose coordinates are not integers, so it is handled — but `grep -v '^#'` alone would let a junk record through.

Arms are assigned by coordinate, which assumes short-arm-first orientation.

## 4. `repeat` — EDTA transposons

```bash
curl -L -o EDTA.TEanno.gff3.gz \
  https://raw.githubusercontent.com/schatzlab/Col-CEN/main/v1.2/t2t-col.20201227.fasta.mod.EDTA.TEanno.gff3.gz

sha256sum EDTA.TEanno.gff3.gz
# d209f344edd10e9e37da498b6625454a0472e24ff975c0979d3d12763cb6bcd6
```

```bash
karyoscope prep-bed edta \
    --input EDTA.TEanno.gff3.gz \
    --output ath_repeat.bed \
    --hierarchy ath_repeat.tsv

sha256sum ath_repeat.bed
# 0ec9610c801e5edd4785f0ad9949e3092e09a033f7d4e35b72347c2e7f9eafd4   (42,927 rows)
```

EDTA's `Classification=` vocabulary aliases each superfamily under both spelled-out and Wicker three-letter names — `DNA/HAT` and `DNA/DTA` are both hAT — so `prep-bed` normalises them to one leaf per superfamily, keeping the BED labels and the hierarchy leaves consistent by construction.

No colours file is emitted: there is no established KaryoScope palette for the EDTA vocabulary, and `build`'s automatic assignment beats an invented one.

## 5. `gene` — Araport11

```bash
curl -L -o Col-CEN_v1.2_genes.araport11.gff3.gz \
  https://raw.githubusercontent.com/schatzlab/Col-CEN/main/v1.2/Col-CEN_v1.2_genes.araport11.gff3.gz

sha256sum Col-CEN_v1.2_genes.araport11.gff3.gz
# baa6ad8aa24f098e4e957682ba945da15fa3736eaacceb02b9a2f5de4ec86e8f
```

```bash
karyoscope prep-bed gff-gene \
    --input Col-CEN_v1.2_genes.araport11.gff3.gz \
    --lengths Col-CEN_v1.2.fasta.gz.fai \
    --output ath_gene.bed \
    --hierarchy ath_gene.tsv

sha256sum ath_gene.bed
# 52e8059b9275f65135495d6406b5d7bb5be3758935754cc58ced5aa49fbcd142   (284,461 rows)
```

Two notes:

- Despite the `.gff3` extension this file uses **GTF-style attributes** (`transcript_id "AT1G01010";`). `prep-bed` reads both syntaxes, so it needs no flag. TAIR publishes the same annotation with true GFF3 attributes as `ColCEN_GENES_Araport11.gff3.gz`; either produces this identical output.
- Introns are derived **per transcript**, from that transcript's own consecutive exons. Deriving them from a merged exon list instead makes every gap between neighbouring genes an intron, which on a compact genome like this one is the difference between ~21% and ~42% intron. The result tiles every sequence, so no gap-fill.

## 6. Build

```yaml
# build.yaml
id: HKS_arabidopsis_ColCEN
version: "1.1.0"
sequence: Col-CEN_v1.2.fasta.gz
kmer:
  s: 31
build:
  threads: 16
  mem_gigas: 8
# Organelles are real sequence but not karyotype chromosomes. Excluding them
# from the whole build means no feature set covers them, so they read as "none"
# uniformly rather than appearing as chromosomes.
exclude: [ChrM, ChrC]
feature_sets:
  - name: chromosome
    bed: ath_chromosome.bed
    hierarchy: colcen_chromosome.hierarchy.txt
    background: null            # already tiles every base
  - name: region
    bed: ath_region.bed
    hierarchy: ath_region.tsv
    background: null            # the CEN180 + arm split tiles every base
  - name: repeat
    bed: ath_repeat.bed
    hierarchy: ath_repeat.tsv
    background: nonrepeat
  - name: gene
    bed: ath_gene.bed            # no hierarchy: the three labels are a flat star,
    background: null             # which build derives on its own
roles:
  chromosome_assignment: chromosome
  region_assignment: region
  centromere_detection: region
smoothing:
  recommended_window_bp: 1000
```

```bash
karyoscope build --spec build.yaml
karyoscope info HKS_arabidopsis_ColCEN
```

This recipe has been run end to end from the URLs above, on nothing but the
downloaded files: all four BEDs match the checksums given here, and all five
index files — the shared base index and one per feature set — come out
bit-identical to the database in use.

The `region` set's arm labels are load-bearing: centromere detection reads anything that is not `p_arm`/`q_arm`/`arm`/`telomere`/`novel` as the centromere catch-all, so `CEN180` and `cen_gap` must be leaves.

## 7. Render a karyotype

```bash
karyoscope karyotype -i Col-CEN_v1.2=Col-CEN_v1.2.fasta.gz \
    --db HKS_arabidopsis_ColCEN \
    --sample-label Col-CEN_v1.2 --telo-motif CCCTAAA \
    --sex unknown --format svg --format png -t 16 -o out/
```

`--telo-motif CCCTAAA` matters: plant telomeres are `TTTAGGG`, and the human default finds none. See [`examples/karyotypes/`](../../examples/karyotypes/) for what this produces.

## See also

- [`karyoscope prep-bed`](../commands/prep-bed.md)
- [`karyoscope build`](../commands/build.md)
- [Recipe: human CHM13v2](human-chm13v2.md)
