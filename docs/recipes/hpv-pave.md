# Recipe: `HKS_hpv_PaVE`

Build a human papillomavirus database from the [PaVE](https://pave.niaid.nih.gov/) reference genomes: a `type` set naming the HPV type and a `gene` set naming the ORF.

224 genomes totalling 1.68 Mb, mean length 7,514 bp. The build takes about 2 s at `-t 8` and produces a 10.5 MB database.

**Requirements:** `karyoscope`, `samtools`, `hks`, and Python 3 for the download.

One file in [`files/`](files/) is used below — the download script. Every hierarchy, priority and colours file here is written by `prep-bed`.

## What this database answers

Each k-mer of a query gets two labels: which papillomavirus type it came from, and which viral gene. Sequence absent from all 224 genomes reads as `novel`, so a host assembly carrying an integrated virus shows the insert against a `novel` background.

## 1. Genome records

PaVE serves one record per genome from `/api/genome/{id}`, carrying the sequence, the ORF coordinates and the ICTV lineage. There is no bulk endpoint that returns all three, so the download walks the reference list:

```bash
cp docs/recipes/files/fetch_pave_human_ref.py .
python3 fetch_pave_human_ref.py

sha256sum pave_human_ref.json
# 2f808c3bfee41aadc06b274419941aee12a0f797946dca1e748bd317cd90348f   (224 records, 3.9 MB)
```

The script sorts records by locus id and serialises with fixed formatting, so two runs against an unchanged database produce byte-identical output. Takes about 100 s.

**PaVE adds genomes.** This checksum pins the database as of 2026-08-07. A mismatch alongside a count other than 224 means the reference set has changed rather than that the download failed; re-pin both.

The list is filtered to `is_ref` human genomes. PaVE also publishes 230 non-reference human genomes — the `nr` putative novel types — which are annotated just as completely but whose species is `Unclassified <Genus>` for all but two, so they enter the taxonomy as one flat bucket per genus. To include them, edit the script in two places: set `includeNonRef=true` in the listing URL and drop the `g["is_ref"]` term from the filter on the next line.

## 2. `gene`, plus the FASTA and the taxonomy

One command reads the records three ways:

```bash
karyoscope prep-bed pave --input pave_human_ref.json \
    --output hpv_gene.bed --hierarchy hpv_gene.tsv \
    --priority hpv_gene.priority.txt --colors hpv_gene.colors.tsv \
    --fasta HPV.fasta --taxonomy hpv_type.hierarchy.txt

sha256sum hpv_gene.bed hpv_gene.tsv hpv_gene.priority.txt hpv_gene.colors.tsv
# 2a9352acaff21e12aa8390c038eafb54c05eaf6552cc405a25de52c89428759b  hpv_gene.bed   (1,757 rows)
# 2b5c25b898e03efeea844fa779c76603e7b45c569122d61ded713d207625c6b5  hpv_gene.tsv   (16 edges)
# c0bfb4780b830781fcfe4cdb9a7c113f2fb88b86c8c56c77e13982ca67304e60  hpv_gene.priority.txt
# cc9c2df8f4841fc874525c53d1a8b37158fe037a184e6772d6fe0534b6583ec2  hpv_gene.colors.tsv
```

Leaves are `E6`, `E7`, `E1`, `E2`, the genus-specific `E5_*`, `E10`, `L2`, `L1` and `URR`, grouped `early` and `late`.

Three features PaVE annotates are not in the set. `E1^E4`, `E8^E2` and `E6*` are spliced from the primary ORFs and lie wholly inside them — 100% of their bases across all 224 genomes — so a ranking that prefers the primary ORF leaves them empty. `E1BS` and `E2BS` are 12–20 bp, shorter than any usable k.

The reading frames overlap: 10.2% of bases carry more than one annotated feature, 1.2% once the spliced transcripts are dropped. The set is therefore built in priority mode rather than flattened, and the priority file ranks the ORFs in genome order — which resolves each overlap toward the earlier ORF (`E6` over `E7`, `E7` over `E1`, `E1` over `E2`, `L2` over `L1`, and `E10` over `URR` through `early` outranking it).

1.56% of bases carry no ORF and go to the `intergenic` gap-fill.

## 3. Genome FASTA

`--fasta` above wrote the sequences the records carry. Re-compress with bgzip, which `build` and `samtools` need for random access:

```bash
bgzip HPV.fasta
samtools faidx HPV.fasta.gz
```

Verify by **content**, not by file — bgzip changes the container:

```bash
zcat HPV.fasta.gz | sha256sum
# 8bc37f684a209f1045c350b5eb3a96276af1cb082a0a799adfb0b6f8f38d3a2c
```

Six genomes carry a single IUPAC ambiguity code each (`Y` ×4, `W`, `N`), in HPV142REF, HPV171REF, HPV172REF, HPV175REF, HPV223REF and HPV228REF. Each one invalidates the 31 k-mers spanning it, so 186 positions are unresolvable before smoothing.

**Papillomavirus genomes are circular, and this FASTA is linear.** A feature spanning the origin becomes two BED records with the same label, but the k-mers crossing the junction cannot be represented — about 30 bp per genome.

## 4. `type`

One record per genome, labelled with its own name, derived from the `.fai`:

```bash
karyoscope prep-bed fai --lengths HPV.fasta.gz.fai --output hpv_type.bed --name type

sha256sum hpv_type.bed hpv_type.hierarchy.txt
# ddb8573f141c8053c9f9aa37188bd0e3b600af95989d2143e71757166e319a18  hpv_type.bed  (224 rows)
# 0a08444b08a6ddcd3ba8cc7272994e361127e9bd01e00fe190d37e83c2fda0ae  hpv_type.hierarchy.txt  (280 edges)
```

The hierarchy is the one `--taxonomy` wrote in step 2: genus → species → type, over 5 genera and 51 species. A `.fai` cannot supply it, but PaVE states it, so it is derived rather than hand-written. The constant `Papillomaviridae` level is dropped, and spaces in taxon names become underscores because hierarchy files are whitespace-separated.

This set needs no priority file. A k-mer shared by several types resolves to their lowest common ancestor, which is the species or genus that contains them.

## 5. Build

```yaml
# build.yaml
id: HKS_hpv_PaVE
version: "1.0.0"
sequence: HPV.fasta.gz
kmer:
  s: 31
build:
  threads: 8
  mem_gigas: 4
feature_sets:
  - name: type
    bed: hpv_type.bed
    hierarchy: hpv_type.hierarchy.txt
    background: null            # one record per genome already tiles everything
  - name: gene
    bed: hpv_gene.bed
    hierarchy: hpv_gene.tsv
    priority: hpv_gene.priority.txt   # the reading frames overlap; genome order wins
    colors: hpv_gene.colors.tsv
    background: intergenic
roles:
  chromosome_assignment: type
  region_assignment: gene
smoothing:
  recommended_window_bp: 100
```

```bash
karyoscope build --spec build.yaml
karyoscope info HKS_hpv_PaVE
```

There is no `centromere_detection` role: a papillomavirus genome has no centromere, and `karyoscope centromeres` does not apply to this database.

This recipe has been run end to end from nothing but the script above: `pave_human_ref.json` and all six derived files match the checksums here, and all five index files — the shared base index and one per feature set — come out bit-identical between runs.

## 6. Annotate

```bash
karyoscope annotate -i sample.fasta.gz --db HKS_hpv_PaVE -o out/ -t 8
```

Self-annotating the reference genomes gives, on the smoothed output, 99.98% of positions labelled with their own type and 99.99% with their own ORF. No position is assigned to a different type, and 16 of the 224 genomes have any non-exact type call at all — each resolving to the correct species or genus where a close relative shares the sequence.

62 positions across 5 genomes take an adjacent ORF's label. Each one's 31-mer also occurs in another genome where PaVE places the ORF boundary at a different offset in the shared sequence; the index stores one label per distinct k-mer, and the priority file breaks the tie.

## Choosing `k`

`s: 31` sets how much sequence a k-mer must match exactly, which trades type resolution against sensitivity to divergence. Measured over the 224 genomes:

| `k` | k-mers resolving to one type | to a species | to a genus | to the root |
| --- | --- | --- | --- | --- |
| 21 | 95.51% | 3.03% | 1.31% | 0.153% |
| 25 | 97.19% | 2.09% | 0.69% | 0.027% |
| 31 | 98.41% | 1.25% | 0.33% | 0.001% |
| 41 | 99.22% | 0.62% | 0.16% | 0.000% |

A divergent isolate matches fewer k-mers exactly at any `k`, and lands on the species or genus rather than the type. Lowering `k` shifts calls from the type toward the clade; raising it does the reverse.

`variable_k` would let one database be queried across that range, but it cannot be combined with priorities, and the `gene` set needs them — so a database offering the sweep would have to carry `type` alone.

## See also

- [`karyoscope prep-bed`](../commands/prep-bed.md)
- [`karyoscope build`](../commands/build.md)
- [Recipe: human CHM13v2](human-chm13v2.md)
- [Recipe: Arabidopsis Col-CEN](arabidopsis-colcen.md)
