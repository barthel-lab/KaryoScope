# Database recipes

Complete, runnable recipes for rebuilding KaryoScope databases from published source annotations — starting from the download.

| Recipe | Database | Feature sets |
| --- | --- | --- |
| [human-chm13v2.md](human-chm13v2.md) | `HKS_human_CHM13_v2` | `chromosome`, `region`, `repeat`, `gene` |
| [arabidopsis-colcen.md](arabidopsis-colcen.md) | `HKS_arabidopsis_ColCEN` | `chromosome`, `region`, `repeat`, `gene` |
| [hpv-pave.md](hpv-pave.md) | `HKS_hpv_PaVE` | `type`, `gene` |

Each recipe gives the exact download URL for every input, a checksum to verify it, the [`prep-bed`](../commands/prep-bed.md) command that converts it, and the [`build`](../commands/build.md) spec that assembles the result.

## Verify your downloads

Every input below is pinned by SHA-256. Check one with:

```bash
sha256sum <file>
```

For **gzipped** files, compare the checksum of the *decompressed* content as well as the file:

```bash
zcat <file.gz> | sha256sum
```

Re-compressing a file changes its bytes without changing a base, so a `.gz` checksum mismatch does not by itself indicate different content, while a matching content checksum does indicate the same content. The CHM13 genome is such a case: the copy these databases were built from was re-bgzipped locally, so its `.gz` does not match upstream.

## Before substituting an input

**The same annotation is often published more than once, in more than one form.** Two examples, each of which produces a *different* feature set rather than an error:

- The T2T bucket's `chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed` and UCSC's `chm13v2.0_rmsk.bb` are the same RepeatMasker run. The first lists one row per *fragment*; the second joins fragments into blocked records. The CHM13 `repeat` set uses the second.
- The Col-CEN gene annotation exists as GTF-style attributes on GitHub and GFF3-style on TAIR. Same content; only the punctuation differs.

This is why every input here is pinned by checksum rather than by name.

**A file's coordinates belong to a specific assembly.** An annotation built on a different version of "the same" genome produces coordinates that are in-bounds and plausible but subtly wrong. If you substitute an input, verify the assembly identity by comparing per-sequence sequence md5, not by spot-checking that a feature lands somewhere sensible.

## Supporting files

[`files/`](files/) holds the few files a recipe needs that nothing derives — curation, a naming convention, or a download that has no single URL:

| File | Used by |
| --- | --- |
| `chm13v2_chromosome.hierarchy.txt` | which human chromosomes are metacentric / submetacentric / acrocentric / sex |
| `colcen_chromosome.hierarchy.txt` | which Col-CEN sequences are nuclear |
| `chm13v2_repeat.priority.txt` | the published order in which repeat classes win an overlap |
| `chm13v2_repeat.order.txt` | which repeat class each base goes to when annotations overlap (the `flatten_order:` for the repeat set) |
| `chm13v2_gene.priority.txt` | how `exon`, `intron` and `intergenic` rank when a k-mer is claimed by more than one |
| `chm13v2_refseq_to_chr.tsv` | RefSeq accession → chromosome name |
| `fetch_pave_human_ref.py` | walks PaVE's per-genome API into the one JSON input the HPV recipe converts |

Every other hierarchy and priority file in these recipes is written by `prep-bed` as it converts, so there is nothing to supply.

## A note on exact reproduction

Where a recipe rebuilds a shipped database, the BEDs it produces are byte-identical to the ones that database was built from — with one deliberate exception, described in each recipe: the shipped BEDs label the mitochondrion with a literal `exclude` *label*, a convention that predates `build`'s `exclude:` list. The recipes use `exclude:`. A label claims the sequence for the feature set; `exclude:` leaves it uncovered, so it reads as `none` everywhere.

The HPV recipe builds a database that is not shipped, so it pins its inputs and outputs by checksum instead.
