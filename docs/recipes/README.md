# Database recipes

Complete, runnable recipes for rebuilding KaryoScope databases from published source annotations — starting from the download.

| Recipe | Database | Feature sets |
| --- | --- | --- |
| [human-chm13v2.md](human-chm13v2.md) | `HKS_human_CHM13_v2` | `chromosome`, `region`, `repeat`, `gene` |
| [arabidopsis-colcen.md](arabidopsis-colcen.md) | `HKS_arabidopsis_ColCEN` | `chromosome`, `region`, `repeat`, `gene` |

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

This matters more than it looks. Re-compressing a file changes its bytes without changing a single base, so a `.gz` checksum mismatch may mean nothing at all — while a matching *content* checksum is proof. The CHM13 genome is exactly this case: the copy these databases were built from was re-bgzipped locally, so its `.gz` never matched upstream.

## Two traps worth knowing about

**The same annotation is often published more than once, in more than one form.** Two real examples, both of which silently produce a *different* feature set rather than an error:

- The T2T bucket's `chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed` and UCSC's `chm13v2.0_rmsk.bb` are the same RepeatMasker run. The first lists one row per *fragment*; the second joins fragments into blocked records. The CHM13 `repeat` set uses the second.
- The Col-CEN gene annotation exists as GTF-style attributes on GitHub and GFF3-style on TAIR. Same content; only the punctuation differs.

This is why every input here is pinned by checksum rather than by name.

**A file's coordinates belong to a specific assembly.** An annotation built on a different version of "the same" genome produces coordinates that are in-bounds and plausible but subtly wrong. If you substitute an input, verify the assembly identity by comparing per-sequence sequence md5, not by spot-checking that a feature lands somewhere sensible.

## Supporting files

[`files/`](files/) holds the few inputs that nothing derives, because they encode curation or a naming convention rather than data:

| File | Used by |
| --- | --- |
| `chm13v2_chromosome.hierarchy.txt` | which human chromosomes are metacentric / submetacentric / acrocentric / sex |
| `colcen_chromosome.hierarchy.txt` | which Col-CEN sequences are nuclear |
| `chm13v2_repeat.priority.txt` | the published order in which repeat classes win an overlap |
| `chm13v2_refseq_to_chr.tsv` | RefSeq accession → chromosome name |

Every other hierarchy and priority file in these recipes is written by `prep-bed` as it converts, so there is nothing to supply.

## A note on exact reproduction

These recipes rebuild the feature sets, and the BEDs they produce are byte-identical to the ones the shipped databases were built from — with one deliberate exception, described in each recipe: the shipped BEDs label the mitochondrion with a literal `exclude` *label*, a convention that predates `build`'s `exclude:` list. The recipes use `exclude:`, which is better: a label claims the sequence for the feature set, whereas `exclude:` leaves it genuinely uncovered so it reads as `none` everywhere.
