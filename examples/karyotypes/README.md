# Reference karyotypes

Karyotype plots for six assemblies, to compare your own output against. If a
plot of your assembly looks unlike anything here, that is usually informative —
either about the assembly or about how the run was configured.

Every plot was produced by a single `karyotype` command, which runs the whole
cascade (annotate → scaffold → centromere detection → bin → render):

```bash
karyoscope karyotype -i hap1=asm.hap1.fa.gz -i hap2=asm.hap2.fa.gz \
    --db HKS_human_CHM13_v2 --sample-label MYSAMPLE \
    --format svg --format png -t 16 -o out/
```

All human plots use the shipped `HKS_human_CHM13_v2` database. Arabidopsis uses
a database built locally with [`karyoscope build`](../../docs/commands/build.md).

## The three views

Each sample has three panels, chosen because they answer different questions.

| Panel | What it answers |
| --- | --- |
| **`genome` view, `chromosome` feature set** | Which reference chromosome does each contig come from? Every bar should be one uniform colour. |
| **`centromere` view, `region` feature set** | Is the centromeric satellite organisation intact — α-satellite HOR arrays, the HSat blocks, the flanking monomeric α-satellite? |
| **`subtelomere` view, `subtelomeric` feature set** | Is every chromosome arm capped by a canonical telomere? |

The `genome`/`chromosome` panel is the one to look at first. Because each
chromosome has its own colour, **a bar that is not one solid colour is a contig
carrying sequence from somewhere else** — a translocation, a misjoin, or a
misassembly. That single property is what makes these plots comparable at a
glance.

The `subtelomere` panels are mostly grey (`nonsubtelomeric`), and that is the
point: the coloured caps at the extreme ends are the telomeres, so a *missing*
cap stands out.

## The samples

| Sample | Assembly | What it demonstrates |
| --- | --- | --- |
| [CHM13](chm13.genome.chromosome.png) | T2T-CHM13 v2.0 | Haploid, essentially complete. The reference baseline. |
| [HG002](hg002.genome.chromosome.png) | HG002 v1.1 | Diploid, near-T2T. Two haplotypes per chromosome. |
| [HG008N](hg008n.genome.chromosome.png) | HG008N v6.3 | Diploid normal — the control half of a tumour/normal pair. |
| [HG008T](hg008t.genome.chromosome.png) | HG008T v3.2 | Matched tumour. Extensive somatic rearrangement. |
| [NA19185](na19185.genome.chromosome.png) | HPRC release 2 | A routine population assembly, not curated to T2T. |
| [Arabidopsis](arabidopsis.genome.chromosome.png) | Col-CEN v1.2 | Non-human, non-mammalian, own database. |

### CHM13 — the baseline

`chm13v2.0` is the haploid assembly the human k-mer database is built from, so
this is the best case the method can produce: 24 chromosomes, each a single
uniform colour, a telomere marker at both ends of every one. Use it to
calibrate your eye. Any departure in your own plot is a real difference, not a
rendering artefact.

Because the database is derived from this assembly, this plot also doubles as a
self-consistency check on the database.

### HG002 — diploid

A near-T2T diploid, drawn as two haplotype columns per chromosome. It is male,
so `chrX` has only `h2` and `chrY` only `h1` — with `--sex unknown`, KaryoScope
draws whichever haplotypes carry data rather than assuming.

Its legend picks up two extra entries, `autosome` and `acrocentric`. Those are
*internal* hierarchy nodes, not chromosomes: where a k-mer is ambiguous between
several chromosomes, it resolves to their common ancestor instead of guessing.
Small amounts of this near acrocentric short arms are normal and expected.

### HG008N and HG008T — a tumour/normal pair

Look at these two together; the contrast is the point.

**HG008N** (normal) is what an unrearranged diploid genome looks like: every
chromosome one clean colour on both haplotypes, with only small flecks near the
acrocentric p-arms of chr13/14/15/21/22, where rDNA and satellite are shared
between chromosomes and genuinely ambiguous.

**HG008T** (tumour, same individual) is visibly rearranged. Multiple
chromosomes carry long segments in another chromosome's colour, and some
columns are absent altogether. Because normal and tumour were annotated against
the same database with the same settings, differences between the two plots are
differences between the genomes.

### NA19185 — a routine population assembly

An HPRC release-2 diploid that had no part in building the database, so it
shows how the method behaves on an ordinary sample rather than a reference.
Most chromosomes are clean, but the acrocentrics are fragmented into many short
contigs and some chromosomes are less complete than in the curated assemblies
above. **This is the most representative example of what your own data will
probably look like** — the curated assemblies are the exception, not the norm.

### Arabidopsis — a non-human genome and a user-built database

Five chromosomes, annotated against `HKS_arabidopsis_ColCEN`, built from a
genome plus BED files with [`karyoscope build`](../../docs/commands/build.md).
It demonstrates that nothing in the method is human-specific.

Two things differ from the human examples, both instructive:

- Its third panel is `subtelomere`/**`region`**, not `subtelomeric`, because
  that database has no `subtelomeric` feature set. **Feature sets are a
  property of the database, not of KaryoScope** — you build the ones your
  organism needs.
- Plant telomeres are `TTTAGGG`, so the run needs `--telo-motif CCCTAAA`. The
  human default finds no telomeres in a plant genome.

The centromere panel shows CEN180 satellite arrays at the expected positions,
and chr2 and chr4 have short NOR-bearing p-arms.

## Regenerating these

`ks_examples_karyotype.sh` in the KaryoScope working notes drives all six.
Rendering is deterministic given the same assembly, database, and KaryoScope
version, so a rebuild should reproduce these images — with the caveat that
changes to the renderer or to binning will legitimately change them.
