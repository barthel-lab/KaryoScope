<div align="center">

<!-- TODO: replace with the project logo once designed -->
<!-- <img src="assets/logo.svg" alt="KaryoScope" width="400"/> -->

# KaryoScope

**KaryoScope: rapid, alignment-free sequence annotation for the pangenome era.**

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2026.05.15.725544-red)](https://doi.org/10.64898/2026.05.15.725544)
<!-- Temporarily hidden while the org's GitHub Actions minutes are
     exhausted for the month; un-comment when they refresh:
[![CI](https://github.com/barthel-lab/KaryoScope/actions/workflows/ci.yml/badge.svg)](https://github.com/barthel-lab/KaryoScope/actions/workflows/ci.yml)
-->
<!-- Add after publication:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.PLACEHOLDER.svg)](https://doi.org/10.5281/zenodo.PLACEHOLDER)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/karyoscope?label=bioconda)](https://bioconda.github.io/recipes/karyoscope/README.html)
-->

</div>

> 🚧 **KaryoScope is under active development in preparation for journal submission.** The user-facing API and command set are being finalized. Expect breaking changes until v1.0.0. Watch [releases](https://github.com/barthel-lab/KaryoScope/releases) for stable versions.

---

## Overview

KaryoScope is an alignment-free annotation tool that assigns each *k*-mer in a query assembly or sequencing read to a feature drawn from one or more user-defined hierarchical feature sets, producing a base-pair resolution annotation in a single pass. Because a feature set is simply any tiling of a reference with labelled regions, KaryoScope is extensible to arbitrary annotation sources, from satellite catalogs and repeat libraries to cytobands, FISH-probe coordinates, and structural-variant breakpoints.

A pre-built database for the human genome is distributed alongside the tool, derived from T2T-CHM13v2.0 with six feature sets covering chromosome of origin, satellite composition, interspersed repeats, subtelomeric structure, gene boundaries, and acrocentric-specific features. From these annotations, KaryoScope produces karyotype visualizations and cytogenetic reports without ever performing read alignment. Additional databases can be built for any reference genome or community-curated annotation source.

<!-- TODO: hero figure of a KaryoScope output, e.g. the HG008T karyotype -->
<!-- <p align="center"><img src="assets/hero_karyotype.png" alt="Example KaryoScope karyotype" width="800"/></p> -->

### Why alignment-free?

- **Pangenome-scale throughput.** Annotates a single feature set on a complete human haplotype in ~8 minutes on a standard workstation, or the full six-feature-set pipeline for a diploid sample in ~30 minutes at 16 threads — scaling to cohorts of hundreds of phased assemblies. The in-progress migration to the [HKS](https://github.com/jnalanko/HKS) *k*-mer indexing backend will further reduce runtime and memory footprint.
- **Base-pair resolution across the entire genome.** Performs well in the satellite-dense centromeres, subtelomeres, and acrocentric short arms where alignment-based pipelines suffer from reference bias and ambiguous mappings.
- **Multiple feature classes in a single pass.** The same *k*-mer can carry labels across feature sets simultaneously, so a single position can be annotated as belonging to a specific chromosome, satellite family, repeat class, and gene at once.
- **Extensible.** Any annotation that tiles a reference of interest can serve as a feature set.

## Installation

> Installation via Bioconda is planned. For now, install from source.

KaryoScope requires Python ≥3.10 and several external tools (`KMC`, `bgzip`, `tabix`, `seqtk`, and `cairo` for PDF/PNG karyotype output). The simplest setup is a dedicated conda environment:

```bash
git clone https://github.com/barthel-lab/KaryoScope.git
cd KaryoScope

# Create a dedicated environment with Python and the bioinformatics tools.
# `samtools` is only needed if you plan to annotate BAM inputs; drop it
# if you only work with FASTA or FASTQ.
conda create -n karyoscope -c conda-forge -c bioconda \
    python=3.12 pip kmc htslib samtools seqtk cairo zlib compilers
conda activate karyoscope

# Install KaryoScope
pip install -e .

# Build the bundled C++ helper (`get_featureIDs`).
# `pip install` is Python-only and does NOT compile the C++ tree.
cd native/get_featureIDs && make && cd ../..
```

The build produces `native/get_featureIDs/build/get_featureIDs`; the
Python wrapper finds it automatically. See [`native/README.md`](native/README.md)
for build-system details (CXX selection, `pkg-config`-driven zlib lookup,
and the macOS + conda `-Wl,-rpath,$CONDA_PREFIX/lib` shim).

## Quick start

This walkthrough uses the [HG002 v1.1 T2T diploid assembly](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz) as input, but any FASTA will work. Substitute your own with `--input <path>` throughout.

```bash
# 1. Download the recommended human reference database (~17 GB, one-time)
karyoscope download

# 2. Download the HG002 v1.1 diploid assembly (~3 GB, one-time)
#    Skip if you already have your own assembly to annotate.
curl -O https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz

# 3. Annotate the assembly. Recommended: at least 16 threads and 50 GB
#    of RAM for human-scale inputs (HG002 runs in ~30 min at -t 16).
#    --no-bgzip keeps the per-feature-set BEDs as plain text for easy
#    inspection; drop it to get the default bgzipped outputs.
#    Accepts FASTA, FASTQ (plain or .gz), and BAM. For BAM, samtools
#    must be on PATH (it's invoked as `samtools fasta` to stream into
#    get_featureIDs). For read-level inputs also pass --no-preserve-order
#    for substantially faster writes.
karyoscope annotate --input hg002v1.1.fasta.gz --outdir results/ --threads 16 --no-bgzip

# 4. Render the three primary karyotype views.
#    --no-scaffolding skips the per-feature-set scaffolded BED rewrite
#    (the expensive step of scaffolding); the scaffold map is still
#    applied at bin time so the renders are equivalent.
#    The first invocation runs the full scaffold + bin + render cascade;
#    the next two reuse the cached intermediates and finish much faster.
COMMON="--input hg002v1.1.fasta.gz --outdir results/ --threads 16 --sex male --no-scaffolding --no-bgzip"
karyoscope karyotype $COMMON --mode genome      --feature-set chromosome
karyoscope karyotype $COMMON --mode centromere  --feature-set region
karyoscope karyotype $COMMON --mode subtelomere --feature-set subtelomeric
```

This produces three SVGs under `results/`:

| File | View | Feature set |
|---|---|---|
| `hg002v1.1.KS_human_CHM13_v2.genome.chromosome.smoothed.karyotype.svg` | Genome view | chromosome |
| `hg002v1.1.KS_human_CHM13_v2.centromere.region.smoothed.karyotype.svg` | Centromere view | region |
| `hg002v1.1.KS_human_CHM13_v2.subtelomere.subtelomeric.smoothed.karyotype.svg` | Subtelomere view | subtelomeric |

Pass `--format pdf` or `--format png` (repeatable) to additionally produce those formats from the SVG.

## Commands

| Command | Purpose |
|---|---|
| [`karyoscope download`](docs/commands/download.md) | Acquire pre-built databases |
| [`karyoscope register`](docs/commands/register.md) | Register a locally-placed database so commands can use it |
| [`karyoscope annotate`](docs/commands/annotate.md) | Annotate sequences with *k*-mer features |
| [`karyoscope scaffold`](docs/commands/scaffold.md) | Order, orient, and rename assembly contigs |
| [`karyoscope bin`](docs/commands/bin.md) | Aggregate base-pair annotations into larger bins |
| [`karyoscope centromeres`](docs/commands/centromeres.md) | Extract centromere coordinates |
| [`karyoscope karyotype`](docs/commands/karyotype.md) | Render karyotype visualization |
| [`karyoscope info`](docs/commands/info.md) | Inspect databases, files, installation |
| [`karyoscope version`](docs/commands/version.md) | Print version and environment info |

Run `karyoscope <command> --help` for full options on any command.

## Documentation

Full documentation is being built. In the meantime, the `--help` output for each command is the authoritative reference.

## Databases

KaryoScope works with pre-built databases distributed via the [KaryoScope registry](https://github.com/barthel-lab/KaryoScope-registry). The current default is `KS_human_CHM13_v2` (~17 GB), built from the T2T-CHM13v2.0 reference.

Browse and download available databases:

```bash
karyoscope download --list
```

Building your own database is supported via `karyoscope build` (coming in v1.0).

## Pre-computed annotations

KaryoScope outputs for the [HPRC Release 2](https://humanpangenome.org/) pangenome samples are hosted by the Human Pangenome Reference Consortium at the [`TGen_HPRCv2_KaryoScope` S3 bucket](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/B5FC8EB1-5B2A-49FE-8421-6D938943DFC9--TGen_HPRCv2_KaryoScope/). Use these to explore HPRC karyotypes without running the pipeline yourself, or as references for downstream analysis.

> **Note**: the currently hosted annotations were generated against a previous version of the KaryoScope database. Updated annotations using the current release will be uploaded as they become available.

Per sample, the bucket contains:

| Path | Contents |
|---|---|
| `<sample>/bed/` | Per-feature-set presmoothed annotations (`<sample>.KaryoScope.v2.0.<feature_set>.bed.gz`) |
| `<sample>/igv/` | Per-feature-set, per-haplotype IGV-ready BEDs with tabix index (`<sample>.KaryoScope.v2.0.<feature_set>.hap<i>.IGV.bed.gz` + `.tbi`) |
| `<sample>/plots/` | Karyotype SVGs: genome view (chromosome feature set), centromere view (region), subtelomere view (subtelomeric) |

## Citation

If you use KaryoScope in your work, please cite our preprint:

> Ranallo-Benavidez TR, Chen YA, Potapova T, Alanko J, Loucks H, Lucas J, Human Pangenome Reference Consortium, Guarracino A, Puglisi SJ, Marchet C, Miga K, Gerton JL, Barthel FP. *KaryoScope: rapid, alignment-free sequence annotation for the pangenome era.* bioRxiv (2026). doi: [10.64898/2026.05.15.725544](https://doi.org/10.64898/2026.05.15.725544)

A `CITATION.cff` file in this repository provides machine-readable citation metadata.

## License

KaryoScope is licensed under [GPL-3.0-or-later](LICENSE) due to its dependency on the GPL-3.0 KMC library. A future release will switch to MIT once we migrate to [HKS](https://github.com/jnalanko/HKS) for *k*-mer indexing.

## Contributing

Contributions, bug reports, and feature requests are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) to get started, and our [Code of Conduct](CODE_OF_CONDUCT.md) for community norms.

## Acknowledgments

Developed in the [Barthel Lab](https://www.barthel-lab.com/) at the Translational Genomics Research Institute (TGen), in collaboration with Jarno Alanko, Simon Puglisi, and Camille Marchet.
