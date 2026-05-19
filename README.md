<div align="center">

<!-- TODO: replace with the project logo once designed -->
<!-- <img src="assets/logo.svg" alt="KaryoScope" width="400"/> -->

# KaryoScope

**KaryoScope: rapid, alignment-free sequence annotation for the pangenome era.**

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![CI](https://github.com/barthel-lab/KaryoScope/actions/workflows/ci.yml/badge.svg)](https://github.com/barthel-lab/KaryoScope/actions/workflows/ci.yml)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2026.05.15.725544-red)](https://doi.org/10.64898/2026.05.15.725544)
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

- **Pangenome-scale throughput.** Annotates a complete human haplotype in ~8 minutes on a standard workstation, scaling to cohorts of hundreds of phased assemblies.
- **Base-pair resolution across the entire genome.** Performs well in the satellite-dense centromeres, subtelomeres, and acrocentric short arms where alignment-based pipelines suffer from reference bias and ambiguous mappings.
- **Multiple feature classes in a single pass.** The same *k*-mer can carry labels across feature sets simultaneously, so a single position can be annotated as belonging to a specific chromosome, satellite family, repeat class, and gene at once.
- **Extensible.** Any annotation that tiles a reference of interest can serve as a feature set.

## Installation

> Installation via Bioconda is planned. For now, install from source.

KaryoScope requires Python ≥3.10 and several external tools (`KMC`, `bgzip`, `tabix`, `seqtk`). The simplest setup is a dedicated conda environment:

```bash
git clone https://github.com/barthel-lab/KaryoScope.git
cd KaryoScope

# Create a dedicated environment with Python and the bioinformatics tools
conda create -n karyoscope -c conda-forge -c bioconda \
    python=3.12 pip kmc htslib seqtk
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

```bash
# 1. Download the recommended human reference database (~17 GB, one-time)
karyoscope download

# 2. Download HG002 example data (~3 GB, one-time) — or skip if you have your own assembly
bash examples/quickstart_hg002.sh

# 3. Annotate the assembly
karyoscope annotate --input HG002.maternal.fa.gz --outdir results/

# 4. Render a karyotype
karyoscope karyotype --annotation results/HG002.maternal.chromosome.smoothed.bed.gz \
                     --output HG002.maternal.chromosome.smoothed.genome.karyotype.svg
```

## Commands

| Command | Purpose |
|---|---|
| [`karyoscope download`](docs/commands/download.md) | Acquire pre-built databases |
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
