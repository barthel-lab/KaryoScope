<div align="center">

<!-- TODO: replace with the project logo once designed -->
<!-- <img src="assets/logo.svg" alt="KaryoScope" width="400"/> -->

# KaryoScope

**KaryoScope: rapid, alignment-free sequence annotation for the pangenome era.**

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2026.05.15.725544-red)](https://doi.org/10.64898/2026.05.15.725544)
[![CI](https://github.com/barthel-lab/KaryoScope/actions/workflows/ci.yml/badge.svg)](https://github.com/barthel-lab/KaryoScope/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20657816.svg)](https://doi.org/10.5281/zenodo.20657816)
<!-- Add after publication:
[![Bioconda](https://img.shields.io/conda/dn/bioconda/karyoscope?label=bioconda)](https://bioconda.github.io/recipes/karyoscope/README.html)
-->

</div>

> ℹ️ **KaryoScope follows [semantic versioning](https://semver.org); v1.0.0 is the first stable release.** The command-line interface is stable: deprecations ship with warnings and a back-compatible transition, and any breaking change will come with a major-version bump. The project is being prepared for journal submission, and new features continue to land between releases — see the [CHANGELOG](CHANGELOG.md) and [releases](https://github.com/barthel-lab/KaryoScope/releases).

---

## Overview

KaryoScope is an alignment-free annotation tool that assigns each *k*-mer in a query assembly or sequencing read to a feature drawn from one or more user-defined hierarchical feature sets, producing a base-pair resolution annotation in a single pass. Because a feature set is simply any tiling of a reference with labelled regions, KaryoScope is extensible to arbitrary annotation sources, from satellite catalogs and repeat libraries to cytobands, FISH-probe coordinates, and structural-variant breakpoints.

A pre-built database for the human genome is distributed alongside the tool, derived from T2T-CHM13v2.0 with six feature sets covering chromosome of origin, satellite composition, interspersed repeats, subtelomeric structure, gene boundaries, and acrocentric-specific features. From these annotations, KaryoScope produces karyotype visualizations and cytogenetic reports without ever performing read alignment. Additional databases can be built for any reference genome or community-curated annotation source.

<!-- TODO: hero figure of a KaryoScope output, e.g. the HG008T karyotype -->
<!-- <p align="center"><img src="assets/hero_karyotype.png" alt="Example KaryoScope karyotype" width="800"/></p> -->

### Why alignment-free?

- **Pangenome-scale throughput.** Annotates a single feature set on a complete human haplotype in ~8 minutes on a standard workstation, or the full six-feature-set pipeline for a diploid sample in ~30 minutes at 16 threads — scaling to cohorts of hundreds of phased assemblies. The [HKS](https://github.com/jnalanko/HKS) *k*-mer indexing backend, now available alongside KMC, annotates all six feature sets for a human haplotype in ~7–9 minutes at a ~10 GB memory peak — roughly 2.5–3× faster than KMC and about a third of the RAM (measured at 16 threads on the T2T-CHM13v2.0, HG008N, and HG008T assemblies).
- **Base-pair resolution across the entire genome.** Performs well in the satellite-dense centromeres, subtelomeres, and acrocentric short arms where alignment-based pipelines suffer from reference bias and ambiguous mappings.
- **Multiple feature classes in a single pass.** The same *k*-mer can carry labels across feature sets simultaneously, so a single position can be annotated as belonging to a specific chromosome, satellite family, repeat class, and gene at once.
- **Extensible.** Any annotation that tiles a reference of interest can serve as a feature set.

## System requirements

**Operating systems.** Linux (x86-64) and macOS (Apple Silicon and Intel). KaryoScope is exercised in continuous integration on the current Ubuntu and macOS GitHub Actions runners and developed day to day on macOS (Apple Silicon). There is no native Windows build; on Windows, run it under WSL2.

**Python.** Python ≥ 3.10 (tested on 3.10, 3.11, 3.12, and 3.13).

**Python dependencies.** Installed automatically by `pip`. Minimum versions are pinned in [`pyproject.toml`](pyproject.toml); the versions the current release is tested against are shown alongside.

| Package | Minimum | Tested |
|---|---|---|
| click | ≥ 8.1 | 8.4.0 |
| drawsvg | ≥ 2.4 | 2.4.1 |
| cairosvg | ≥ 2.7 | 2.9.0 |
| requests | ≥ 2.31 | 2.34.2 |
| pyyaml | ≥ 6.0 | 6.0.3 |
| tqdm | ≥ 4.66 | 4.67.3 |
| jsonschema | ≥ 4.21 | 4.26.0 |

**External tools.** These are not Python packages; the recommended way to obtain them is the conda environment in [Installation](#installation).

| Tool | Required for | Version used in testing |
|---|---|---|
| C++20 compiler | building the bundled `get_featureIDs` helper | GCC ≥ 11 or Clang ≥ 13 (Apple Clang) |
| `bgzip`, `tabix` (htslib) | compressing and indexing BED output | 1.22.1 |
| seqtk | telomere detection (`scaffold`, `centromeres`, `karyotype`) | 1.5 |
| KMC | building databases (the bundled helper queries the resulting index; not needed to *use* a pre-built database) | 3.2.x (vendored API 3.2.4) |
| libcairo | rendering `--format pdf` / `--format png` | any recent release |
| samtools | only for BAM input to `annotate` | 1.22.1 |

**Hardware.** No non-standard hardware is required — KaryoScope runs on a standard CPU and has no GPU dependency. Resource needs scale with the input:

- **Demo and small inputs:** run on any laptop in seconds (see [Demo](#demo)).
- **Human whole-genome inputs (KMC backend):** the pre-built human database is ~17 GB on disk, and loading its KMC index during `annotate` needs roughly 30 GB of RAM (we measured ~30–35 GB peak at 16 threads). We recommend ≥ 50 GB RAM and ≥ 16 CPU cores for human-scale runs. A single human haplotype's six-feature-set run takes ~17–22 minutes at 16 threads.
- **Human whole-genome inputs (HKS backend):** the HKS database is ~24 GB on disk (~13 GB compressed). `annotate` holds the index (~10 GB, fixed) plus per-query buffering that grows with how much sequence one lookup processes at once, so the memory you should request depends on the input shape. A **single haplotype** peaks at **~10 GB** (request ≥ 16 GB) and finishes in **~7–9 minutes** at 16 threads. A **combined diploid assembly** — both haplotypes in one file, e.g. HG002 v1.1 — peaks at **~17 GB** (request ≥ 24 GB); annotating each haplotype separately keeps the peak at ~10 GB. Measured at 16 threads on T2T-CHM13v2.0, HG008N, HG008T, and HG002 v1.1.

## Installation

> Installation via Bioconda is planned. For now, install from source.

KaryoScope requires Python ≥3.10 and several external tools (`KMC`, `bgzip`, `tabix`, `seqtk`, `cairo` for PDF/PNG karyotype output, and `librsvg`/`rsvg-convert` for the SVG→PNG export used by `karyoplot` and `karyoscope-iscn zoom --png`). The simplest setup is a dedicated conda environment:

```bash
git clone https://github.com/barthel-lab/KaryoScope.git
cd KaryoScope

# Create the shared environment from environment.yml (Python, the
# bioinformatics tools, and the C/C++ + Rust toolchains for the query
# helpers). This is the same environment used across the KaryoScope
# ecosystem, so it works unchanged if you later add the sibling repos or ISCN.
conda env create -f environment.yml
conda activate karyoscope

# Install KaryoScope
pip install -e .

# Build the bundled C++ helper (`get_featureIDs`) for the KMC backend.
# `pip install` is Python-only and does NOT compile the C++ tree.
cd native/get_featureIDs && make && cd ../..
```

To use an **HKS** database (`index.type: hks`), also build the `hks` query
binary (the Rust toolchain is already in the environment). Clone HKS with its
submodules and install it onto `PATH`:

```bash
git clone --recurse-submodules https://github.com/jnalanko/HKS.git
cargo install --path HKS --root "$CONDA_PREFIX"   # installs $CONDA_PREFIX/bin/hks
```

The build produces `native/get_featureIDs/build/get_featureIDs`; the
Python wrapper finds it automatically. See [`native/README.md`](native/README.md)
for build-system details (CXX selection, `pkg-config`-driven zlib lookup,
and the macOS + conda `-Wl,-rpath,$CONDA_PREFIX/lib` shim).

**Typical install time.** On a normal desktop with a good network connection, creating the conda environment is the bulk of the time — usually ~5–10 minutes to resolve and download the bioinformatics tools. `pip install -e .` then takes under a minute, and building the C++ helper takes a few seconds. Reinstalls into an existing environment are much faster.

## Demo

This demo runs the full annotation step end to end in a few seconds, using a tiny synthetic database that ships with the repository (~5 MB unpacked). It downloads nothing and needs no special hardware, so it is the fastest way to confirm that your installation — the Python package **and** the compiled `get_featureIDs` helper — works. It is a smoke test on constructed inputs, not a biological example; the meaningful workflow on real data is in [Quick start](#quick-start) below.

> **Note.** The synthetic database is a real, structurally complete KaryoScope database (a genuine KMC index built from three short constructed sequences, plus the usual manifest, hierarchy, features, and colors files). Its k-mers are chosen so the annotation output is exactly predictable. It is the same fixture the integration test suite runs against. Running the demo does **not** require `kmc` to be installed — `kmc` is only needed to *build* a database; querying an existing one uses the bundled `get_featureIDs` helper.

Run it from the repository root with the `karyoscope` environment active:

```bash
# Use a throwaway database root so the demo never touches your real one.
export KARYOSCOPE_DB="$(mktemp -d)/db" && mkdir -p "$KARYOSCOPE_DB"

# 1. Install the bundled synthetic database and register it.
tar -xzf tests/data/dummy_db.tar.gz -C "$KARYOSCOPE_DB"
karyoscope register KS_dummy_test_v1

# 2. Build a tiny query: one sequence whose k-mers are in the database,
#    and one whose k-mers are not.
printf '>seq_with_features\nACGTGCTAGCTAGGCTATCGTAC\n>seq_novel\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n' > demo.fa

# 3. Annotate. --no-bgzip keeps the BEDs as plain text so you can read them.
karyoscope annotate -i demo.fa -o demo_out --no-bgzip
```

A convenience wrapper that runs exactly these steps (in a throwaway database root, cleaning up after itself) is provided at [`examples/run_demo.sh`](examples/run_demo.sh):

```bash
bash examples/run_demo.sh
```

**Expected run time:** a few seconds; the annotation step itself completes in under a second.

**Expected output.** Four BEDs are written under `demo_out/` (a presmoothed and a smoothed BED for each of the database's two feature sets). The smoothed chromosome track, `demo_out/demo.KS_dummy_test_v1.chromosome.smoothed.bed`, is:

```
seq_with_features	0	2	chr1
seq_with_features	2	3	chr2
seq_novel	0	10	novel
```

and the smoothed region track, `demo_out/demo.KS_dummy_test_v1.region.smoothed.bed`, is:

```
seq_with_features	0	1	rA
seq_with_features	1	2	rB
seq_with_features	2	3	rC
seq_novel	0	10	novel
```

The synthetic database maps three 21-mers to features (`chr1`/`rA`, `chr1`/`rB`, `chr2`/`rC`), so `seq_with_features` is annotated base by base; the poly-T `seq_novel` contains no indexed k-mers, so every position is labelled `novel`. This is the same end-to-end path that human-scale runs take — only the database and inputs differ.

## Quick start

This walkthrough uses the [HG002 v1.1 T2T diploid assembly](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz) as input, but any FASTA will work. Substitute your own with `--input <path>` throughout.

```bash
# 1. Download the recommended human reference database (~17 GB, one-time)
karyoscope download

# 2. Download the HG002 v1.1 diploid assembly (~3 GB, one-time)
#    Skip if you already have your own assembly to annotate.
curl -O https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz

# 3. Annotate the assembly. Recommended: at least 16 threads and, for the
#    HKS backend, ≥ 24 GB of RAM — HG002 v1.1 is a combined diploid
#    assembly and its annotate peaks at ~17 GB (a single haplotype peaks
#    at ~10 GB and fits 16 GB). The KMC backend needs ≥ 50 GB. HG002 runs
#    in ~20-30 min at -t 16.
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
| [`karyoscope build`](docs/commands/build.md) | Build an HKS database from a genome and per-feature-set BEDs |
| [`karyoscope annotate`](docs/commands/annotate.md) | Annotate sequences with *k*-mer features |
| [`karyoscope scaffold`](docs/commands/scaffold.md) | Order, orient, and rename assembly contigs |
| [`karyoscope remap-bed`](docs/commands/remap-bed.md) | Apply an existing scaffold map to a separately-annotated BED |
| [`karyoscope bin`](docs/commands/bin.md) | Aggregate base-pair annotations into larger bins |
| [`karyoscope centromeres`](docs/commands/centromeres.md) | Extract centromere coordinates |
| [`karyoscope karyotype`](docs/commands/karyotype.md) | Render karyotype visualization |
| [`karyoscope info`](docs/commands/info.md) | Inspect databases, files, installation |
| [`karyoscope version`](docs/commands/version.md) | Print version and environment info |

Run `karyoscope <command> --help` for full options on any command.

## Documentation

Per-command reference pages live under [`docs/commands/`](docs/commands/) — one page per subcommand, linked from the [Commands](#commands) table above. The `--help` output of each command (`karyoscope <command> --help`) remains the authoritative, always-current reference.

## Databases

KaryoScope works with pre-built databases distributed via the [KaryoScope registry](https://github.com/barthel-lab/KaryoScope-registry). The current default is `KS_human_CHM13_v2` (~17 GB), built from the T2T-CHM13v2.0 reference.

Browse and download available databases:

```bash
karyoscope download --list
```

You can also build your own database from a genome and per-feature-set BED annotations with [`karyoscope build`](docs/commands/build.md).

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

To cite the software itself — for example, to record the exact version you ran — use its Zenodo archive. The concept DOI [10.5281/zenodo.20657816](https://doi.org/10.5281/zenodo.20657816) always resolves to the latest release, and each release also receives its own version-specific DOI on the same record.

A `CITATION.cff` file in this repository provides machine-readable citation metadata.

## License

KaryoScope is licensed under [GPL-3.0-or-later](LICENSE) due to its dependency on the GPL-3.0 KMC library. A future release will switch to MIT once we migrate to [HKS](https://github.com/jnalanko/HKS) for *k*-mer indexing.

## Contributing

Contributions, bug reports, and feature requests are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) to get started, and our [Code of Conduct](CODE_OF_CONDUCT.md) for community norms.

## Acknowledgments

Developed in the [Barthel Lab](https://www.barthel-lab.com/) at the Translational Genomics Research Institute (TGen), in collaboration with Jarno Alanko, Simon Puglisi, and Camille Marchet.
