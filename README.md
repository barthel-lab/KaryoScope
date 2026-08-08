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

> 🚀 **New in v2.2.0.** [`karyoscope prep-bed`](docs/commands/prep-bed.md) converts source annotations — RepeatMasker, EDTA, GFF3/GTF gene models, UCSC cytoband tables, CenSat, satellite catalogs, a plain `.fai` — into the feature-set BEDs [`karyoscope build`](docs/commands/build.md) consumes, so building a database no longer needs a script per dataset. [`docs/recipes/`](docs/recipes/) rebuilds the shipped human and *Arabidopsis* databases end to end, from download URL to built index. [`examples/karyotypes/`](examples/karyotypes/) collects reference karyotype plots for six assemblies to compare your own output against. Legends can now be collapsed with a `legend_group` column. See the [CHANGELOG](CHANGELOG.md).

---

## Overview

KaryoScope is an alignment-free annotation tool that assigns each *k*-mer in a query assembly or sequencing read to a feature drawn from one or more user-defined hierarchical feature sets, producing a base-pair resolution annotation in a single pass. Because a feature set is any tiling of a reference with labelled regions, KaryoScope is extensible to arbitrary annotation sources, from satellite catalogs and repeat libraries to cytobands, FISH-probe coordinates, and structural-variant breakpoints.

A pre-built database for the human genome is distributed alongside the tool, derived from T2T-CHM13v2.0 with six feature sets covering chromosome of origin, satellite composition, interspersed repeats, subtelomeric structure, gene boundaries, and acrocentric-specific features. From these annotations, KaryoScope produces karyotype visualizations and cytogenetic reports without ever performing read alignment. Additional databases can be built for any reference genome or community-curated annotation source. [`karyoscope prep-bed`](docs/commands/prep-bed.md) converts the annotation into labelled feature-set BEDs and [`karyoscope build`](docs/commands/build.md) indexes them; [`docs/recipes/`](docs/recipes/) works both steps through for the human and *Arabidopsis* databases. Databases built this way can be published to, and installed from, the [KaryoScope registry](https://github.com/barthel-lab/KaryoScope-registry).

<!-- TODO: hero figure of a KaryoScope output, e.g. the HG008T karyotype -->
<!-- <p align="center"><img src="assets/hero_karyotype.png" alt="Example KaryoScope karyotype" width="800"/></p> -->

### Why alignment-free?

- **Pangenome-scale throughput.** On the [HKS](https://github.com/jnalanko/HKS) *k*-mer indexing backend, one `annotate` run covering all six feature sets finishes a complete human haplotype in **~7–9 minutes at a ~10 GB memory peak**, at 16 threads — scaling to cohorts of hundreds of phased assemblies. The original KMC backend does the same work in ~13–14 minutes at a ~36 GB peak, so HKS is roughly 1.5–1.8× faster on under a third of the RAM. (Measured at 16 threads on T2T-CHM13v2.0 and HG002 v1.1 hap1, as a single sequential run of the whole pipeline — assignment plus smoothing for every feature set. See [System requirements](#system-requirements) for the machine.)
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
| [`hks`](https://github.com/jnalanko/HKS) | building databases with `karyoscope build`, and annotating against HKS-backend databases (e.g. `HKS_human_CHM13_v2`) | **0.4.0 or newer** |
| KMC | building a *KMC-backend* database (legacy — `karyoscope build` produces HKS databases). Not needed to *use* a pre-built KMC database; the bundled `get_featureIDs` helper queries its index directly | 3.2.x (vendored API 3.2.4) |
| libcairo | rendering `--format pdf` / `--format png` | any recent release |
| samtools | BAM input to `annotate`, and `karyoscope build` when the genome has no `.fai` index yet | 1.22.1 |

**Hardware.** No non-standard hardware is required — KaryoScope runs on a standard CPU and has no GPU dependency. Every runtime and memory figure in this README was measured on a cluster node with two Intel Xeon Gold 6342 CPUs (Ice Lake, 2.80 GHz, 24 cores each, SMT disabled) and 503 GiB RAM, running the job on 16 threads. Resource needs scale with the input:

- **Demo and small inputs:** run on any laptop in seconds (see [Demo](#demo)).
- **Human whole-genome inputs (HKS backend, recommended):** `annotate`'s peak is set by the index, not by the input. It holds the shared base index (6.0 GiB) plus one feature set's labeling at a time (the largest is 3.0 GiB), so the peak is **~10 GB whatever you annotate** — a single haplotype and a combined diploid assembly measure the same. **Request ≥ 16 GB.** A single haplotype's six-feature-set run takes **~7–9 minutes** at 16 threads; a combined diploid is roughly twice that (~14–16 minutes). Measured on T2T-CHM13v2.0 and HG002 v1.1 (hap1 and combined diploid). See [Disk space](#disk-space) for storage.

  Memory barely moves with `--threads`, so there is no reason to hold threads back to save RAM. What *does* vary by machine is the index load: 9 GB is read per feature set, and if there is enough free RAM for the page cache to hold it, the second and later feature sets load in ~2 s instead of ~25 s. On a memory-tight machine expect every load to be cold.
- **Human whole-genome inputs (KMC backend):** `annotate` peaks at ~36 GB on a single haplotype and ~46 GB on a combined diploid, so we recommend ≥ 96 GB RAM and ≥ 16 CPU cores. A single haplotype's six-feature-set run takes ~13–14 minutes at 16 threads; a combined diploid ~23–27 minutes.

### Disk space

The archive is **not** deleted until extraction finishes, so installing needs room for the download *and* the extracted database at the same time. That peak, not the download size, is the number to check against `df -h`:

| Database | Backend | Download (`.tar.gz`) | On disk after install | **Free space to install** |
|---|---|---|---|---|
| `HKS_human_CHM13_v2` | HKS | 13.3 GB | 22.7 GB | **~36 GB** |
| `KS_human_CHM13_v2` (default) | KMC | 16.3 GB | 17.2 GB | **~34 GB** |
| `HKS_arabidopsis_ColCEN` | HKS | 0.46 GB | 0.68 GB | **~1.1 GB** |

Size tracks the genome, not the tool: the 135 Mb *Arabidopsis* database installs in about a gigabyte, while the human references need tens. Between the two human databases, the HKS archive compresses far better than the KMC one, so its download is the *smaller* while its installed footprint is the larger. Once the install completes the archive is removed and only the "on disk" column remains occupied.

`karyoscope download` checks free space before it starts transferring anything and refuses with the shortfall spelled out rather than failing after a 25-minute download. `karyoscope download --info <ID>` prints all three figures for any database. Sizes here are decimal GB (10⁹ bytes); `df -h` reports binary GiB, so 36 GB shows there as 34 GiB.

**Annotation output** is also substantial, and `--bgzip` does not reduce the peak — every BED is written in full before the compression pass runs. Budget roughly **0.8 GB per feature set per Gbp of input**:

| Input | Feature sets | Output BEDs (uncompressed) | Peak (HKS backend) | Peak (KMC backend) |
|---|---|---|---|---|
| Single human haplotype (~3.1 Gbp) | 6 | ~15 GB | ~15 GB | ~18 GB |
| Diploid assembly, e.g. HG002 v1.1 (~6.0 Gbp) | 6 | ~29 GB | ~29 GB | ~34 GB |

Measured on HG002 v1.1 against `HKS_human_CHM13_v2`. An HKS run needs no headroom beyond its outputs: `hks lookup` writes the presmoothed BED directly, so no intermediate copy of it ever exists. (Under `--no-keep-presmoothed` that file becomes a temp one, and the KMC column applies instead.) The KMC backend still writes a combined intermediate BED alongside.

Restricting `--feature-set`, or dropping one output variant with `--no-keep-presmoothed` / `--no-smooth`, scales this down proportionally. `annotate` estimates the footprint from the input size and checks `--outdir` before starting; pass `--no-resource-check` to override the estimate.

## Installation

> Installation via Bioconda is planned. For now, install from source.

KaryoScope requires Python ≥3.10 and several external tools (`bgzip`, `tabix`, `seqtk`, and `cairo` for PDF/PNG karyotype output). It also needs a *k*-mer backend query helper — see [k-mer index backends](#k-mer-index-backends) below. The simplest setup is a dedicated conda environment:

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
```

### k-mer index backends

KaryoScope queries a database through one of two interchangeable *k*-mer backends;
each database declares which one it uses in its manifest (`index.type`). Install the
backend(s) your databases need — most users need only one. `pip install` is
Python-only and builds neither, so each backend's query helper is a separate step.

#### HKS

[HKS](https://github.com/jnalanko/HKS) powers the `HKS_human_CHM13_v2` database, any
other `index.type: hks` database, and databases you create with
[`karyoscope build`](docs/commands/build.md). It annotates a human haplotype
~1.5–1.8× faster than KMC at under a third of the memory. The Rust toolchain is already
in the environment, so clone HKS with its submodules and install the `hks` binary
onto `PATH`:

```bash
git clone --recurse-submodules https://github.com/jnalanko/HKS.git
cargo install --path HKS --root "$CONDA_PREFIX"   # installs $CONDA_PREFIX/bin/hks
```

KaryoScope finds `hks` on `PATH` automatically (or set `$KARYOSCOPE_HKS` to its
path). A future release will bundle `hks` via conda so this step goes away.

**`hks` 0.4.0 or newer is required.** KaryoScope asks `hks` to write its output
in the exact shape it needs — headerless, with `novel` for k-mers absent from
the index — which removes a full read and a full write of a multi-gigabyte file
per feature set. The options that do this (`--miss-label`, `--no-header`,
`--report-label-ids` on `smooth`) landed in 0.3.0. KaryoScope reads the version
`hks` logs on startup and refuses to begin against an older one, rather than
failing partway through with an unrecognised-argument error.

#### KMC

KMC powers the `KS_human_CHM13_v2` database — the current default from
`karyoscope download`. Its query helper is a bundled C++ tool, `get_featureIDs`,
that you compile in place:

```bash
cd native/get_featureIDs && make && cd ../..
```

This produces `native/get_featureIDs/build/get_featureIDs`; the Python wrapper finds
it automatically. See [`native/README.md`](native/README.md) for build-system details
(CXX selection, `pkg-config`-driven zlib lookup, and the macOS + conda
`-Wl,-rpath,$CONDA_PREFIX/lib` shim). The `KMC` tool itself is only needed to *build*
a KMC database, not to query one.

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

<!-- TODO(hks-default): this Quick start uses the default `karyoscope download`,
     which currently installs the KMC KS_human_CHM13_v2 database (hence the KMC
     name in the output filenames below and the KMC RAM figures). When
     the registry default flips to HKS_human_CHM13_v2 (after KMC removal + the
     conda recipe guarantees `hks` on every install), switch this walkthrough to
     HKS: the output names become HKS_human_CHM13_v2 and the RAM guidance drops
     to the HKS figure (~16 GB, since the peak is the index and not the
     input). Keep it on KMC until then. -->

```bash
# 1. Download the recommended human reference database (one-time).
#    16.3 GB download, 17.2 GB on disk once installed -- but you need
#    ~34 GB free to install, because the archive isn't deleted until
#    extraction finishes. See "Disk space" above.
karyoscope download

# 2. Download the HG002 v1.1 diploid assembly (~3 GB, one-time)
#    Skip if you already have your own assembly to annotate.
curl -O https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz

# 3. Annotate the assembly. This runs against the KMC database installed
#    in step 1, so budget for the KMC backend: it peaks at ~46 GB on
#    this diploid input, so request >= 96 GB of RAM, with at least 16
#    threads. HG002 v1.1 is a combined diploid and takes ~23-27 min at
#    -t 16; a single haplotype takes ~13-14 min.
#    Allow ~34 GB free in results/ for six feature sets over this 6 Gbp
#    assembly, plus headroom for the combined-BED intermediate the KMC
#    backend writes alongside them. See "Disk space" above.
#    --no-bgzip keeps the per-feature-set BEDs as plain
#    text for easy inspection; drop it to get the default bgzipped
#    outputs (which shrink the result but not the peak, since compression
#    runs after every BED has been written).
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
| [`karyoscope prep-bed`](docs/commands/prep-bed.md) | Convert a source annotation into a feature-set BED |
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

Beyond the per-command pages:

| Where | What |
|---|---|
| [`docs/recipes/`](docs/recipes/) | Rebuild the shipped human and *Arabidopsis* databases from published sources, and build an HPV database from [PaVE](https://pave.niaid.nih.gov/) reference genomes — every download URL and checksum, the `prep-bed` command for each feature set, and the build spec. |
| [`examples/karyotypes/`](examples/karyotypes/) | Reference karyotype plots for six assemblies (CHM13, HG002, an HG008 tumour/normal pair, an HPRC population sample, and *Arabidopsis*), with notes on what each shows and the commands that produced them. |

## Databases

KaryoScope works with pre-built databases distributed via the [KaryoScope registry](https://github.com/barthel-lab/KaryoScope-registry). The current default is `KS_human_CHM13_v2` (17.2 GB installed; see [Disk space](#disk-space) for the free space needed to install it), built from the T2T-CHM13v2.0 reference.

The registry also carries `HKS_arabidopsis_ColCEN` (0.68 GB installed), built from the *Arabidopsis thaliana* Col-CEN v1.2 telomere-to-telomere assembly, with `chromosome`, `region`, `repeat` and `gene` feature sets. Two differences from the human databases: its structural set is `region` rather than `subtelomeric`, and rendering a karyotype from it needs `--telo-motif CCCTAAA`, because plant telomeres are `TTTAGGG` and the default motif finds none.

Browse and download available databases:

```bash
karyoscope download --list           # download and on-disk size for each
karyoscope download --info <ID>      # ...plus the free space needed to install
```

You can also build your own database from a genome and per-feature-set BED annotations with [`karyoscope build`](docs/commands/build.md). `build` starts from a final labelled BED; [`karyoscope prep-bed`](docs/commands/prep-bed.md) produces one from the formats annotation usually arrives in — RepeatMasker and EDTA tables, GFF3/GTF gene models, UCSC cytoband tables, CenSat and satellite catalogs, PaVE genome records, and a plain `.fai` — and prints the build-spec stanza to go with it.

For complete worked examples — every download URL, checksum, conversion and build spec — see [docs/recipes/](docs/recipes/), which rebuilds the human CHM13v2 and *Arabidopsis* Col-CEN databases from published sources and builds an HPV database from PaVE reference genomes.

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

KaryoScope is licensed under [GPL-3.0-or-later](LICENSE) because the KMC backend links the GPL-3.0 KMC library. The [HKS](https://github.com/jnalanko/HKS) backend (MIT) is now available alongside KMC; a future release will remove the KMC dependency and relicense to MIT.

## Contributing

Contributions, bug reports, and feature requests are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) to get started, and our [Code of Conduct](CODE_OF_CONDUCT.md) for community norms.

## Acknowledgments

Developed in the [Barthel Lab](https://www.barthel-lab.com/) at the Translational Genomics Research Institute (TGen), in collaboration with Jarno Alanko, Simon Puglisi, and Camille Marchet.
