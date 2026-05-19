# `native/` — KaryoScope C++ helpers

This directory holds the C++ code that KaryoScope uses for its
performance-critical k-mer querying. The Python package (`src/karyoscope/`)
invokes the resulting binaries via subprocess; there is no Python C
extension and no Python-side compilation step.

## Layout

```
native/
├── README.md                       # this file
├── external/                       # vendored third-party code (read-only)
│   ├── kmc_api/                    # KMC 3.2.4 API headers + sources (GPL-3)
│   └── cxxopts/cxxopts.hpp         # cxxopts 3.2.0 single-header (MIT)
└── get_featureIDs/                 # one binary per subdirectory
    ├── Makefile
    ├── get_featureIDs.cpp
    └── build/                      # gitignored, created by `make`
        └── get_featureIDs
```

## Building

Prerequisites:

- A C++20-capable compiler (gcc 10+, clang 13+, Apple Clang 15+).
- `zlib` development headers and library. On Linux, `apt install
  zlib1g-dev`; on macOS, the Xcode Command Line Tools include it; in a
  conda environment, `conda install zlib`.
- `make`.

To build:

```bash
make -C native/get_featureIDs           # release build, -O3
make -C native/get_featureIDs -j4       # parallel
make -C native/get_featureIDs clean     # remove build/
```

For debug builds: `make -C native/get_featureIDs OPTIMIZE='-O0 -g'`.

The Makefile honours `CXX` to pick a non-default compiler (e.g.
`CXX=clang++` or `CXX=g++-13`). zlib is located via `pkg-config` when
available, with a fallback to plain `-lz`. If `CONDA_PREFIX` is set in
the environment, its `include/` and `lib/` directories are added to the
search paths so conda installs of zlib work without further
configuration.

## Vendored code

### `external/kmc_api/`

Sourced verbatim from [KMC](https://github.com/refresh-bio/KMC) version
3.2.4, with **one minor patch** to support macOS:

- `kmer_defs.h`: the unconditional `#include <ext/algorithm>` is wrapped
  in `#ifdef __GLIBCXX__` so it's skipped when building with libc++
  (Apple Clang's default standard library on macOS). The header is a
  libstdc++ extension and doesn't exist in libc++; KMC doesn't actually
  use any symbols from it in the files we vendor, so this is safe.

The KMC API is licensed under GPL-3, which is why KaryoScope is
GPL-3-or-later overall. A planned migration to
[HKS](https://github.com/jnalanko/HKS) will let us relicense to MIT.

### `external/cxxopts/cxxopts.hpp`

[cxxopts](https://github.com/jarro2783/cxxopts) version 3.2.0, single
header, MIT-licensed. Used by `get_featureIDs.cpp` for command-line
parsing.

## How the Python side finds the binary

`src/karyoscope/core/io/kmc.py` resolves the binary in this order:

1. `$KARYOSCOPE_GET_FEATUREIDS` — explicit override (most useful for
   tests or custom installs).
2. `get_featureIDs` on `$PATH` — the path Bioconda installs will take.
3. `<repo>/native/get_featureIDs/build/get_featureIDs` — found by
   walking up from the package source location, which works for
   editable installs (`pip install -e .`).

If none resolve, the user gets an actionable error pointing at this
README's "Building" section.

## Why not Bioconda yet?

Decided in the Stage 5 design discussion: Bioconda is the right
distribution channel for v1.0 (it's where the rest of our users'
toolchain — KMC, bgzip, samtools — comes from), but setting up the
recipe is a polish step rather than something that blocks core feature
work. For Stage 5a we just need the source-build path, which is what
developers and CI need anyway.

## Adding another binary

If you add a second program (e.g., `native/build_kmc_index/`), the
pattern is:

1. Create the subdirectory with its own `Makefile` modelled on this one.
2. Put the source under that directory; keep shared dependencies in
   `native/external/`.
3. Add a Python wrapper in `src/karyoscope/core/io/` and an integration
   test marked `@pytest.mark.integration`.
4. Add a build step to `.github/workflows/ci.yml`'s
   `cpp_build_and_integration` job.

The `.gitignore` already ignores `native/*/build/`, so no changes there.

## Troubleshooting

### macOS + conda: `dyld: Library not loaded: @rpath/libz.1.dylib`

The full error looks like:

```
dyld[...]: Library not loaded: @rpath/libz.1.dylib
  Referenced from: <...>/build/get_featureIDs
  Reason: no LC_RPATH's found
```

This happens when the build picks up zlib from a conda environment (via
`$CONDA_PREFIX/lib`) but the resulting binary has no `LC_RPATH` entry
telling dyld where to find dylibs at runtime. The Makefile in this repo
sets `-Wl,-rpath,$(CONDA_PREFIX)/lib` whenever `CONDA_PREFIX` is set, so
this should be handled automatically. If you hit this anyway:

- Confirm `CONDA_PREFIX` was set when you ran `make`. The Makefile only
  adds the rpath if it sees it. `echo $CONDA_PREFIX` should print your
  active env's path.
- If you're in a conda env that doesn't set `CONDA_PREFIX` (some manual
  activations don't), set it explicitly: `make CONDA_PREFIX=$(conda
  info --base)/envs/your-env`.
- Inspect the binary: `otool -l build/get_featureIDs | grep -A 2
  LC_RPATH` should list your env's `lib` directory. If it doesn't,
  re-run `make clean && make` from inside the activated env.

Background: conda-forge dylibs on macOS use the relocatable
`@rpath/libfoo.dylib` install-name pattern. The consuming binary has to
embed an `LC_RPATH` entry to tell dyld where to search. Conda's own
compilers do this automatically; the system `c++` doesn't, so we do it
explicitly. On Linux the same flag adds `DT_RUNPATH`, which is benign.
