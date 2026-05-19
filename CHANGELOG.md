# Changelog

All notable changes to KaryoScope will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial project scaffold with command-line interface skeleton.
- `karyoscope` console entry point with subcommand dispatch.
- Foundational documentation: `README.md`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`,
  and issue/PR templates.
- GPL-3.0 license.
- Python packaging via `hatchling`, requiring Python ≥3.10.
- GitHub Actions CI workflow running lint, format, and tests on every push.
- Integration with the `KaryoScope-registry` for pre-built database discovery.
- Full implementation of the `karyoscope download` command, supporting installing
  databases by id or as the registry default, listing available databases
  (`--list`) with `--organism` / `--tag` / `--community` filters, inspecting
  individual entries (`--info ID`), showing locally installed databases
  (`--status`), and uninstalling them (`--remove ID`).
- Registry fetching with a 24-hour TTL cache at `<db_root>/registry_cache.yaml`
  and graceful fall-back to a stale cache on transient network errors. The
  default registry URL points at `barthel-lab/KaryoScope-registry`; private or
  test registries can be supplied via `--registry-url`.
- SHA-256 verification of downloaded tarballs, with `--no-checksum` for
  debugging. HTTP downloads stage to a `.part` file and resume via `Range`
  requests when interrupted.
- Safe tarball extraction that refuses entries containing `..`, symlinks or
  other special types, or files outside the database's expected top-level
  directory.
- Per-database manifest schema (`manifest.yaml`) and on-disk layout validator.
- `installed.json` for tracking which databases are installed locally, with
  atomic writes and recovery from corrupted files.
- Test fixtures: a 22 KB dummy database tarball containing a real (tiny) KMC
  index, plus a build script for regenerating it.
- Top-level verbosity flags ``-v`` / ``-vv`` (raise log level to INFO / DEBUG)
  and ``-q`` (lower to ERROR-only). Logging is reserved for diagnostics
  (URL fetches, cache hits, SHA-256 verifications, subprocess invocations);
  program output remains visible regardless of verbosity. Default is WARNING.
- ``karyoscope.core.external``: small wrapper around :mod:`subprocess`
  with consistent error reporting (:class:`ExternalToolError` includes the
  command, exit code, and a tail of stderr) and a :func:`require_tool`
  helper for "binary not found on PATH" cases with actionable error messages.
- ``karyoscope.core.io.hierarchy``: parser for the ``hierarchy.tsv``
  format used inside KaryoScope databases.
- The ``karyoscope info`` command now provides a real implementation
  rather than a stub. With no argument it lists installed databases
  with version, install date, size, and source. Given a database id, it
  prints the parsed manifest plus feature-set counts derived from the
  database's ``hierarchy.tsv``. Given a filesystem path, it probes
  whether the path is a KaryoScope database directory.
- `native/` directory holding the C++ helper code that KaryoScope uses
  for the performance-critical k-mer querying step. Contains a vendored
  copy of the [KMC](https://github.com/refresh-bio/KMC) 3.2.4 API
  (`native/external/kmc_api/`, GPL-3) with a one-line patch wrapping the
  libstdc++-only `<ext/algorithm>` include in `#ifdef __GLIBCXX__` so
  the source compiles under macOS libc++. Also vendors
  [cxxopts](https://github.com/jarro2783/cxxopts) 3.2.0 as a single
  header (MIT-licensed).
- `native/get_featureIDs/` — the first C++ binary: a multi-threaded
  FASTA/FASTQ reader that queries every k-mer in each input sequence
  against a KMC database and emits run-length-encoded BED records of
  per-position feature ids. Supports plain FASTA/FASTQ, gzipped input,
  and stdin. Adapted from the version in the KaryoScope archive repo;
  the only changes are include paths driven by the new Makefile's `-I`
  flags rather than the previous `../external/...` relative paths.
- `native/get_featureIDs/Makefile` — cross-platform build with
  incremental compilation. Defaults `CXX` to `c++` (the system's
  default C++ compiler — g++ on Linux, Apple Clang on macOS). Locates
  zlib via `pkg-config` with a fall-back to plain `-lz`, and picks up
  `$CONDA_PREFIX` paths when set. Requires C++20.
- `karyoscope.core.io.kmc` — Python wrapper around `get_featureIDs`.
  Resolves the binary in precedence order:
  `$KARYOSCOPE_GET_FEATUREIDS` → `shutil.which("get_featureIDs")` →
  walking up from the package source to find
  `<repo>/native/get_featureIDs/build/get_featureIDs` (works in
  editable installs). Raises `ToolNotFoundError` with actionable
  instructions if none resolve.
- Integration tests (`tests/test_kmc.py`, marked
  `@pytest.mark.integration`) verifying the full pipeline:
  Python wrapper → subprocess → BED output → schema validation. CI
  builds the binary on Linux and macOS in a separate
  `cpp_build_and_integration` job and runs `pytest -m integration`.
- `native/README.md` documenting the layout, vendored code, build
  prerequisites, and the binary-lookup logic.
- The `karyoscope annotate` command, replacing the previous stub. Takes
  a FASTA (plain or gzipped) and produces one BED per feature set,
  named `<input>.<dbid>.<feature_set>.presmoothed.bed[.gz]`. Each
  record's 4th column is the human-readable feature name in that set,
  obtained by translating ``get_featureIDs``' integer ids through
  ``features.tsv``. Feature id 0 (k-mer not in the KMC index) renders
  as ``"novel"``; any other id absent from features.tsv is a hard
  error rather than a silent fallback (real KaryoScope feature sets
  can include ``"Unknown"`` as a legitimate name, so mapping missing
  ids to a placeholder string would be ambiguous and would hide
  database / index mismatches). Adjacent records that translate to
  the same feature name in a given set are merged in a streaming pass
  — important because two different ids can collapse to the same name
  (e.g., ids for ``(chr1, rA)`` and ``(chr1, rB)`` both have ``chr1``
  in the chromosome BED).
- Flags: ``-i/--input`` (required), ``-o/--outdir``, ``--db ID``,
  ``--feature-set NAME`` (repeatable; default: all sets in the
  manifest), ``-t/--threads``, ``--keep-intermediates`` (default: off
  — the combined BED from the C++ step is deleted after splitting),
  and ``--bgzip/--no-bgzip`` (default: ``--bgzip``).
- ``karyoscope.core.io.features`` — parser for ``features.tsv`` in the
  one-row-per-id, one-column-per-feature-set schema. Validates header,
  rejects featureID 0 in the file (reserved sentinel) and duplicate ids.
  Exposes :func:`render_feature` which raises :class:`FeaturesError`
  on positive ids absent from the table.
- ``karyoscope.core.annotate`` — orchestration module for the annotate
  pipeline (resolve db → validate layout → run get_featureIDs → split
  combined BED → bgzip). Importable for callers other than the CLI.
- Integration tests for the annotate command exercising the default
  path, ``--no-bgzip``, ``--keep-intermediates``, feature-set
  filtering, and error paths (unknown set, no database installed,
  featureID present in BED but absent from features.tsv).
  Marked ``@pytest.mark.integration``.
- Hierarchy-aware smoothing of per-feature-set BED tracks: the
  ``karyoscope annotate`` command now produces a second BED per
  feature set (``<input>.<dbid>.<feature_set>.smoothed.bed[.gz]``)
  alongside the existing ``.presmoothed`` BED. Smoothing promotes
  noisy intermediate intervals (especially short ``novel`` runs
  flanked by specific features) to the lowest common ancestor of
  their flankers in the feature set's hierarchy, generalising from
  "I don't know exactly what this is" to "I know it's at least
  this".
- Two new flags on ``annotate``: ``--smooth/--no-smooth`` (default
  ``--smooth``) and ``--keep-presmoothed/--no-keep-presmoothed``
  (default ``--keep-presmoothed``). The combination
  ``--no-smooth --no-keep-presmoothed`` would produce no output and
  is rejected with a clean error.
- ``karyoscope.core.smooth`` — a faithful port of the archive's
  ``smooth_features.py``. Public API: :class:`HierarchyIndex`
  (ancestor / LCA queries memoised per feature set),
  :func:`smooth_intervals` (the LCA-promotion + fixed-point algorithm,
  with ``max_gap = 1000``), :func:`merge_adjacent` (interval
  coalescing that preserves the novelty boundary at the root),
  :func:`process_seq_chunk` (per-worker entry point for
  ``multiprocessing.Pool``), and :func:`chunked_seq_reader` (yields
  chunks of BED lines flushed only at sequence boundaries, so
  workers always see complete sequences). Workers are initialised
  via :func:`worker_initializer` with the shared :class:`HierarchyIndex`
  and id-to-name table. The archive's ``_specific``-suffix stripping
  has been removed (it was a leftover from an older implementation),
  and the ``Phylogeny`` class has been renamed to :class:`HierarchyIndex`
  to reflect the more general terminology.
- Sequence chunking respects sequence boundaries: ``chunked_seq_reader``
  holds back the current sequence's lines until the chunk has at
  least ``chunk_size`` lines AND the next line is from a new
  sequence. This guarantees that :func:`smooth_intervals` always
  sees the complete flanking context for each sequence.
- ``karyoscope.core.io.hierarchy`` rewritten for the production
  schema: three columns (``feature_set``, ``child``, ``parent``)
  with each feature set forming an independent rooted tree.
  Required root name is ``"categorized"`` for v0.1 (documented
  restriction). Each feature set's rows are interpreted only against
  that feature set's own rows — parent names must refer to nodes
  within the same set, not cross-set ancestors.
- :func:`karyoscope.core.io.hierarchy.validate_hierarchy` — a
  validation pass separate from the structural parser. Per feature
  set independently, it checks (1) no child has multiple parent
  rows, (2) exactly one root exists and is named ``"categorized"``,
  (3) parent-pointer walks terminate (no cycles, all connected),
  and (4) every name in features.tsv's columns has a corresponding
  hierarchy node (optional, opt-in via ``feature_columns=``).
  Returns a list of human-readable issue strings.
- ``karyoscope info`` runs hierarchy validation when inspecting a
  database and prints any issues as warnings (does not exit
  non-zero). When ``annotate`` (or any command that actually
  consumes the hierarchy) encounters the same issues, it raises a
  hard error instead.
- Multiprocessing: smoothing uses ``multiprocessing.Pool`` sized by
  the ``-t/--threads`` flag (same flag that's passed to the C++
  k-mer-query step). The pool initializer pickles the
  :class:`HierarchyIndex` and feature table once per worker rather
  than once per chunk.
- Comprehensive unit tests for the smoothing algorithm
  (``tests/test_smooth.py``): :class:`HierarchyIndex` correctness
  (ancestors, LCA on siblings / across subtrees / for the root /
  for unknown nodes, caching, malformed-hierarchy rejection),
  :func:`merge_adjacent` cases (same-feature coalescing, gap
  preservation, root + novelty boundary, non-root novelty ignored),
  :func:`smooth_intervals` (promotion to sibling LCA, deeper LCA
  promotion, no-op when flankers are missing, ``max_gap``
  enforcement, no-demotion guarantee, empty input), and
  :func:`chunked_seq_reader` (flushing only at sequence
  boundaries).
- Integration tests for the smoothed-output paths in
  ``tests/test_command_annotate.py``: default flags produce both
  variants; ``--no-smooth`` produces only presmoothed;
  ``--no-keep-presmoothed`` produces only smoothed;
  ``--no-smooth --no-keep-presmoothed`` errors cleanly; end-to-end
  smoothing of a constructed ``rA → novel → rB`` query produces a
  region BED where the novel run has been promoted to ``aSat`` (the
  LCA of ``rA`` and ``rB`` in the dummy db's hierarchy).
- ``karyoscope bin`` command, replacing the previous stub. Aggregates
  a sorted per-base-pair (or run-length-encoded) BED into
  fixed-bin-size windows, labelling each bin with the feature that
  wins a deterministic three-rule selection: (1) sum bp of overlap per
  feature within the bin, (2) if a leaf-feature set is in play, leaf
  features compete first (fall back to all features when no leaf is
  present), (3) the ``novel`` sentinel only wins if it covers a strict
  majority of the bin — otherwise the best non-novel feature wins.
  Adjacent bins with the same winning feature are coalesced in a
  streaming pass, so output is coordinate-sorted and free of
  same-label neighbours.
- ``karyoscope.core.bin`` — the binner as a library. Public API:
  :func:`bin_records` (pure in-memory binning, useful for tests and
  in-process callers), :func:`bin_features` (file-to-file with
  transparent ``.gz`` in/out and ``-`` for stdin/stdout), and
  :func:`leaves_for` (compute the leaf set of a feature set from a
  parsed :class:`Hierarchy`). Importable in-process by
  :mod:`karyoscope.commands.scaffold` (next stage), which avoids
  shelling out to bin its own intermediate BEDs.
- ``BinError(KaryoscopeError)`` for malformed-input / malformed-flag
  conditions surfaced by the binner.
- ``bin`` CLI flags: ``-i/--input``, ``-o/--output``, ``-b/--bin-size``
  (required), and optionally ``--db ID`` plus ``--feature-set NAME``
  for leaf prioritisation. Passing ``--db`` without ``--feature-set``
  is an error (the leaf set is per-feature-set; one without the other
  is ambiguous). Bare invocation (no ``--db`` / ``--feature-set``) is
  supported and skips leaf prioritisation entirely.
- :func:`karyoscope.core.annotate.resolve_database` promoted from the
  underscore-private ``_resolve_database`` so the same db-id
  resolution logic is available to :mod:`karyoscope.commands.bin_cmd`
  and (next stage) :mod:`karyoscope.commands.scaffold` without
  duplicating it.
- Unit + integration tests for the binner in ``tests/test_bin.py``:
  the selection-rule ladder (overlap-largest-wins, novel-majority,
  novel-minority-falls-back, equal-half-novel-loses, alphabetic
  tie-break, all-novel-keeps-novel), leaf-prioritisation
  (leaf-beats-internal, fall-back-when-no-leaf-present,
  competition-between-leaves), pure binning (single record, coalescing,
  no-coalesce-across-chroms, trailing-bin truncation), file I/O (plain
  and gzip in/out, malformed-row error), and CLI invocation against
  the dummy database for the leaf-prioritisation path.
- Deliberate diff from the archive's ``bin_features.py``: the
  ``--specific`` suffix path is dropped. v0.1's ``hierarchy.tsv``
  makes leaf detection structural (a child that is never a parent),
  so the old "feature name ending in ``_specific``" hack is no
  longer needed. Callers that want leaf prioritisation pass an
  explicit ``leaf_set`` (in-process) or ``--feature-set`` (CLI).
- ``karyoscope scaffold`` command, replacing the previous stub.
  Takes one or more FASTA inputs (one per haplotype typically) and
  produces per-input scaffolded outputs: a ``scaffold_map.tsv`` (the
  authoritative source-of-truth mapping from encoded contig name back
  to source), a ``scaffold_stats.tsv`` (legacy 2-column format kept
  for back-compat with archive scripts), and one rewritten BED per
  feature set
  (``<input_stem>.<dbid>.<feature_set>.smoothed.scaffolded.bed[.gz]``).
  The encoded contig-name format is always
  ``<chrom>_<hap>_<contig>[_rc]``; the map file is the contract that
  downstream stages parse, so the encoded name can change between
  releases without breaking the pipeline.
- Topology-preservation principle: each ``-i`` input produces its own
  set of outputs (one map file, one stats file, N scaffolded BEDs).
  The Snakemake pipeline always collapsed everything into one combined
  per-sample BED; the CLI does not. Joint information that must be
  considered across files (orientation, category ordering, chromosome
  cell grouping) is computed in-memory across all inputs at once,
  but the artefacts on disk mirror what came in. This works cleanly
  for the pangenome convention (one file per hap), the HG002
  distributed-as-one-file convention, and the haploid case.
- ``-i [NAME=]PATH`` form for inputs: explicit names take precedence,
  otherwise hap labels are inferred from filename stems with a
  positional ``input1`` / ``input2`` fallback. The reserved label
  ``unassigned`` is only ever produced by explicit
  ``-i unassigned=PATH`` — auto-inference never produces it.
- ``karyoscope.core.hap_inference``: pattern library + per-contig
  classifier. Case-insensitive built-in regexes match hifiasm
  (``h[12]tg``), explicit ``hap1`` / ``hap2`` tags, and verkko-style
  ``maternal`` / ``paternal`` (and short forms ``mat`` / ``pat``)
  against contig names and filename stems. ``maternal`` and
  ``paternal`` are kept as distinct labels rather than mapped to
  ``hap1`` / ``hap2`` (the karyotype renderer can map them at render
  time). For a single combined-file input the rules split into
  multiple haps if patterns split the contigs cleanly; otherwise all
  contigs become ``hap1`` with a warning. ``--split-haps REGEX``
  overrides the built-in patterns with a user-supplied capture group.
- ``karyoscope.core.scaffold``: a faithful port of the archive's
  ``scaffold_stats.py`` (the newer version with the chrY-centroid
  fix). Public API: pure helpers (``get_simple_region``,
  ``assign_main_chromosome``, ``find_largest_contiguous_region``,
  ``half_region_totals``, ``need_to_flip``, ``flip_bins``,
  ``category_index``), the high-level entry point
  ``classify_and_orient`` (takes pre-binned BEDs + telomere flags,
  returns a list of :class:`MapRow` in canonical order), and the BED
  rewriter ``rewrite_bed`` (uses the explicit map file rather than
  parsing the encoded name).
- ``karyoscope.core.scaffold_run``: the high-level orchestrator that
  the CLI sits on top of. Implements the auto-derive cascade:
  per-input, missing annotation BEDs trigger ``annotate``, missing
  telomere files trigger ``seqtk telo``, missing 1 Mb binned BEDs
  trigger ``bin_features`` in-process. ``--no-auto`` turns missing
  prerequisites into hard errors. The cascade is structured as a
  list of independent per-input units so the future move to batched
  parallel execution is a mechanical loop change rather than a
  rewrite.
- ``karyoscope.core.io.scaffold_map``: ``MapRow`` dataclass plus
  ``write_map`` / ``read_map`` for the 8-column TSV format
  (``new_name``, ``original_name``, ``input_file``, ``hap``,
  ``chromosome``, ``flipped``, ``length``, ``stats``) and
  ``write_legacy_stats`` for the 2-column archive format.
- ``karyoscope.core.io.telo``: parser for the seqtk-telo 3-column
  output and ``run_seqtk_telo`` shellout via the existing
  ``require_tool`` / ``run_tool`` machinery.
- Manifest-driven role resolution: scaffold reads
  ``roles.chromosome_assignment`` and ``roles.region_assignment``
  from ``manifest.yaml`` to pick which feature sets to use for
  chromosome classification and orientation respectively. Falls back
  to literal feature-set names ``"chromosome"`` and ``"region"`` with
  a warning when the manifest omits them. Errors when the resolved
  set is not declared in the manifest's ``feature_sets``.
- ``ScaffoldError(KaryoscopeError)`` for problems classifying,
  orienting, or rewriting per-input outputs.
- ``--acrocentric CHROM`` flag (repeatable, accepts comma-separated
  lists) overrides the default human acrocentric set
  (``chr13``, ``chr14``, ``chr15``, ``chr21``, ``chr22``). The
  default-set fallback emits a single warning suggesting that
  non-human assemblies set it explicitly.
- ``--telo [NAME=]PATH`` flag (repeatable) lets users supply a
  precomputed ``seqtk telo`` output to skip the auto-run. Telo
  entries match inputs by name; unmatched ``--telo`` names error
  cleanly.
- Deliberate diffs from the archive's ``scaffold_stats.py`` and
  ``scaffold_features_cli.py``: leaf-chromosome detection is
  structural (``leaves_for(hierarchy, chromosome_feature_set)`` rather
  than ``startswith("chr")``), so non-human / non-``chr``-prefixed
  databases work; the BED rewriter does not round-trip through the
  encoded contig name (which the archive did via string-parsing,
  fragile when original contig names contained ``_``); the
  ``--priority`` file path is gone in favour of ``--acrocentric``.
- Unit tests for the scaffold algorithm (``tests/test_scaffold.py``):
  every pure helper (``get_simple_region``, ``chromosome_sort_key``,
  ``assign_main_chromosome``, ``find_largest_contiguous_region``,
  ``half_region_totals``, ``scaffold_region_majority``,
  ``need_to_flip`` across all branches of the boolean ladder,
  ``flip_bins``, ``category_index``), ``classify_and_orient``
  end-to-end (basic keep/drop, RC suffix on flips, ordering across
  chromosomes and haps), and the ``rewrite_bed`` reverse-mirror
  semantics. 44 tests.
- Unit tests for the map format (``tests/test_scaffold_map.py``, 9
  tests), seqtk telo parser (``tests/test_telo.py``, 6 tests), and
  hap inference (``tests/test_hap_inference.py``, 29 tests including
  built-in pattern coverage, single-input split inference, and
  multi-input filename-stem heuristic with collision fall-through).
- Integration tests for the CLI
  (``tests/test_command_scaffold.py``): argument-parsing unit tests
  (named-pair form, telo without name, unmatched telo) plus three
  end-to-end tests marked ``@pytest.mark.integration`` that exercise
  the full auto-derive cascade (annotate -> seqtk telo -> bin ->
  scaffold) against the dummy database, plus the ``--no-auto`` error
  path and the map-file format round-trip.
- ``karyoscope scaffold --mode {fasta, bed, both}``. Stage 5d-1
  shipped only the BED-rewriting path; this stage adds FASTA
  rewriting and flips the default mode to ``fasta``. Rationale:
  when a human invokes ``karyoscope scaffold`` directly they want
  a scaffolded FASTA; the pipeline use (``karyoscope karyotype``
  driving scaffold under the hood) will explicitly pass
  ``--mode bed``. Combining ``--mode fasta`` with ``--feature-set``
  is a clean ``UsageError`` (no per-feature-set BEDs are written
  in fasta mode, so feature-set filtering would have no effect).
- FASTA-mode output:
  ``<input_stem>.<dbid>.scaffolded.fa[.gz]`` per input. The map and
  legacy stats files are still written in every mode. ``--mode both``
  produces both the FASTA and the per-feature-set BEDs.
- ``--keep-unscaffolded / --drop-unscaffolded`` flag (default
  ``--keep-unscaffolded``): in FASTA mode, contigs that the
  classify-and-orient step dropped (too short, or no leaf-chromosome
  hit) are by default appended at the end of the output FASTA under
  their original names. This matches the archive's
  ``scaffold_hap_assembly.py`` behaviour and keeps the assembly
  complete; pass ``--drop-unscaffolded`` to emit only the
  scaffolded contigs.
- ``karyoscope.core.io.fasta``: minimal FASTA reader / writer +
  IUPAC-aware reverse-complement. ``read_fasta_records`` returns
  an :class:`collections.OrderedDict` so source order can be
  preserved on round-trips. ``write_fasta_records`` supports plain
  and gzipped output with optional line wrapping (default:
  unwrapped, matching the archive). ``reverse_complement`` handles
  the full IUPAC table (ACGTU + RYSWKMBDHVN) and preserves case so
  soft-masked sequences round-trip cleanly. The lighter
  ``read_fasta_contig_names`` lives here too;
  :mod:`karyoscope.core.hap_inference` re-exports it for back-compat.
- ``karyoscope.core.scaffold.rewrite_fasta``: walks the per-input
  map in stats order, reverse-complementing flipped contigs, and
  optionally appending unscaffolded contigs at the end under their
  original names.
- FASTA-mode optimisation in ``scaffold_run``: when
  ``mode='fasta'``, the auto-derive cascade only requests the role
  feature sets from ``annotate`` (rather than every set declared
  in the manifest), since per-feature-set BEDs aren't going to be
  rewritten anyway. ``mode='bed'`` and ``mode='both'`` request all
  user-asked-for feature sets as before.
- Unit tests for the new fasta module
  (``tests/test_fasta.py``, 21 tests covering reverse-complement
  including IUPAC ambiguity codes and case preservation, FASTA
  read/write with plain + gzip + CRLF handling, line wrapping,
  insertion-order preservation), plus 7 ``rewrite_fasta`` tests
  added to ``tests/test_scaffold.py`` (no-flip rename, RC on
  flip, keep-unscaffolded vs drop semantics, map-order emission,
  map contigs absent from FASTA, gzip in/out).
- Integration tests in ``tests/test_command_scaffold.py``:
  ``--mode fasta`` writes the FASTA and no BEDs; ``--mode both``
  writes both; ``--feature-set`` with ``--mode fasta`` is
  rejected with a clean ``UsageError``;
  ``--drop-unscaffolded`` omits leftover contigs from the FASTA.

### Changed
- Test fixture ``tests/data/dummy_db.tar.gz`` rebuilt from scratch.
  ``features.tsv`` now uses the correct one-row-per-id schema
  (``featureID`` plus one column per feature set), and the KMC index
  was regenerated using a deliberately-designed three-sequence FASTA
  whose 21-mer occurrence counts are exactly 1, 2, and 3 — so the
  KMC counters themselves become the featureIDs in ``features.tsv``.
  See the docstrings of ``DUMMY_SEED`` and ``_run_kmc`` in
  ``tests/data/build_dummy_db.py`` for the construction. The
  ``hierarchy.tsv`` was subsequently rewritten in the new
  ``(feature_set, child, parent)`` schema with a three-level region
  hierarchy (``rA → aSat → centromeric → categorized``,
  ``rC → HSat → centromeric → categorized``) so smoothing tests can
  exercise non-trivial LCA promotion. SHA-256 updated in
  ``tests/data/dummy_db.sha256`` and ``tests/data/dummy_registry.yaml``.
- Terminology standardisation: smoothing code uses ``sequence``,
  ``seq_name``, ``seq_id`` throughout rather than ``read`` /
  ``contig`` / ``scaffold``. Domain terms are kept where they're
  domain-correct (e.g., the ``scaffold`` subcommand operates on
  contigs).
- Lint and format tooling is now managed by [`pre-commit`](https://pre-commit.com).
  Tool versions are pinned in `.pre-commit-config.yaml` (currently ruff 0.15.13
  and pre-commit-hooks v6.0.0); CI runs `pre-commit run --all-files`, so
  developers and CI use identical tool versions. See `CONTRIBUTING.md` for
  setup instructions.
- `ruff` removed from `pyproject.toml`'s `[dev]` extras since pre-commit
  manages its own isolated installations.
- `tests/conftest.py` now surfaces a clear actionable error when pytest is
  launched from a Python interpreter where `karyoscope` is not installed
  (a common pitfall on macOS where multiple Python installations coexist).
  `CONTRIBUTING.md` updated to recommend `python -m pytest` over bare
  `pytest` to sidestep the issue.
- Default `pytest` invocations now deselect tests marked
  `@pytest.mark.integration`. CI runs them explicitly via
  `pytest -m integration` after building the C++ helper. This keeps the
  unit-test workflow fast and friendly for Python-only contributors,
  while still ensuring the full pipeline is exercised in CI.

### Fixed
- `paths.default_db_root()` now resolves `~` at call time rather than module
  import time, so `$HOME` changes during a process lifetime are honored.
- `CITATION.cff` and `CODE_OF_CONDUCT.md` normalized to a single trailing
  newline (auto-fixed by the new `end-of-file-fixer` pre-commit hook).
- Test fixtures (`unpacked_dummy_db`, `populated_db_root`) no longer crash
  on Python versions older than the `tarfile.TarFile.extractall(filter=...)`
  backport (3.10.12 / 3.11.4 / 3.12+). A shared `_extractall_compat`
  helper in `conftest.py` wraps the call in a `try`/`except TypeError`
  fall-through, mirroring the same pattern that production code in
  `download.py` already uses.
- `native/get_featureIDs/Makefile` now adds `-Wl,-rpath,$(CONDA_PREFIX)/lib`
  to the link flags when `CONDA_PREFIX` is set. Without this, on
  macOS + conda the binary built fine but failed at runtime with
  "dyld: Library not loaded: @rpath/libz.1.dylib ... no LC_RPATH's
  found" — conda-forge dylibs use `@rpath/libfoo.dylib` install-names
  and require the consuming binary to embed an `LC_RPATH`. Conda's own
  compilers do this automatically; the system `c++` doesn't.
  `native/README.md` documents the gotcha for future reference.

<!--
Use the following sections when adding entries:

### Added
### Changed
### Deprecated
### Removed
### Fixed
### Security

For releases, copy the [Unreleased] section to a new heading like:

## [1.0.0] - 2026-MM-DD
-->

[Unreleased]: https://github.com/barthel-lab/KaryoScope/compare/HEAD...HEAD
