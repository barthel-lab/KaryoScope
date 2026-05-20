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
- ``karyoscope centromeres`` command, replacing the previous stub.
  Per scaffolded contig, identifies the centromere's start/end
  coordinates from the binned scaffolded region BED produced by
  ``scaffold``. Coarse-pass (1 Mb bins by default) finds
  ``min(start)`` / ``max(stop)`` of bins classified as
  ``"centromere"`` via :func:`karyoscope.core.scaffold.get_simple_region`;
  an optional fine-refinement pass (100 kb bins by default) tightens
  the call within a ``+/- 1 Mb`` window around the coarse range.
- Output: ``<input_stem>.<dbid>.centromeres.bed[.gz]`` per input.
  3-column BED (contig, start, end) -- coordinate-only, no feature
  column, so the filename has no ``<feature_set>`` segment (matches
  the project convention that a feature-set segment in a filename
  means the file *contains* records from that set). Contigs with no
  centromeric content are omitted; centromere coordinates are in
  the scaffolded (post-flip) coordinate system so
  ``karyoscope karyotype --mode centromere`` can consume them
  directly.
- Like ``karyoscope scaffold``, the command auto-derives its
  prerequisites: missing scaffolded BED triggers ``scaffold_run``
  with ``mode="bed"`` and only the centromere feature set, which
  itself cascades through annotate, seqtk telo, and bin; missing
  binned scaffolded BEDs trigger ``bin_features`` in-process.
  ``--no-auto`` turns missing inputs into hard errors.
- ``karyoscope.core.centromeres`` exposes
  :func:`find_centromere_ranges` (pure algorithm: takes per-contig
  binned intervals, returns ``OrderedDict[contig, (start, end)]``
  with insertion order preserved from the input) and
  :func:`centromeres_run` (the orchestrator).
- ``CentromereError(KaryoscopeError)`` for problems extracting
  centromere coordinates.
- New CLI flag ``--centromere-feature-set`` overrides the manifest
  role (rather than reusing the ``--feature-set`` name from other
  commands, which conventionally means "what gets written"; here it
  means "what drives detection").
- Manifest role chain: ``roles.centromere_detection`` ->
  ``roles.region_assignment`` -> the literal name ``"region"``. Each
  fallback step emits a warning so the user can add the role to
  ``manifest.yaml`` to silence it.
- ``--coarse-bin-size`` (default 1 Mb) and ``--fine-bin-size``
  (default 100 kb; pass 0 to disable the fine pass) flags.
- Unit tests in ``tests/test_centromeres.py`` (13 tests):
  ``find_centromere_ranges`` coarse-only (single block, contig
  without centromere omitted, insertion-order preservation,
  non-contiguous outer-bounds semantics), coarse + fine
  refinement (fine narrows coarse, out-of-window fine ignored,
  no-fine-signal falls back to coarse, missing contig in fine,
  ``refinement_buffer=0`` strict-inside behaviour), plus CLI
  parse tests. Two integration tests marked
  ``@pytest.mark.integration`` exercise the full
  FASTA -> scaffold -> bin -> centromeres cascade against the
  dummy database and the ``--no-auto`` error path.
- ``karyoscope karyotype`` command, replacing the previous stub.
  Top of the cascade: given one or more FASTA inputs, runs the
  full pipeline through to a karyotype SVG (annotate -> seqtk
  telo -> bin -> scaffold mode='bed' -> centromeres for the
  centromere mode -> render). Existing intermediate files are
  reused; ``--no-auto`` turns missing inputs into hard errors.
- Three render modes:

  * ``full`` (default): whole-chromosome view, 1 Mb bin size,
    10 Mb scale bar.
  * ``subtelomere``: zoomed p/q-arm telomere view, 100 bp bin
    size, 10 kb scale bar. Only contigs flagged with at least
    one telomere appear.
  * ``centromere``: zoomed per-contig centromere view, 100 kb
    bin size, 1 Mb scale bar. Only contigs with centromere
    coordinates appear.
- One SVG written per requested feature set
  (``--feature-set``, repeatable; default: every feature set in
  the manifest). Output naming:
  ``<base>.<dbid>.<mode>.<feature_set>.karyotype.svg``, where
  ``<base>`` is ``"karyotype"`` by default or the basename of
  ``--output PATH`` when given.
- Sex-determination support (port of the archive's logic):
  ``--sex {male, female, reference, unknown}`` and
  ``--sex-determination-system {XY, X0, ZW, ZO}``. Per chromosome,
  decides which (chrom, hap) cells get rendered, drawn empty, or
  skipped. ``unknown`` draws sex-chromosome cells only where data
  is present (no expectations imposed).
- ``--background-color {white, black}``: white draws sequence
  outlines, black omits them and uses light text + scale bar for
  contrast.
- ``--subtelomere-boundary`` flag for the subtelomere window
  size; ``--bin-size`` to override the mode-appropriate default.
- ``--no-human-chroms`` to skip seeding the chromosome list with
  the standard human set (chr1..chr22, chrX, chrY) -- for
  non-human assemblies where the chromosomes come entirely from
  the data.
- ``karyoscope.core.io.colors``: parser for ``colors.txt`` (3-column:
  feature_set / feature / color). Returns a nested
  ``{feature_set: {feature: hex}}`` mapping; ``colors_for_set``
  flattens it for one feature set and pre-populates the
  ``novel -> #ffffff`` sentinel. Missing colours render as white
  with a one-time warning per feature so a single bad colour
  entry doesn't break the whole render.
- ``karyoscope.core.karyotype``: the renderer. Faithful port of
  ``KaryoScope_assembly.py`` adapted for v0.1's per-input + map-file
  topology. ``render_karyotype(inputs, *, colors, mode, sex, ...)``
  takes :class:`RenderInput` records (each carrying per-input
  ``map_rows``, binned BED, and optional centromere coordinates)
  and writes the SVG using ``drawsvg``. Pure-ish: no subprocess,
  no in-place mutation of arguments.
- ``karyoscope.core.karyotype_run``: the orchestrator that the
  CLI sits on top of. Cascades through ``scaffold_run`` (for
  scaffolded BEDs), ``bin_features`` (for the per-mode binned
  scaffolded BEDs), and ``centromeres_run`` (for centromere
  mode). One SVG written per requested feature set.
- ``KaryotypeError(KaryoscopeError)`` for render-time problems
  (unknown mode, missing centromere ranges, etc.).
- Unit tests in ``tests/test_karyotype.py`` (19 tests):
  ``get_expected_haps`` across every (sex, sex-system, chromosome)
  branch including custom dict systems and the unknown-system
  error path, ``render_karyotype`` smoke tests per mode
  (full / subtelomere skips non-telomere contigs / centromere
  requires ranges), unknown-feature-renders-white, CLI parse
  errors (missing input, --outdir/--output conflict). Two
  integration tests marked ``@pytest.mark.integration`` exercise
  the full FASTA -> render cascade against the dummy database
  for ``--mode full`` and ``--mode centromere``.
- Unit tests for the colors parser
  (``tests/test_colors.py``, 10 tests): basic parsing,
  same-feature-different-sets, missing/empty/bad-header files,
  wrong column count, blank values rejected, novel-sentinel
  pre-population, user-override of novel, unknown feature set
  yields novel-only.
- ``test_cli.py`` no longer has the ``test_stub_subcommands_exit_cleanly``
  test -- every subcommand is now a real implementation. The
  Stage 5d series is complete.

### Changed
- ``karyoscope`` now treats missing colour entries as malformed-
  database errors rather than silently rendering white. Surfaced
  during dogfooding when ``repeat`` (a top-level grouping in the
  production hierarchy) was missing from colors.tsv and the
  legend showed a white swatch -- visually indistinguishable from
  the ``novel`` sentinel. The fallback was wrong by design.
- New :func:`karyoscope.core.io.colors.validate_colors` cross-checks
  hierarchy.tsv against colors.tsv: every hierarchy node (children
  + parents, so including the ``categorized`` root) must have an
  entry in colors.tsv for its feature set. The ``novel`` sentinel
  is exempt -- it's always rendered white.
- ``karyoscope info <db>`` now runs ``validate_colors`` and prints
  any issues under a ``Colors warnings:`` heading (informational,
  doesn't fail the command).
- ``karyoscope download`` runs ``validate_colors`` after extracting
  the tarball and **refuses to register** the install (no entry
  added to installed.json) if any colours are missing. The
  extracted directory is left on disk for inspection; the error
  message lists every missing ``(feature_set, node)`` pair.
- ``karyoscope karyotype`` runs ``validate_colors`` before any
  render starts, raising :class:`KaryotypeError` if any colour is
  missing. Catches malformed databases that bypassed
  ``download`` validation (e.g. our manual installed.json
  bootstrap on the cluster).
- :func:`render_karyotype._color_for` now raises
  :class:`KaryotypeError` on missing colours instead of falling
  back to white. ``karyotype_run`` validates upstream so the
  raise is defence-in-depth for direct callers.
- Test fixture ``tests/data/dummy_db.tar.gz`` rebuilt to add the
  missing internal nodes (``categorized``, ``autosome``,
  ``centromeric``, ``aSat``, ``HSat``) to colors.txt so the dummy
  db itself passes the new validator. New SHA-256 propagated to
  ``dummy_db.sha256`` and ``dummy_registry.yaml``.
- ``karyoscope karyotype`` title-band wording: each metadata segment
  now carries an explanatory noun for consistency with ``genome view``.
  Format: ``<sample>  |  <dbid> database  |  <mode> view  |
  <feature_set> feature set  |  smoothed``. The plain ``<dbid>`` /
  ``<feature_set>`` form from the first dogfooding pass was
  ambiguous to scan.
- ``karyoscope karyotype`` legend layout: gap between karyotype and
  legend is now ``chrom_gap`` (same as between adjacent chromosomes,
  so the legend reads as just another column) rather than the
  previous double-padded ``x_border + 18`` (~82 px of wasted space).
  Legend column width is computed dynamically from the longest
  drawn label rather than a fixed 170 px, removing the wasted right
  margin in the SVG.
- ``karyoscope karyotype`` legend sort: features are now ordered by
  the ``hierarchy.tsv`` ``child`` column in file order (the
  database author's intended ordering -- internal nodes first, then
  the subtree under each). This groups related features together
  visually. The previous natural-chr-then-alpha sort is still the
  fallback when no order is given. For chromosome feature sets the
  hierarchy.tsv naturally lists ``chr1``, ``chr2``, ...,
  ``chrX``, ``chrY`` in order so the user-visible behaviour
  matches the previous natural sort. No special case needed.
- ``karyoscope karyotype`` SVGs now carry a title band at the top
  (centred over the karyotype area) listing sample, database, mode,
  feature set, and a ``smoothed`` flag, plus a colour legend in the
  right margin listing only features that actually appear in the
  rendered data (so the legend stays compact for large feature
  sets). Three new flags: ``--sample-label TEXT`` overrides the
  auto-derived label (default: joined input-FASTA stems);
  ``--no-title`` and ``--no-legend`` suppress the new bands when
  not wanted. Layout impact: title adds ~35 px to the top, legend
  adds ~170 px to the right; visible differences in image
  dimensions vs prior renders.
- ``karyoscope karyotype --format {svg,pdf,png}`` (repeatable;
  default ``svg``). PDF and PNG are produced by converting the SVG
  via ``cairosvg`` (already a Python dep), which requires the
  native ``libcairo`` library at runtime -- install with
  ``conda install -c conda-forge cairo`` in the active env. When
  requested formats include only non-SVG, the intermediate SVG is
  written to its conventional path, converted, then deleted; pass
  ``--format svg --format pdf`` to keep both.
- ``KaryotypeResult.output_paths`` replaces ``svg_path`` (now a
  property pointing at ``output_paths[0]`` for back-compat). One
  ``KaryotypeResult`` per (mode, feature_set) carries the full list
  of files written for that combination.
- ``karyoscope bin`` gained per-sequence-chunk parallelism via
  :class:`multiprocessing.Pool`. New CLI flag ``-t / --threads``
  (default 1; 0 = auto via ``os.cpu_count()``). Threaded output is
  byte-for-byte identical to the single-threaded path (verified by
  ``tests/test_bin.py::TestBinFeaturesThreaded``). Stdin/stdout I/O
  remains single-threaded regardless of the flag (the pool path
  requires a real on-disk input). The worker pool reuses the same
  ``chunked_seq_reader`` from :mod:`karyoscope.core.smooth`: chunks
  always end at a sequence boundary so per-sequence binning
  semantics are preserved without needing a bin-boundary-aware
  reader. ``scaffold_run``, ``centromeres_run``, and
  ``karyotype_run`` propagate their ``threads`` parameter to
  ``bin_features`` so the per-sequence parallelism flows through
  the auto-derive cascade end-to-end.
- ``karyoscope karyotype`` renamed its ``full`` render mode to
  ``genome``. The CLI now accepts ``--mode genome`` (case-insensitive)
  and the output filename segment is ``.genome.`` rather than
  ``.full.``. The internal API (:func:`render_karyotype`,
  :data:`ALL_MODES`, :data:`DEFAULT_BIN_SIZE_BY_MODE` keys, the
  ``Mode`` literal type) was updated to match. Hard rename with no
  back-compat: callers that used ``mode="full"`` need to switch to
  ``mode="genome"``.
- ``karyoscope karyotype --mode`` is now repeatable, and the default
  (no ``--mode`` flag) renders **every** mode rather than just one.
  Combined with ``--feature-set`` (already repeatable), the default
  ``karyoscope karyotype -i hap1.fa -i hap2.fa --sex male`` invocation
  now produces ``len(modes) x len(feature_sets)`` SVGs. Restrict
  either axis to subset.
- ``karyoscope karyotype --bin-size`` is now only valid when exactly
  one ``--mode`` is specified, since different modes have different
  natural bin sizes (1 Mb / 100 kb / 100 bp for genome / centromere /
  subtelomere respectively). Combining ``--bin-size`` with multiple
  modes raises ``UsageError`` rather than silently applying the
  override to all of them.
- The auto-derive cascade in ``karyotype_run`` only calls
  ``centromeres_run`` when ``"centromere"`` is in the requested
  modes (it's expensive: an extra bin pass per input). Previous
  behaviour ran it unconditionally because the caller could only
  specify a single mode at a time.
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

### Changed
- ``karyoscope scaffold -v``, ``karyoscope centromeres -v``, and
  ``karyoscope karyotype -v`` now bookend every major phase with
  start + end INFO lines, matching the depth of the annotate logging.
  Per-pipeline additions:
    * **scaffold**: opening banner (``scaffolding N input(s) ...``);
      ``classifying + orienting ... contigs`` / ``classified ... in Xs``;
      per-(input, FS) ``rewriting scaffolded BED`` / ``wrote ... in Xs``;
      per-input ``wrote <scaffolded.fa> in Xs`` (paired with the existing
      ``writing scaffolded FASTA`` start log); closing
      ``scaffold complete in Xs``.
    * **centromeres**: opening banner; per-input ``finding centromere
      ranges`` / ``found N centromere range(s) in Xs``; closing
      ``centromeres complete in Xs``.
    * **karyotype**: opening banner with the total ``N renders x M
      formats`` count; per-(mode, FS) ``rendering karyotype`` /
      ``rendered ... in Xs``; per-format ``converting`` / ``converted
      ... in Xs`` (PDF / PNG via cairosvg can take several seconds
      each); closing ``karyotype complete in Xs``.
  Also adds completion lines to two reusable building blocks:
  ``bin_features`` now logs ``binned ... in Xs`` after the existing
  start line, and ``run_seqtk_telo`` logs ``ran seqtk telo on ... in
  Xs``. The redundant ``binning X -> Y (bin_size=Z)`` line previously
  emitted by every ``_ensure_binned*`` helper (in scaffold_run,
  centromeres, karyotype_run) is dropped -- ``bin_features``'
  own start line carries strictly more information (leaf_set + thread
  count), so the helper-level announcement was just noise.
- ``karyoscope annotate -v`` now logs at INFO level for every major
  pipeline phase: ``running get_featureIDs`` (start) /
  ``ran get_featureIDs in Xs (combined BED: N GB)`` (end),
  ``smoothing pass: N feature set(s)`` /
  ``smoothing pass complete in Xs``, ``bgzip pass: N BED(s)`` /
  per-file ``bgzipping ...`` / ``bgzipped ... (input -> output) in Xs``
  / ``bgzip pass complete in Xs``, and a closing
  ``annotate complete in Xs (N output BED(s))``. Previously the C++
  k-mer query and bgzip steps were completely silent at INFO; long
  whole-genome runs looked indistinguishable from a hung process.
- ``karyoscope annotate`` smoothing is now driven by **one worker pool
  for all feature sets**, replacing the previous per-feature-set loop.
  Each worker is initialised once with the full
  ``{feature_set: HierarchyIndex}`` /
  ``{feature_set: FeaturesForWorker}`` state and processes every
  requested set per chunk in a single invocation. Wins:
    * One pool spawn instead of ``N`` (saves the spawn-context init
      cost ``N-1`` times -- on whole-genome inputs this is ~30 s ×
      ``N-1`` of fresh Python process startup plus initargs unpickle).
    * One pass over the combined BED instead of ``N`` (saves
      ``N-1`` × multi-GB sequential I/O).
    * No inter-feature-set idle gaps (previously ~2 min per
      transition on HG002 while the old pool drained, per-sequence
      temp files were concatenated, and the next pool initialised).
  In assembly mode workers write each ``(feature_set, sequence)``
  pair's output directly to a per-feature-set temp directory using a
  sanitised sequence-name filename; the main process concatenates the
  per-sequence files in input-FASTA order at the end. Reads mode
  (``--no-preserve-order``) keeps the IPC return path. Output BEDs
  are byte-equivalent to the previous per-feature-set code.

### Fixed
- ``karyoscope annotate`` no longer hangs silently for hours when a
  smoothing worker is killed by an external signal (most commonly
  SIGKILL from the kernel OOM-killer on memory-constrained nodes).
  ``multiprocessing.Pool`` does not reassign tasks whose worker died,
  so ``pool.imap_unordered`` would block forever waiting on results
  that will never arrive -- observed on whole-genome HG002 runs
  against a 50 GB SLURM allocation, where ``-t 16`` blew through the
  memory ceiling on chr1/chr2 ``region`` tasks. A daemon watchdog
  thread now polls the pool every 2 s; on detecting a dead or
  pool-replaced worker it writes an actionable stderr message
  (likely cause + remediation: reduce ``--threads`` or increase RAM)
  and ``os._exit``s the process with code 137. The hang turns into
  a loud, immediate failure within seconds. Uses ``pool._pool`` (a
  documented-private attribute stable since Python 3.4) for worker
  inspection.
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
