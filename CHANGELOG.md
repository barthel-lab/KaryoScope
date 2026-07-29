# Changelog

All notable changes to KaryoScope will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **`annotate` reads CRAM.** `--input` accepts `.cram`, and a new `--reference`
  option supplies the FASTA the alignment was encoded against. CRAM stores
  bases as a diff against its reference, so `--reference` is **required** for
  CRAM input and the CLI refuses the run up front without it. That check is
  deliberately loud, because the quiet failure is worse than the noisy one:
  left to itself, htslib resolves each contig's M5 checksum through
  `$REF_PATH`/`$REF_CACHE` and can decode against a *different* build of the
  same genome, producing plausible sequence that is simply wrong.

  Note `--reference` is the **alignment** reference and is unrelated to `--db`.
  `annotate` is alignment-free, so annotating GRCh38-aligned reads against a
  CHM13 database is a normal thing to do, not a mismatch.

  Both backends are covered. KMC streams `samtools fasta` straight into
  `get_featureIDs` with no temp file; HKS cannot (`hks lookup` needs a seekable
  path) and so materialises a temp FASTA in `$TMPDIR` first. That file is
  **full size** — a 64x human WGS CRAM measured 28.9 GB in, 254 GB out — so
  `$TMPDIR` must point at node-local scratch, never at a shared filesystem.

- **`--query-names/--no-query-names`** controls which identifier the output BED
  carries, and **`--query-names` is now refused outright for read-level input**
  (FASTQ/BAM/CRAM). It is not a tuning knob there — it is a way to lose a run.

  `hks --report-query-names` is not streaming: `load_seq_names` reads the whole
  query file and builds a `Vec<String>` of every sequence name before the first
  k-mer is looked up. Fine for an assembly's few thousand contigs; fatal for
  reads. Measured on a 64x human WGS CRAM — 1.315 G names at ~68 bytes each is
  ~90 GB on top of the ~12.3 GB index. The run was OOM-killed ~20 minutes in,
  having already spent 15 of them decoding the CRAM, with nothing written and
  nothing to resume from.

  Raising the memory limit is not a fix: it would take a ~160 GB allocation per
  sample to buy nothing but longer identifiers. So this is an error rather than
  a warning, and it fires in a second instead of after the decode.

  No information is lost. Rank N is the Nth record of the query file, and for
  an alignment input that file is a deterministic function of the source, so
  the mapping regenerates on demand:
  `samtools fasta -F 0x900 -N --reference <ref> <input> | grep '^>'`.

  Assemblies are unaffected: they still default to names, and
  `--no-query-names` switches them to ranks to save space.

- **`annotate` now checks available memory before it starts**, not just disk.
  For the HKS backend the requirement is knowable exactly rather than
  estimated: `hks` holds the shared base index plus one feature set's labeling
  at a time, so the peak is the sum of two file sizes read off disk (measured
  9.73 GB predicted vs 10.04 GB observed on the human database). It resolves
  the limit from `$SLURM_MEM_PER_NODE`/`_PER_CPU`, then the cgroup limit, then
  `MemAvailable` — and fails open if either the index or the limit cannot be
  determined, so it can only stop a run that would have been killed anyway.
- **A process killed by a signal now says so.** `subprocess` reports `-N` for
  signal N, so an out-of-memory kill surfaced as `command failed with exit code
  -9` — not an exit code at all, and unsearchable. Failures now read `killed by
  SIGKILL (signal 9)`, including the shell-style `137` that SLURM and Docker
  report, with backend-specific advice on how much memory to request. This
  previously existed for the KMC backend only; it now applies to every external
  tool KaryoScope runs.

- **`karyoscope prep-bed pave` — a papillomavirus `gene` feature set from PaVE
  genome records.** Leaves are the ORFs — `E6`, `E7`, `E1`, `E2`, the
  genus-specific `E5_*`, `E10`, `L2`, `L1` — plus the `URR`, grouped `early` and
  `late`.

  The spliced transcripts `E1^E4`, `E8^E2` and `E6*` are excluded: they lie
  wholly inside the ORFs they are spliced from — 100% of their bases across
  PaVE's 224 human reference genomes — so a ranking that prefers the primary ORF
  leaves them holding no sequence. The `E1BS` and `E2BS` binding motifs are
  excluded as 12–20 bp, shorter than any usable k. A feature name with no leaf
  raises rather than being dropped to the gap-fill.

  The reading frames overlap, so `--priority` writes a ranking in genome order
  rather than the set being flattened: 10.2% of bases carry more than one
  annotated feature, 1.2% once the spliced transcripts are dropped. A feature
  spanning the origin of a circular genome becomes two records with one label.

  A PaVE record also carries the genome sequence and the ICTV lineage, so
  `--fasta` writes the FASTA for `build`'s `sequence:` and `--taxonomy` writes
  the genus → species → type hierarchy for a `type` set whose BED comes from
  `prep-bed fai`.

- **A recipe for an HPV database, `docs/recipes/hpv-pave.md`.** 224 genomes,
  1.68 Mb, a 2 s build and a 10.5 MB database with `type` and `gene` sets. Both
  the input records and every derived file are pinned by checksum, and the
  download script ships in `docs/recipes/files/`.

- **`karyoscope prep-bed asat` — a per-array alpha-satellite feature set from a
  CenSat annotation.** Where `prep-bed censat` collapses CenSat onto its 14
  broad classes, `asat` reads the same file per array: `hor_1_5(S1C1/5/19H1L)`
  becomes a record labelled `S1C1_5_19H1L`. On CHM13 v2: 766 records, 86 arrays,
  85,444,299 bp.

  Both CenSat dialects are read. An interval naming several arrays produces one
  record per array over the full interval, so HKS resolves the shared k-mers to
  their common ancestor. Bare continuation suffixes are expanded first, so
  `hor_1_1(S3C1H2-A,B,C)` names S3C1H2-A, S3C1H2-B and S3C1H2-C rather than
  leaves called `A`, `B` and `C`.

  Every alpha-satellite class is included by default, nested under a shared
  `asat` parent. Excluding one leaves that sequence to the gap-fill, whose leaf
  sits at the hierarchy root, so k-mers the named arrays share with it resolve
  to the root. The shared sequence is the conserved monomer core: on CHM13,
  excluding `mon` and `dhor` puts 36.9% of array bases on the root, against 4.6%
  with them included.

- **`flatten_order:` (build spec) and `--flatten-order` (CLI) set a feature
  set's per-base flatten ranking independently of `priority:`.** The two settings
  answer different questions: `flatten` picks one label per base before k-mers
  are extracted, while `priority` resolves a k-mer claimed by several labels as
  the index is built. A set may want a different ranking for each. Supplying
  `flatten_order:` implies `flatten` for that set.

  Previously `flatten` derived its order from the priority file, so a set that
  wanted a specific per-base partition had to accept the same ranking at index
  time. The `repeat` set in the CHM13 recipe is the case in point: its
  fully distinct 1..15 leaf ranking is what makes the set a hard partition, but
  applied per k-mer it collapses every shared k-mer onto one class rather than
  letting siblings tie and resolve to their common ancestor.

- **`karyoscope info <archive>` validates a database archive's layout.** Given a
  `.tar.gz`, `.tgz`, or `.tar`, `info` now reports the top-level directory, the
  extracted size, whether the layout is valid, and the same manifest and
  feature-set summary it prints for a database directory. It reports
  `Layout valid: NO` for an archive with more or fewer than one top-level
  directory, without a `manifest.yaml`, or missing any file the manifest
  declares, and warns when the manifest `id` and the directory name disagree.

  Previously an archive fell through to the generic file branch and printed only
  `Type: file (435.7 MB)`. The registry's contributor instructions told people to
  run this command and fix what it reported, so the check passed for any tarball.

  The archive is streamed once: the manifest and TSV/TXT sidecars are read into a
  temporary directory and the index files are checked by name, so a multi-gigabyte
  database is verified without being unpacked. Members whose paths escape the
  archive root are skipped rather than written.

  Alongside the layout it reports an `Install check`, which runs the same
  hierarchy-versus-colours validation `download` applies before recording an
  install. A database can be structurally valid and still be refused at install;
  that now shows up before publishing rather than after. Both paths call the new
  `manifest.check_install_readiness()`, so they cannot drift apart.

### Changed

- **`build` now rejects a priority file that omits any hierarchy node**, listing
  the nodes it did not find. The root is exempt, since it is never a child and
  its priority is never compared. Existing specs whose priority files are
  incomplete will fail until the missing nodes are added.

  An omitted node previously took priority 0. HKS treats a lower value as higher
  priority, so a node absent from the file outranked every node present in it.
  The gap-fill leaf was the most common omission: `build` attaches it to the
  root, so it appears in no hierarchy file and had no obvious reason to appear
  in a priority file. Left out, it took precedence over every real feature, and
  a feature sharing all its k-mers with unannotated sequence stopped painting.

  Filling in a default does not avoid this. Zero wins outright, and a losing
  value fails differently: two omitted siblings would share it, so a group whose
  listed members are distinct would satisfy neither half of HKS's
  all-equal-or-all-distinct rule.

- **`prep-bed` writes the gap-fill leaf into the priority files it generates**,
  ranked after the features, so a k-mer a feature shares with unannotated
  sequence stays with the feature. Applies to `repeatmasker` and `edta`. The
  converters whose output already tiles declare no gap-fill and are unchanged.

- **Genome-scale temp files are created next to their output, not in the
  system tempdir.** The scaffold rewriters' spill directories (up to the whole
  assembly with unscaffolded contigs kept), the BAM→FASTA conversion, and the
  `hks smooth` intermediate TSV all previously landed in `/tmp` — on cluster
  nodes often a small or RAM-backed filesystem. They now sit beside the output,
  as annotate's smoothing pass already did.

- **A BAM input to `annotate` (HKS backend) is converted to FASTA once**,
  shared by every feature set's lookup, instead of once per feature set.

- **`hks` now writes KaryoScope's BED files directly, and the conversion pass
  is gone.** `convert_hks_tsv_to_bed` ran twice per feature set — on the raw
  lookup output and again inside `run_hks_smooth` — and never actually
  converted a format: the lookup TSV is already BED-shaped, so all it did was
  drop the header line and rewrite `none` to `novel`. `hks` 0.3.0 can be asked
  for both directly (`--miss-label`, `--no-header`), so both passes were
  deleted. Per feature set that removes **two ~5 GB reads and two ~5 GB
  writes**, on a pass measured between 9.5 s (fast local disk) and 21.0 s
  (cluster `/scratch`).

  It also **roughly halves peak temp disk** for an HKS run: the raw lookup
  output and the presmoothed BED were two near-identical multi-gigabyte copies
  of the same data, and now there is one file that `hks smooth` reads in place.
  A six-set HG002 run peaks at ~29 GB rather than ~34 GB, and the disk-space
  preflight no longer charges HKS runs for an intermediate that no longer
  exists.

  Verified byte-identical against the previous two-step pipeline on a real
  index, for both the presmoothed and the smoothed BED.
- **`annotate` requires `hks` >= 0.3.0**, checked before the run starts. `hks`
  has no `--version` flag, so KaryoScope reads the version it logs on startup.
  The check fails open — an unreadable version is not treated as too old — but
  a readable, too-old one stops the run with an upgrade message instead of
  letting it die partway through on `unexpected argument '--miss-label'`.
- **`annotate` reports per-phase timing, bytes and peak memory.** Each `hks`
  lookup and smooth logs its own wall time, output size and throughput, and
  peak memory is reported from `ru_maxrss` across the `hks` processes — where
  it actually lives. Previously one elapsed time per feature set covered both
  phases, which cannot attribute a change to either. `hks` 0.3.0 likewise
  times its base-index and labeling loads separately rather than reporting one
  combined figure.

### Fixed

- **BAM/CRAM conversion no longer depends on samtools' default flag filter.**
  The `samtools fasta` invocation now states `-F 0x900 -N` explicitly.
  - `-F 0x900` (secondary + supplementary excluded) was already the samtools
    default, but relying on a tool's default for a correctness-critical choice
    is how it changes underneath you. It is the right mask for this job:
    every read has exactly one primary record carrying full-length SEQ, whereas
    supplementary records are hard-clipped *slices* of reads the primary
    already supplied in full. Verified on a DRAGEN WGS CRAM — of 2,995,805
    primary records, zero were hard-clipped, zero lacked SEQ, and zero were
    under full read length. Unmapped reads (`0x4`) and duplicates (`0x400`) are
    *not* excluded: they are genuine distinct reads and dropping them would
    under-count.
  - `-N` forces the `/1`/`/2` mate suffix. This one is a real bug fix, not
    just an explicitness change: aligners routinely strip the suffix from
    QNAME, leaving both mates of a pair sharing a byte-identical name, and
    samtools only restores it where the READ1/READ2 flag bits happen to be set.
    Measured on the same data: with `-n`, 1,996,531 records collapsed onto
    1,002,139 distinct names — every mate pair silently indistinguishable. With
    `-N`, all 1,996,531 stayed distinct.

- **The `samtools | get_featureIDs` BAM pipeline can no longer deadlock on a
  warning-heavy BAM.** samtools stderr was read only after `get_featureIDs`
  finished; more than a pipe buffer's worth of warnings (malformed records
  warn once each) blocked samtools mid-write and stalled the pipeline. stderr
  now drains to a temp file.

### Documentation

- **An accuracy pass across the command pages, README and CONTRIBUTING.**
  `bin`'s `--threads 0` is documented as the CPUs the process may actually use
  (SLURM allocation / affinity), not `os.cpu_count()`; `scaffold`'s output
  filenames carry their `<input_stem>.<dbid>.` prefix; `version`'s tool list
  includes `get_featureIDs` and HKS; `build`'s options table lists
  `--flatten-order`; `annotate`'s HG002 runtime figure is attributed to the
  KMC backend; the recipes file table lists all six supporting files; and
  CONTRIBUTING's development setup uses `environment.yml`, which carries the
  tools the ad-hoc conda command it replaced was missing (samtools and the
  C++/Rust toolchains).

- **`build`: background placement and the hierarchy root.** The gap-fill leaf
  sits at the root, so a k-mer it shares with a real feature resolves to
  `categorized` and paints nothing. Covers when that is harmless (a set that
  already tiles; a background disjoint by construction, as `nonrepeat` is) and
  when it costs a set most of its resolution (a background holding homologs of
  the features), plus the fix: bring the leftovers into the tree under a shared
  parent instead of leaving them to gap-fill.

## [2.2.0] - 2026-08-04

### Changed
- **`--no-space-check` is now `--no-resource-check`** on `annotate` and
  `download`. The flag governs the memory check as well as the disk one, so
  "space" had become the wrong word. `--no-space-check` still works as a
  hidden alias and warns; it will be removed in a future release. Supplying
  both is a usage error.

### Added

- **`docs/recipes/` — reproducible recipes for the shipped databases.** One page
  each for `HKS_human_CHM13_v2` (`chromosome`, `region`, `repeat`, `gene`) and
  `HKS_arabidopsis_ColCEN` (all four sets), giving the exact download URL for
  every input, a SHA-256 for the file and for its decompressed content, the
  `prep-bed` command that converts it, a checksum for the resulting BED, and the
  build spec that assembles them.

  Inputs that cannot be derived ship alongside in `docs/recipes/files/`: the
  curated chromosome groupings, the repeat priority order, and the RefSeq
  accession map.

  Both recipes reproduce their databases. The Arabidopsis one has been run end
  to end from its published URLs and produces bit-identical index files; the
  CHM13 feature-set BEDs are byte-identical to those the shipped database was
  built from.

- **`karyoscope prep-bed`, converting source annotations into feature-set BEDs.**
  `build` starts from a final labelled BED; producing that BED previously needed
  a script per dataset. There is now one subcommand per *source format*:
  `repeatmasker` (native `.out` or the UCSC BED repackaging), `edta`,
  `gff-gene` (GFF3 or GTF), `cytoband` (UCSC golden-path or `cytoBandMapped`),
  `censat`, `fai`, and `satellite`. Formats are the unit rather than feature-set
  names because unrelated formats yield the same kind of set: RepeatMasker
  output and an EDTA GFF3 both produce a `repeat` set but share no parsing.

  Each writes its BED and hierarchy and prints the matching `feature_sets:`
  stanza on **stdout**, with progress and warnings on stderr, so the stanza can
  be appended straight to a build spec. Optional `--colors` writes the reference
  palette with `legend_group` filled in, so a cytoband set arrives with its
  legend already collapsed to nine stain rows.

  `prep-bed` does not gap-fill, flatten overlaps, or drop sequences: `build`
  owns `background:`, `flatten:` and `exclude:`. Sequences a set does not cover
  are reported for the spec's `exclude:` rather than given a placeholder
  `exclude` *label*. `censat` and `satellite` tile, because separating `p_arm`
  from `q_arm` needs the centromere position. A RepeatMasker class the converter
  has no leaf for is labelled `other_repeat`, keeping it distinct from
  RepeatMasker's own `Unknown` class.

  See [prep-bed](docs/commands/prep-bed.md).

- A dependency missing at *import* time now reports itself in a few readable
  lines naming the interpreter in use, instead of a bare `ModuleNotFoundError`
  traceback. `karyoscope.cli` eagerly imports every command module, so one
  missing package took down **every** command including `--help` and `version`
  — the two a user would reach for to diagnose it. That import happens inside
  the pip-generated console script, outside any code KaryoScope controls, so
  the fix required moving the entry point; the dependency preflight added in
  2.1.0 cannot help because it runs long after import.
- An explicit `--threads` above the usable CPU count now logs a warning naming
  the limiting source (e.g. `$SLURM_CPUS_PER_TASK`) and suggesting a value.
  Deliberately a warning, not a cap: oversubscription is sometimes faster on
  heterogeneous CPUs — an Apple M1 Max is 8 performance + 2 efficiency cores,
  and macOS exposes no way to learn that split, so "number of CPUs" is advice
  rather than a limit.

### Changed

- **`build` documentation restructured, and the variable-k description
  corrected.** The "mode A / mode B" naming is gone: BED plus a genome is
  described as the input, and the per-feature-FASTA form now sits in the
  build-spec section where it is actually usable. Hierarchy, priorities and
  colours each have a reference section that the options table links to, and
  overlap resolution is shown worked through — a hierarchy, a BED, and the label
  a shared k-mer ends up with, with and without priorities.

  Two claims were wrong. Priorities are applied when the index is built, not at
  query time. And `kmer.max_size` read as a range one may query within: on a
  fixed-k index it equals `kmer.size`, and that single length is the only one
  `annotate` accepts. Since HKS rejects priorities together with variable-k, a
  priority-resolved database is necessarily fixed-k — now stated in `build`,
  `annotate` and `info`.

- **CLI subcommands are imported on demand, cutting startup from ~290 ms to
  ~70 ms.** Registering the eleven subcommands meant importing all eleven
  command modules, and one of them (`download`) pulls in `requests` — about
  190 ms of that 290 ms, paid by every invocation regardless of whether the
  command used it. `karyoscope --help` still imports everything, because
  Click asks each command for its short help to build the listing; every other
  path now imports only the command being run.

  The broken-install check that `_entry` performs used to be a side-effect of
  those eager imports, so it is now explicit: the declared dependencies are
  verified with `importlib.util.find_spec`, which locates a module without
  executing it (~8 ms for the full set). The check is consequently *stricter*
  than before — it covers every declared dependency, where the old form caught
  only what the eager imports happened to touch — and the list is read from
  installed metadata so it cannot drift from `pyproject.toml`. Since
  `find_spec` proves presence but not importability, the CLI run is also
  wrapped, so a dependency that is installed yet unimportable (a truncated
  install, or `cairosvg` built against the wrong ABI) still reports readably
  instead of as a traceback.

- **`convert_hks_tsv_to_bed` reads binary blocks instead of looping per line.**
  It runs twice per feature set — on the raw lookup TSV and again inside
  `run_hks_smooth` — and on human input those two passes were ~10% of
  `annotate`'s wall time (129 s of a 21-minute HG002 run), single-threaded,
  while the rest of the machine idled. Now 8 MiB blocks through
  `bytes.replace`, with no decode/encode per line. **Measured 2.3x** on a 4 GB
  real `hks` TSV (21.9 s → 9.5 s), 1.7x at a 30% miss rate; output is
  byte-identical. Still streams — memory is one block plus a partial line.
  Blocks are cut at their last newline, so a match can never straddle a
  boundary (the token ends in a newline, and contains no interior one).
- **The `karyoscope` console script now routes through `karyoscope._entry`.**
  Existing editable installs need a `pip install -e .` for the new script to
  take effect; a fresh install is unaffected.

### Fixed

- **`colors.tsv` row order is now reproducible, and follows the hierarchy.**
  `build` wrote one row per node by iterating `Hierarchy.nodes()`, which returns
  a **set** — so Python's per-run string hash randomisation reordered the file,
  and rebuilding a database from byte-identical inputs produced a
  byte-different `colors.tsv`. Row order is emitted in hierarchy order instead
  (each root, then every child in edge order), which is also what makes legend
  groups meaningful: they are ordered by first appearance in `colors.tsv`, so
  for a cytoband set that ordering is now the Giemsa intensity progression
  rather than whatever the hash seed produced. Colour *assignments* are
  unchanged; only the order of the rows.

- **`HKS_arabidopsis_ColCEN` regenerated — its `gene` set labelled intergenic
  regions as intron.** The set had been built by deriving introns from a
  chromosome-wide exon list, which treats every gap between neighbouring genes
  as an intron: 8,325 such spans, reporting 42.0% of the genome as intron
  against 20.6%, and 21.6% intergenic against 43.0%. The database is rebuilt and
  now at 1.2.0. `HKS_human_CHM13_v2` was unaffected.

  `prep-bed gff-gene` derives introns from each transcript's own consecutive
  exons, and where transcripts disagree the more specific label wins
  (`exon` > `intron` > `intergenic`).

- **`colors.tsv` gains an optional 4th column, `legend_group`, and `build`
  carries it through.** A feature set with hundreds of leaves in a handful of
  colours produced a legend that dwarfed the figure and was truncated to
  whatever fit the canvas, with no indication that rows had been dropped. The
  CHM13 cytoband database has 833 features and the legend drew 51 of them.
  Features sharing a `legend_group` now collapse to one legend row: 833 entries
  become 9, labelled by Giemsa stain.

  The column is optional and the header must be exactly 3 or exactly 4 columns,
  so every existing database parses unchanged and keeps its per-feature legend
  (verified against `HKS_human_CHM13_v2`: 0 groups). `build --colors` accepts
  the column and writes it to the database's `colors.tsv`, so the file you
  supply and the file the database ships have the same shape; it is emitted only
  when at least one feature declares a group, leaving ungrouped builds
  byte-identical. Because a legend row is one swatch and one label, `build`
  fails if a group spans two colours, since there is no well-defined swatch for
  such a group. See
  [build → Grouping the legend](docs/commands/build.md#grouping-the-legend).

- **`examples/karyotypes/`** — reference karyotype plots for six assemblies
  (CHM13, HG002, the HG008 tumour/normal pair, an HPRC population sample, and
  Arabidopsis) with notes on what each shows, for comparing your own output
  against.

- **Karyotype outlines now follow the theme, and the legend lists only what is
  visible.** The sequence-outline guard read `if background_color == "white"`,
  so a black-background plot got no border at all — and the cytoband palette
  contains pure `#000000` (the gpos100 bands), which then merged into the
  backdrop and vanished; 159 fills in a real HG002 plot were affected. Legend
  swatches were stroked black unconditionally, so on a dark background a
  black-filled swatch disappeared entirely. Outlines are now drawn on every
  background and derive from the same `text_color` the module already computed;
  there are zero hardcoded `stroke="black"` left.

  Separately, the legend now omits features whose entire drawn extent is under
  half a pixel, since a legend row for such a feature names a colour that is not
  visibly present in the figure.

- **`bin` no longer emits a runt trailing bin.** A trailing partial bin cast a
  label vote of equal standing to a full one, so a handful of bases could
  outvote hundreds of kilobases. This was not a rare edge case: the karyotype
  genome view chooses `bin_size = round(longest_sequence / 250)`, which by
  construction leaves the longest sequence with a remainder of order tens of
  bases every time. On CHM13 it produced a lone **48 bp `categorized` row** at
  the end of chr1 — the last k-mer starts before the trailing `k-1` bases,
  which stay ambiguous up to the hierarchy root — and that row then earned a
  full karyotype legend entry for something occupying about 1/5000 of a pixel.
  A trailing bin shorter than half `bin_size` is now folded into the preceding
  row. No bases are lost, and a runt that is the only bin on its contig is
  still emitted rather than dropped.

- **A failed `build` no longer blocks the re-run.** It left the
  partially-written database directory behind, so the obvious next move — fix
  the offending BED, colours, or priority file and run the same command again —
  hit `database directory already exists ... Pass --force to overwrite`, about
  a directory that never contained a working database. A build that fails now
  removes what it created, restoring the state it started from.
  `--keep-intermediates` retains it for inspecting a build that failed late.
  (Note that `--force` still deletes an existing database *before* building, so
  a failed `--force` rebuild does not bring the old one back.)

- **`--threads 0` (the default) now sizes its worker pool from the CPUs the
  process may actually use**, not the machine's core count. It read
  `os.cpu_count()`, which ignores cgroups, CPU affinity, and SLURM allocations.
  Measured on one of our own nodes: `os.cpu_count()` reported 36 while
  `sched_getaffinity` reported 1, so the default spawned 36 workers to contend
  for a single allocated CPU. New `karyoscope.cpus` resolves, in order,
  `$SLURM_CPUS_PER_TASK`, `sched_getaffinity`, `os.process_cpu_count` (3.13+),
  then `os.cpu_count`.

## [2.1.0] - 2026-07-27

### Added
- **Progress output for the long-running commands.** `download` and `build`
  already announced themselves and their results, but `annotate` (7-22 min)
  and `karyotype` (tens of minutes) printed nothing at all until their closing
  `Wrote:` block — so the two slowest commands were the two silent ones, and a
  user could not tell a working run from a hung one. Both now print an opening
  summary and a milestone line as each unit of work completes:

  ```
  Annotating hg002v1.1.fasta.gz against HKS_human_CHM13_v2
    6 feature set(s), 16 thread(s), ~34 GB estimated output
    [1/6] chromosome    4m05s
    ...
    bgzip (12 file(s))  1m31s
  Wrote:
    ...
  ```

  The shape follows what each backend actually knows: HKS reports per feature
  set, while KMC reports named phases (`k-mer query`, `smoothing 6 feature
  set(s)`) because it smooths every set in a single streaming pass. `karyotype`
  reports per `mode/feature_set` render and indents the cascade's nested
  `annotate` one level, so its headline reads as a step of the run rather than
  a separate command. It also names any feature set the cascade pulls in
  beyond `--feature-set` — scaffolding needs the chromosome- and
  region-assignment sets to place and orient contigs regardless of what is
  being plotted, so asking for two sets and watching `annotate` report three
  was otherwise baffling. Detailed per-step timings stay at INFO, so `-v` is
  unchanged and nothing is printed twice.
- `-q/--quiet` now suppresses this narration as well as logging. It previously
  only lowered the log level, which would have left no way to silence a run.
  The `Wrote:` block still prints — it is the command's result, not narration.
- **Free-space preflight.** `download` and `annotate` now verify the target
  filesystem can hold what they are about to write, and fail immediately with
  the required / available / shortfall figures instead of dying on
  `OSError: [Errno 28]` after the work is done. Previously a user with 12 GB
  free could spend 25 minutes downloading `HKS_human_CHM13_v2` before extraction
  ran out of room, and a six-feature-set `annotate` of a diploid assembly could
  fill a disk 20 minutes in. Both checks are overridable with `--no-space-check`.
  A new `karyoscope.diskspace` module also translates any `ENOSPC` that escapes
  mid-run — from any command — into a message naming the filesystem that filled
  up, rather than a traceback.
- **Dependency preflight.** Commands now resolve the external tools they will
  need before starting, and report *every* missing one at once with an install
  hint, instead of failing at the point of use. `annotate` checks its k-mer
  backend binary (plus `samtools` for BAM input and `bgzip` unless `--no-bgzip`);
  `karyotype` additionally checks `cairosvg` when `--format pdf/png` is requested
  and `seqtk` when telomere detection will be auto-run — both of which used to
  surface only after the entire cascade had already run. Resolution goes through
  each backend's own lookup order, so `$KARYOSCOPE_HKS`,
  `$KARYOSCOPE_GET_FEATUREIDS`, and the source-tree `get_featureIDs` of an
  editable install are all honoured.
- Registry entries may declare `download_size_gb` (the `.tar.gz`) alongside
  `size_gb`, which is now defined as the size of the *extracted* database.
  Installing needs the sum of the two, since the archive is not deleted until
  extraction succeeds — ~34 GB for `KS_human_CHM13_v2` and ~36 GB for
  `HKS_human_CHM13_v2`. Entries without `download_size_gb` fall back to
  `size_gb` for both and are labelled as estimates.

### Fixed
- `download` no longer discards a completed archive when the install fails.
  The archive was unlinked in a `finally`, so a run that finished a 25-minute
  transfer and then hit ENOSPC during extraction left nothing behind and the
  retry re-downloaded everything. It is now kept (its SHA-256 proves the bytes
  are identical to a fresh fetch), re-verified on the next run, and extracted
  directly — a failed install costs one extraction, not a second transfer.
  Conversely, a partially-extracted database directory is now *removed* on
  failure: it is unusable, is deleted at the start of the next attempt anyway,
  and when the failure was ENOSPC it occupies exactly the space the retry
  needs. Both outcomes are reported at WARNING so they are visible by default.
  A staged archive whose checksum doesn't match is discarded and re-fetched.

### Changed
- `download --list` and `--info` report the download size and the on-disk size
  separately, and `--info` also reports the free space needed to install. The
  single unlabelled size they used to print was the archive size for
  `HKS_human_CHM13_v2` and the extracted size for `KS_human_CHM13_v2`, so it
  understated the requirement by ~9 GB for the former.
- `karyoscope version` reports `get_featureIDs` (previously omitted entirely)
  and resolves it and `hks` the same way the commands do, so an editable install
  or an environment override is no longer reported as "not found on PATH".

### Documentation
- README and the `download` / `annotate` command pages gained Disk space
  sections covering archive vs extracted database size, the peak during install,
  and annotation output footprint (~0.8 GB per feature set per Gbp of input).
  Documents that `--bgzip` shrinks the result but not the peak, because
  compression runs only after every BED has been written.

## [2.0.0] - 2026-07-23

### Added
- **HKS index backend** alongside the existing KMC backend. Databases declare
  `index.type: hks` in their manifest; `annotate` dispatches to it, querying with
  `hks lookup` and smoothing with `hks smooth` (both threaded via `--threads`).
  HKS (Hierarchical K-mer Sets) resolves multiply-labelled k-mers through the
  label hierarchy at query time. On the CHM13 feature sets it runs ~2.5-3x faster
  and at ~1/3 the peak RAM of the KMC backend. Ships the new `HKS_human_CHM13_v2`
  database, which re-indexes the same k-mers as the KMC `KS_human_CHM13_v2` (the
  KMC database is unchanged). `features.tsv` is not needed by this backend (see
  Changed). Requires the `hks` binary on `PATH` (or `$KARYOSCOPE_HKS`).
- `environment.yml` for the suite's conda environment (includes `rust`, needed to
  build the `hks` binary from source), shared with KaryoScope-ISCN.
- `build --exclude SEQID` (and spec `exclude:`): drop sequences (e.g. organelles
  `ChrM`/`ChrC`) from the whole build — removed from every feature BED and the
  gap-fill index, so no set covers them and they read as `none` everywhere
  (uniform across sets) and never appear as karyotype chromosomes. The
  `chromosome` feature set thus declares the karyotype chromosomes.
- `karyotype` / `scaffold` / `centromeres` `--telo-motif`: telomere repeat motif
  for the auto-run `seqtk telo` (its `-m`). Default is seqtk's `CCCTAA` (vertebrate
  `TTAGGG`); non-vertebrates need their own (e.g. Arabidopsis/plants `CCCTAAA`,
  `TTTAGGG`), for which the human default detects nothing.
- `karyotype --pixels-per-mb`: pin the vertical zoom (px per Mb) to compare plots
  across assemblies at a fixed scale.
- `build` command: construct a complete, registry-ready HKS database from a
  genome and per-feature-set BED annotations (4th column = leaf label), then
  register it. Runs the HKS `build-base` / `add-feature-set` steps, gap-fills
  unannotated regions with a named `background` leaf (distinct from HKS's `none`
  novel sentinel), derives the label hierarchy (flat star by default) and
  colours (auto palette by default), and writes `manifest.yaml` / `hierarchy.tsv`
  / `colors.tsv`. Overlapping BEDs are allowed — HKS resolves multiply-labelled
  k-mers via its hierarchy, and an optional per-set priority file makes the
  higher-priority label win (the per-k-mer equivalent of pre-flattening; a
  `--flatten` fallback remains). Two entry forms feed one pipeline: the simple
  `--id`/`--sequence`/`--feature-set NAME=bed` flags, or a `--spec build.yaml`
  file for multi-feature-set databases. `--variable-k` builds an index queryable
  at any k ≤ s from a single build (e.g. a k-sweep to validate the k=31 default),
  supported even from BED input by building the base from the generated
  per-feature FASTAs; the manifest's `kmer.type` becomes `variable`. New
  per-command docs at `docs/commands/build.md`.
- `annotate --k INTEGER`: query at a chosen k-mer length. Defaults to the
  database's k; a value other than the database's k is accepted only on a
  variable-k HKS index, enabling a k-sweep (e.g. re-annotating the same input at
  k=21/25/31 from one index). Outputs are tagged `.k<k>` so a sweep into one
  directory doesn't overwrite itself; fixed-k indexes reject a differing `--k`
  with a message pointing at `build --variable-k`.
- `karyotype --colors PATH`: render with a custom colour file (same
  `feature_set`/`feature`/`color` format as a database `colors.tsv`) instead of
  the database default — e.g. the cytoband database's `colors_chromosome.tsv` to
  colour bands by chromosome. The colour file's stem is appended to the output
  filename (`...smoothed.<colors-stem>.karyotype.svg`) so a custom-colour render
  never clobbers the default-colour one in the same directory.
- `remap-bed` command: apply an existing `scaffold_map.tsv` to a BED that was
  annotated separately -- possibly against a *different* database than the one
  used to derive the map -- rewriting it into the scaffolded coordinate system.
  This is the standalone, file-producing counterpart of the in-pipeline rewrite
  `scaffold` performs, and of the in-memory remap `karyotype --scaffold-db`
  does while rendering. It enables a two-database workflow (lay out with a
  roles-bearing database, annotate a feature set from a plot-only database such
  as a cytoband database, then remap onto the layout). Before rewriting it
  validates that the BED and map describe the same assembly: hard errors on
  zero contig-name overlap and on BED intervals exceeding the contig lengths
  recorded in the map; an advisory warning (promotable via `--strict`) on a
  filename-stem mismatch. Wraps the shared
  `karyoscope.core.scaffold.remap_bed_with_map` helper.
- Zenodo archive DOI for the software (concept DOI
  ``10.5281/zenodo.20657816``) added to the README badge, the citation
  section, and ``CITATION.cff`` (top-level ``doi`` / ``version`` /
  ``date-released`` fields).

### Fixed
- `annotate`'s "no databases installed" error now mentions `--db-root`, making it
  clear the databases may simply be under a non-default root.
- The default karyotype output filename for multi-input runs (no `--output`) is
  no longer named after only the first input. It now collapses the input stems to
  their common prefix (`GM04890.haplotype1` + `GM04890.haplotype2` -> `GM04890`),
  so a both-haplotype plot is no longer misleadingly named `...haplotype1...`. The
  title band already showed both stems; the filename now agrees. Single-input runs
  are unchanged, and runs whose stems share no separator-delimited prefix fall
  back to the first stem. (Also corrected the stale module docstring, which claimed
  the default base was the literal `"karyotype"`.)
- Genome-karyotype haplotype columns are now assigned and ordered by the *true*
  haplotype encoded in each contig's name (`infer_hap_from_contig`) rather than
  the file-level `hap` label recorded during scaffolding. Combined-FASTA
  assemblies (one file, hap-tagged contig names such as `haplotype1-*` /
  `haplotype2-*`) are labelled with a single file-level hap by `scaffold`, which
  previously collapsed every contig into one column drawn by contig size -- so a
  larger hap2 contig could read as hap1. The renderer now derives the column
  from the contig name, keeping `hap1` left of `hap2` regardless of size, and
  labels each column (`h1` / `h2`). Surfaced by ISCN validation on GM00392
  chr16.
- Unassigned contigs (`*_unassigned-*` fragments in combined-FASTA assemblies)
  are now segregated into their own labelled (`u`) column at the right of each
  chromosome, instead of being silently mixed into a haplotype column. Makes a
  small unassigned fragment (e.g. GM00392 `chr2_*_unassigned-*`, a 1-bin piece)
  visually distinguishable from a real haplotype sequence.
- Acrocentric short-arm contigs (chr13/14/15/21/22 p-arms: satellite, stalk,
  rDNA, plus the p-ter telomere) are now consistently oriented telomere-first
  (p-ter at the top) during scaffolding. Previously the flip decision for these
  fragments was gated on the p-ter telomere being *continuous* with the
  chromosome-assignment region block; satellite/novel gaps routinely break that
  continuity, so otherwise-identical short arms were oriented inconsistently
  (some ``CT``, some ``TC``). The acrocentric short-arm rule now keys off the raw
  telomere flags instead, flipping a single-telomere non-q-arm fragment so its
  telomere sits at the top regardless of continuity. Full-chromosome and q-arm
  body contigs (with both telomeres or a q-ter telomere) are unaffected.
  Surfaced by ISCN validation on GM03417.

### Changed
- With a second reference database now shipped (`HKS_human_CHM13_v2` alongside the
  KMC `KS_human_CHM13_v2`), `annotate` / `karyotype` / `scaffold` require `--db`
  when more than one database is installed, erroring with the list of installed
  databases. A single-database install still selects it automatically, and
  existing annotations are reused per database (every intermediate filename is
  tagged with its database id), so upgrading and switching backends never mixes them.
- `karyotype` rendering now scales to the data instead of fixed human-tuned
  constants: the longest chromosome (genome view) / longest centromere (centromere
  view) fills a fixed pixel height, the genome-view bin size is derived from the
  longest sequence (≈ longest / 250, restoring feature diversity that a fixed 1 Mb
  bin washed out on small genomes), and the scale bar is a "nice" round length for
  the zoom. Human output is unchanged; small genomes (e.g. Arabidopsis) no longer
  render tiny. Override with `--pixels-per-mb` / `--bin-size`. The centromere view
  is also taller (fills the plot like the genome view).
- `karyotype --no-human-chroms` / layout seeding: the karyotype layout is now
  seeded from the database's own `chromosome` feature-set leaves rather than a
  hardcoded human list (identical for the human database; organism-correct
  elsewhere). `--no-human-chroms` now suppresses that seeding.
- `karyotype` title band auto-widens the canvas so a long title over few
  chromosomes is no longer clipped.
- Memory: the `bin`, `scaffold`, and `centromeres` stages and the HKS TSV->BED
  conversion now stream per contig / in bounded batches instead of loading whole
  assemblies or files into memory. Peak RAM on a whole-genome cascade drops from
  ~32-46 GB to ~10 GB (now bounded by `annotate`), so runs fit on much smaller
  nodes; outputs are byte-identical. Read inputs also emit integer query ranks
  rather than names, and a redundant scaffold pass was dropped.
- `features.tsv` is now optional for HKS databases (`index.type: hks`). It maps
  integer feature ids to names for the KMC backend only; the HKS backend reads
  label names from the index and never consults it. The manifest parser and
  layout validator no longer require it for HKS (KMC still does), `annotate`
  tolerates its absence, and databases built by `build` omit it. Existing HKS
  databases that still list `features.tsv` continue to validate.
- With ``--combine-chromosomes`` but not ``--combine-acrocentrics``, the
  acrocentric contigs that stay as separate records are now renamed in
  canonical order to ``<chrom>_<hap>_<A|B|C...>`` (e.g. ``chr14_hap2_A``,
  ``chr14_hap2_B``), dropping the original contig name and any ``_rc``
  suffix -- matching how the combined ``<chrom>_<hap>`` records are named.
  Previously these contigs kept their encoded ``<chrom>_<hap>_<contig>[_rc]``
  names. Affects the combined FASTA, BED, and AGP outputs of both
  ``scaffold`` and ``karyotype``; the per-contig scaffold map TSV is
  unchanged.

## [1.1.0] - 2026-06-11

### Added
- Per-command reference documentation under ``docs/commands/`` (one page
  per subcommand: ``download``, ``register``, ``annotate``, ``scaffold``,
  ``bin``, ``centromeres``, ``karyotype``, ``info``, ``version``), linked
  from the README's Commands table (those links were previously dead).
- README ``System requirements`` section (operating systems, Python and
  dependency versions, external tools, and hardware/RAM guidance), a
  typical install-time note, and a ``Demo`` section that runs the
  annotation pipeline end to end in seconds against the bundled synthetic
  database -- addressing the Nature code-submission checklist.
- ``examples/`` directory with ``run_demo.sh`` (a self-contained install
  smoke test that annotates ``examples/demo.fa`` in a throwaway database
  root and cleans up after itself), ``demo.fa``, and a README. The tiny
  demo FASTA is committed via a ``.gitignore`` exception.
- ``karyoscope scaffold`` can now concatenate the contigs of each
  chromosome+haplotype into a single sequence, producing a less
  fragmented assembly. New flags (FASTA output only, i.e. ``--mode
  fasta`` / ``--mode both``):
  - ``--combine-chromosomes`` (default off): join all contigs of one
    ``(chromosome, hap)`` into one ``<chrom>_<hap>`` record.
  - ``--scaffold-gap-size`` (default 100000): number of ``N`` bases
    inserted between concatenated contigs.
  - ``--combine-acrocentrics`` (default off): also combine the
    acrocentric chromosomes. Off by default because acrocentric p-arms
    recombine and KaryoScope's chromosome assignment there is less
    certain, so their contigs stay as separate records unless this is
    set. The acrocentric set is the same one ``--acrocentric`` controls.

  Combined outputs carry a ``combined_chromosomes`` filename tag so a
  user can run with and without combining in the same directory:
  ``<stem>.<db>.scaffolded.combined_chromosomes.fa[.gz]`` and, in
  ``--mode both``, ``<stem>.<db>.<fs>.smoothed.scaffolded.combined_chromosomes.bed[.gz]``
  (the plain per-contig scaffolded BEDs are not written in a combine
  run). An AGP 2.1 file,
  ``<stem>.<db>.scaffolded.combined_chromosomes.agp``, documents every
  component placement and gap (``gap_type=scaffold``, ``linkage=yes``,
  ``linkage_evidence=align_genus`` -- KaryoScope is alignment-free, but
  the contig order/orientation is asserted from each contig's k-mer
  feature profile against a same-genus human reference). The AGP fully
  describes the output FASTA, including kept unscaffolded leftovers as
  singleton objects.

  The combined BEDs share the FASTA's coordinate system exactly. Each
  per-contig annotation BED tiles ``[0, E)`` where ``E = L - k + 1``
  (``L`` = true contig length, ``k`` = k-mer size); the map's ``length``
  field is ``E``. When contigs are concatenated as ``seq + N*gap + seq``,
  each contig's intervals shift by its true-base offset and a ``novel``
  interval fills ``[offset + E, next_offset)`` between contigs -- the
  literal N gap **plus** the contig's own untiled ``k-1`` tail, whose
  k-mers in the concatenated assembly overlap the Ns. This is
  byte-identical to re-annotating the combined FASTA (smoothing never
  bridges the gap, which is far larger than ``max_gap=1000``). The
  implementation reads true lengths from the FASTA and tiling ends from
  each BED, so it hardcodes no ``k`` and stays correct for a future
  variable-k database (designed-for but, lacking such a database, not
  yet exercised in tests).
- ``karyoscope karyotype`` gained ``--combine-chromosomes`` (plus
  ``--scaffold-gap-size`` and ``--combine-acrocentrics``, mirroring the
  ``scaffold`` flags) to render combined-chromosome karyotypes. When
  set, karyotype cascades ``scaffold`` in ``--mode both`` with
  ``--combine-chromosomes`` so the ``combined_chromosomes`` BEDs (and
  the combined FASTA + AGP side artifacts) are created, then lays out
  the figure from those combined BEDs: each ``<chrom>_<hap>`` is one
  ideogram cell instead of one cell per contig. The renderer is
  unchanged -- karyotype bins the combined BEDs (keyed ``<chrom>_<hap>``)
  and joins them to synthetic combined map rows (one per object, with
  end-telomere flags taken from the first and last component). All three
  modes work: ``genome`` / ``subtelomere`` consume the combined binned
  BEDs directly; ``centromere`` detects centromere ranges in the
  combined coordinate system (keyed ``<chrom>_<hap>``) and writes a
  ``<stem>.<db>.centromeres.combined_chromosomes.bed[.gz]``. The
  SVG / PDF / PNG filenames carry the ``combined_chromosomes`` tag so
  combined and per-contig figures coexist in one directory. Acrocentric
  groups left uncombined (the default) still render as their per-contig
  cells.
- ``karyoscope annotate`` is now resumable across the expensive
  ``get_featureIDs`` (k-mer query) step. On success that step writes a
  small ``<input>.<dbid>.combined.presmoothed.featureIDs.bed.done``
  completion marker recording the combined BED's size and mtime. On a
  rerun, if a complete-and-verified combined BED is already present,
  ``annotate`` skips ``get_featureIDs`` and resumes straight into the
  smoothing pass -- the common recovery path after smoothing workers
  are OOM-killed: rerun with fewer ``--threads`` (or more RAM) and the
  k-mer query is not repeated. Pass the new ``--force`` flag to
  regenerate the intermediate unconditionally.

  A combined BED left behind by a *killed* run (truncated / partial)
  has no matching marker and is **never** silently reused -- it is
  regenerated -- so resume cannot produce wrong annotations from a
  half-written intermediate. Previously a rerun unconditionally re-ran
  ``get_featureIDs``, overwriting an already-completed combined BED in
  place (and, if that rerun was itself OOM-killed mid-write, destroying
  the good file).
- ``karyoscope karyotype`` now accepts ``--smoothed/--presmoothed``
  (default ``--smoothed``). ``--presmoothed`` renders from the raw
  (unsmoothed) annotation BEDs instead of the hierarchy-smoothed
  ones, so users can visually compare the effect of smoothing. Both
  variants are produced by ``karyoscope annotate``.
- ``karyoscope register`` registers a database that is already present
  under the database root (built locally or copied from another machine)
  by writing its ``installed.json`` entry, so the data commands can use
  it without going through ``download``. It validates the layout, derives
  the id and version from ``manifest.yaml``, records ``source: local``,
  and refuses to clobber an existing entry without ``--force``.

### Fixed
- Haplotype inference now recognises contig names of the form
  ``haplotype1...`` / ``haplotype2...`` (e.g. ``haplotype1-0000001``,
  as emitted by some long-read assemblers) and maps them to
  ``hap1`` / ``hap2``. Previously these matched no built-in pattern, so
  a single combined dual-haplotype FASTA collapsed to one haplotype and
  the karyotype drew no h1/h2 split. The trailing boundary keeps
  ``haplotype10`` from being misread as hap1.
- The karyotype no longer reuses a **stale binned-scaffolded BED** when
  the scaffold map has changed since that BED was built (for example
  after the haplotype-inference fix above moves a contig from ``hap1``
  to ``hap2``). The binned-scaffolded BED bakes the map's rename and
  orientation into its contents, so each one now records a signature of
  the scaffold map it was built from in a ``.mapsig`` sidecar; on a
  later run the binned BED is rebuilt when that signature no longer
  matches the current map (with ``--no-auto``, a clear error is raised
  instead of silently rendering a stale layout). Without this, a user
  picking up the inference fix would have had to manually delete the
  derived ``*.scaffolded.binned*.bed.gz`` files to see the corrected
  karyotype.
- The smoothing pass no longer floods stderr with benign
  ``BrokenPipeError`` / ``EOFError`` tracebacks when a worker is
  OOM-killed. Those came from ``multiprocessing.Pool``'s daemon helper
  threads writing to the dead worker's pipe in the window before the
  watchdog fired; they are now suppressed for the duration of the
  smoothing pass so the watchdog's single actionable ``FATAL`` message
  stands on its own.

### Changed
- ``Development Status`` trove classifier in ``pyproject.toml`` updated
  from ``3 - Alpha`` to ``5 - Production/Stable`` to reflect the stable
  release line, the README's pre-1.0 "expect breaking changes" banner
  replaced with an accurate semantic-versioning statement, and the CI
  status badge re-enabled now that the org's GitHub Actions minutes have
  refreshed (runs are green again on ``main``).
- Karyotype SVG filenames now include the annotation variant
  (``smoothed`` or ``presmoothed``) so both can coexist on disk:
  ``<stem>.<dbid>.<mode>.<fs>.smoothed.karyotype.svg`` (was
  ``<stem>.<dbid>.<mode>.<fs>.karyotype.svg``). All intermediate
  BED filenames (scaffolded, binned-scaffolded) also include the
  variant. This is a naming change from v1.0.0; existing scripts
  that glob for ``*.karyotype.svg`` will still match.

### Deprecated
- The database-root override on ``karyoscope info`` and ``karyoscope
  download`` is now spelled ``--db-root``, matching the data commands
  (``annotate``, ``bin``, ``scaffold``, ``centromeres``, ``karyotype``),
  where ``--db`` selects a database *id*. ``--db`` remains a hidden,
  working alias on ``info``/``download`` for one release and prints a
  deprecation warning; it will be removed in a future version. Switch to
  ``--db-root``.

## [1.0.0] - 2026-05-21

### Added
- ``karyoscope download`` now validates registry-entry metadata up
  front (``url``, ``sha256``, ``karyoscope_min_version``) and fires a
  clean, actionable error if the entry is incomplete -- before any
  network or filesystem work. Catches three real-world cases:

    * ``url`` is missing, empty, ``"PLACEHOLDER"``, or doesn't have
      a known URL scheme (``http://``, ``https://``, ``file://``).
      Previously the user saw an opaque ``urllib`` error like
      ``unknown url type: 'PLACEHOLDER'``; now they see a message
      pointing at the KaryoScope-registry repo for publication
      status.
    * ``sha256`` is missing, ``"PLACEHOLDER"``, or not a 64-hex
      digest. Skipped when ``--no-checksum`` was passed.
    * ``karyoscope_min_version`` is missing or ``"PLACEHOLDER"``.
      Previously these silently parsed to ``(0,)`` in
      ``_check_version_compatibility`` and bypassed the compat
      guard entirely; now the hygiene check catches them.

  Critically, the hygiene check fires BEFORE the destructive
  ``shutil.rmtree`` step in ``install_database``, so a malformed
  registry entry can never destroy an existing install just
  because the URL turned out to be bogus. Validates at download
  time (not parse time) so ``karyoscope download --list`` and
  ``--info`` still work for in-progress registry entries -- users
  can see what's coming even if it's not yet downloadable. 33
  new unit tests in ``test_download_core.py``.
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
  and stdin. Adapted from an earlier internal version; the only changes
  are include paths driven by the new Makefile's `-I` flags rather than
  the previous `../external/...` relative paths.
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

### Added
- ``karyoscope annotate`` now appends an actionable hint to the error
  message when ``get_featureIDs`` exits with SIGKILL (-9 or 137).
  SIGKILL on the k-mer query step is overwhelmingly "the kernel
  OOM-killer or the job scheduler killed the process for using too
  much memory" -- the previous bare ``Error: command failed with exit
  code -9`` was opaque to non-cluster-experienced users. The hint
  spells out the likely cause (KMC index needs ~20-30 GB to load),
  recommended fixes (request more memory on SLURM, limit threads, move
  off login nodes), and the threading footgun (default ``-t 0`` auto-
  detects the machine's full core count, which on shared nodes can be
  much higher than the memory allocation supports). Applied to both
  the file-input and BAM-pipe paths via a shared
  ``_augment_with_oom_hint`` helper that only fires for SIGKILL-like
  exit codes -- non-OOM failures (malformed input, missing files, etc.)
  pass through with no hint to avoid misdirection. 11 new unit tests in
  ``test_kmc_wrapper.py``.
- ``karyoscope annotate`` smoothing now has **three** dispatch paths
  selected automatically based on ``--preserve-order`` and the input
  file extension:

    * **Assembly (FASTA) + ``--preserve-order``** (default): per-(fs, seq)
      temp files + concat in input order. Unchanged. Essential for
      whole-chromosome chunks that produce hundreds of MB of BED
      lines per FS that would overload the IPC pipe.
    * **Reads (FASTQ/BAM) + ``--preserve-order``** (NEW): streaming
      dispatch via ``Pool.imap`` (ordered), workers return BED lines
      via IPC, main writes in input order. No per-sequence temp
      files -- per-read temp files would scale catastrophically
      (millions of tiny files, file-descriptor exhaustion). Safe
      because read chunks are uniform-size, so the ordered iterator
      doesn't stall waiting for a single slow chunk and per-chunk
      IPC payloads stay small (a few MB).
    * **Anything + ``--no-preserve-order``**: streaming dispatch via
      ``Pool.imap_unordered``. Fastest path when output order doesn't
      matter downstream.

  Detection is by extension (``.fastq``, ``.fq``, ``.fastq.gz``,
  ``.fq.gz``, ``.bam``); long-read FASTA users with millions of
  sequences should pass ``--no-preserve-order`` explicitly to opt
  out of the per-sequence temp-files path.

  Internally, the previous ``_smooth_streaming_unordered`` helper was
  generalised to ``_smooth_streaming(*, ordered: bool, ...)`` so both
  the ordered (reads + preserve) and unordered (any + no-preserve)
  paths share the same body modulo the choice of ``Pool.imap`` vs
  ``Pool.imap_unordered``.
- ``karyoscope annotate`` now accepts **FASTQ** (``.fastq``, ``.fq``,
  ``.fastq.gz``, ``.fq.gz``) and **BAM** (``.bam``) inputs in addition
  to FASTA. FASTQ is read by the C++ ``get_featureIDs`` binary natively;
  BAM is streamed through ``samtools fasta`` (not ``samtools fastq`` --
  KaryoScope only needs the sequence, not the quality string, and FASTA
  is smaller and slightly faster to write). The BAM path is a streaming
  pipe -- no intermediate file is created and memory stays bounded.
  ``samtools`` is required on ``$PATH`` for BAM inputs; the wrapper
  raises ``ToolNotFoundError`` with an actionable install hint when
  missing. For read-level inputs (many short sequences), also pass
  ``--no-preserve-order`` to skip the per-sequence temp-file machinery
  designed for assembly-scale inputs.
- ``karyoscope scaffold`` / ``karyotype`` / ``centromeres`` now reject
  FASTQ, BAM, and SAM inputs up front with a clean error message
  pointing the user at ``karyoscope annotate``. These pipelines need
  contig names from a FASTA; reads have no chromosome-scale "contig"
  concept, so scaffolding and karyotype rendering don't apply.
  Previously a FASTQ/BAM input would silently produce zero contigs.

### Changed
- ``karyoscope karyotype`` legend sort now pins chromosomes (``chr*``)
  to the very top in natural order (chr1, chr2, ..., chr22, then
  alphabetical for chrX/chrY/chrM), then ``"categorized"`` (the
  hierarchy root), then any hierarchy-order categorical groupings
  (e.g. ``autosome``, ``acrocentric``, ``metacentric``,
  ``submetacentric``, ``sex`` for the production CHM13 chromosome
  feature set), then unranked features alphabetical, then ``"novel"``
  at the very bottom. Previously ``"categorized"`` pinned above the
  chromosomes, which read awkwardly for the chromosome feature set
  -- the parent category appeared before the leaves it parents.
  Helper extracted to a module-level ``_legend_sort_key`` for
  testability; 6 new unit tests cover the layout including the
  production chromosome FS, the "future database with unranked
  features" fallback, and the chr-natural-order subcase.
- ``karyoscope karyotype`` now orders ``paternal`` before ``maternal``
  in the rendered hap columns (matching HPRC's convention of
  ``hap1 = paternal``, ``hap2 = maternal``). Previously the alphabetical
  sort put ``maternal`` first, which placed the maternal column on the
  left of every chromosome cell -- backwards from the standard HPRC
  layout. Numeric ``hapN`` labels still sort numerically (``hap1`` →
  ``hap2`` → ``hap10``), other labels still sort alphabetically,
  ``unassigned`` is still last.
- ``karyoscope karyotype`` legend text now renders at 14 pt (was
  11 pt) to match the chromosome and hap labels. Swatch (14 px) and
  row height (20 px) scaled in proportion; legend-band width estimate
  updated to 8 px/char (was 6 px/char) so the canvas stays tight
  against the longest label.
- ``karyoscope karyotype`` heterogametic-hap detection (which hap
  holds the heterogametic chromosome, e.g. chrY in XY male) now
  infers from data first, falling back to the sort-order convention.
  The archive's logic assumed ``haplotypes[0]`` was the heterogametic
  hap -- correct for ``hap1``/``hap2`` labels where ``hap1`` holds
  chrY by convention, but **backwards for the conventional
  ``maternal``/``paternal`` labels** (chrY is paternal, not maternal).
  The data-driven inference looks at which hap actually has chrY
  records in the scaffold map and uses that as the heterogametic
  hap; ``chrX`` (homogametic chromosome in male) gets the other hap.
  Falls back to ``haplotypes[0]`` only when inference is ambiguous
  (no chrY data, e.g. cancer with chrY loss; or chrY present in
  multiple haps from mis-labelled contigs). Combined with the
  HPRC-aligned sort order above, the fallback now also produces the
  biologically correct labelling for the conventional human case.
  Fixes the rendering bug where HG002 (male) drew empty extra
  columns on chrX and chrY with the wrong hap labels.
- ``karyoscope karyotype`` now defaults its output filename base to the
  first ``--input``'s stem (e.g. ``hg002v1.1.<dbid>.<mode>.<fs>.karyotype.svg``)
  instead of the previous literal ``karyotype`` (e.g.
  ``karyotype.<dbid>.<mode>.<fs>.karyotype.svg``). Tells the user at a
  glance which sample produced the SVGs. Pass ``--output PATH`` to
  override the base name explicitly. Multi-input runs still get the
  first input's stem; the user can disambiguate with ``--output`` if
  needed.
- ``karyoscope karyotype`` now exposes ``--scaffolding/--no-scaffolding``
  (default on). With ``--no-scaffolding`` the cascade still writes
  the scaffold map (so chromosome assignment, hap label, and
  orientation per contig are recorded), but **skips the per-feature-set
  ``rewrite_bed`` step** that materialises full-resolution scaffolded
  BEDs on disk. The scaffold map is instead applied at *bin time*:
  the binner runs against the smoothed (unscaffolded) BED, and a
  post-bin :func:`rewrite_bed` pass on the small binned output
  performs the rename + coordinate mirroring. The final binned
  scaffolded BEDs and karyotype SVGs are equivalent (modulo sub-bin
  alignment at contig length boundaries); the multi-minute smoothed
  rewrite is replaced by a microsecond-scale binned rewrite. Saves
  ~5-10 min on whole-genome HG002 for the default karyotype cascade
  (all 6 feature sets x 3 modes).
- ``scaffold_run`` and ``centromeres_run`` both accept a new
  ``write_scaffolded_beds`` parameter (default ``True``) that
  controls whether full-resolution scaffolded BEDs are materialised.
  ``karyotype_run`` forwards ``write_scaffolded_beds=scaffolding`` to
  both calls. ``_ensure_binned_scaffolded`` in both
  :mod:`karyoscope.core.karyotype_run` and
  :mod:`karyoscope.core.centromeres` learned a fallback path: when
  the scaffolded BED is missing, bin the smoothed BED into a temp
  file and apply the map via :func:`rewrite_bed` to produce the
  binned scaffolded BED directly. The downstream renderer is
  unchanged -- the on-disk format of the binned scaffolded BED is
  identical to the historical path.
- ``karyoscope karyotype`` now exposes ``--bgzip/--no-bgzip`` (default
  on), matching the surface of ``annotate``, ``scaffold``, and
  ``centromeres``. The flag controls compression of the intermediate
  scaffolded BEDs and the centromeres BED produced by the cascade;
  the SVG / PDF / PNG outputs are unaffected. Previously karyotype
  always bgzipped its intermediates with no override, which was
  inconsistent and inconvenient for benchmarking.
- ``scaffold_run._ensure_annotated`` no longer hardcodes ``bgzip=True``
  when auto-deriving missing annotation BEDs; it now inherits the
  scaffold-level setting. This makes ``karyoscope scaffold --no-bgzip``
  on a fresh input produce uncompressed annotation BEDs end-to-end
  (matching what a manual ``karyoscope annotate --no-bgzip`` would
  have produced). Downstream readers already handle ``.bed`` and
  ``.bed.gz`` transparently via ``chunked_seq_reader``, so this only
  affects on-disk storage, not behaviour.
- ``annotate`` now parallelises the per-(feature_set, kind) temp-file
  concatenation phase via a :class:`ThreadPoolExecutor` of
  ``--threads`` workers. Each (fs, kind) target reads a distinct
  per-FS temp directory and writes a distinct output file, so there's
  no cross-task contention -- pure I/O work that benefits from
  /scratch SSD parallel throughput. On HG002 the concat phase was a
  ~3 min single-threaded tail at the end of the smoothing pass; with
  ``-t 16`` it should drop to closer to 30 s. Threads (not processes)
  are right here: ``shutil.copyfileobj`` releases the GIL across the
  blocking read/write so we get real parallelism. New
  ``concat pass: N temp-file group(s) (threads=M)`` /
  ``concat pass complete in Xs`` log lines surface the timing.
- ``_bgzip_file`` now forwards the caller's ``--threads`` to bgzip as
  ``-@ N``, parallelising single-file compression. Wired through every
  bgzip caller: ``annotate`` (per-feature-set presmoothed / smoothed
  BEDs), ``scaffold`` (scaffolded BEDs and the scaffolded FASTA), and
  ``centromeres`` (the centromeres BED). On HG002 the bgzip pass was
  ~1-2 min of single-threaded compression for the 12 per-feature-set
  BEDs; ``-t 16`` should drop that to ~10-20 s. Per-file log lines now
  surface ``threads=N`` for visibility. When ``threads=1`` the ``-@``
  flag is omitted entirely for cleanest subprocess invocation.
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

## [1.1.0] - YYYY-MM-DD
-->

[Unreleased]: https://github.com/barthel-lab/KaryoScope/compare/v2.2.0...HEAD
[2.2.0]: https://github.com/barthel-lab/KaryoScope/releases/tag/v2.2.0
[2.1.0]: https://github.com/barthel-lab/KaryoScope/releases/tag/v2.1.0
[2.0.0]: https://github.com/barthel-lab/KaryoScope/releases/tag/v2.0.0
[1.1.0]: https://github.com/barthel-lab/KaryoScope/releases/tag/v1.1.0
[1.0.0]: https://github.com/barthel-lab/KaryoScope/releases/tag/v1.0.0
