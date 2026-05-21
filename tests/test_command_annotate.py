"""Integration tests for ``karyoscope annotate``.

Marked ``@pytest.mark.integration`` because they need the ``get_featureIDs``
C++ binary (and ``bgzip`` for the default-codepath tests).
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest
from click.testing import CliRunner

from karyoscope.cli import main
from karyoscope.core.annotate import _derive_input_basename, _split_combined_bed
from karyoscope.core.io.features import parse_features

pytestmark = pytest.mark.integration


# --- helpers ----------------------------------------------------------


def _has_bgzip() -> bool:
    """True if bgzip is on PATH (skip-aware for the default test)."""
    import shutil

    return shutil.which("bgzip") is not None


@pytest.fixture
def query_fasta(tmp_path: Path) -> Path:
    """A FASTA whose k-mers exercise both real-feature and 'novel' paths.

    The dummy db's index contains exactly three 21-mers (m1, m2, m3)
    drawn from the seed ``ACGTGCTAGCTAGGCTATCGTAC``, with counts 1, 2, 3
    → featureIDs 1, 2, 3. features.tsv maps these to:

    * featureID 1 → chr1 / rA
    * featureID 2 → chr1 / rB
    * featureID 3 → chr2 / rC

    ``seq_with_features`` is the seed verbatim, so its 3 k-mers will
    resolve to featureIDs 1, 2, 3 in order. ``seq_novel`` is a long
    homopolymer whose 21-mer is guaranteed not to be in the index, so
    all its k-mers will resolve to the 'novel' sentinel.
    """
    fa = tmp_path / "my_query.fa"
    fa.write_text(
        ">seq_with_features\n"
        "ACGTGCTAGCTAGGCTATCGTAC\n"
        ">seq_novel\n"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
    )
    return fa


def _read_bed(path: Path) -> list[tuple[str, int, int, str]]:
    """Read a BED (plain or .bed.gz) and return parsed records."""
    text = gzip.decompress(path.read_bytes()).decode() if path.suffix == ".gz" else path.read_text()
    out: list[tuple[str, int, int, str]] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        seq, start, end, name = line.split("\t")
        out.append((seq, int(start), int(end), name))
    return out


# --- _derive_input_basename (pure function, no binary needed) ---------


def test_derive_basename_strips_known_extensions() -> None:
    assert _derive_input_basename(Path("foo.fa.gz")) == "foo"
    assert _derive_input_basename(Path("path/to/bar.fasta")) == "bar"
    assert _derive_input_basename(Path("baz.fna.gz")) == "baz"
    assert _derive_input_basename(Path("BAR.FA.GZ")) == "BAR"  # case-insensitive
    assert _derive_input_basename(Path("name_only")) == "name_only"
    # Falls back to .stem for unknown extensions
    assert _derive_input_basename(Path("foo.unknown")) == "foo"


# --- _split_combined_bed (unit test with hand-crafted data) -----------


def test_split_combined_bed_translates_and_merges(tmp_path: Path, unpacked_dummy_db: Path) -> None:
    """The splitter should translate ids and merge adjacent same-name runs."""
    # Hand-crafted combined BED using REAL feature ids from the dummy db
    # (features.tsv has rows 1->chr1/rA, 2->chr1/rB, 3->chr2/rC).
    combined = tmp_path / "combined.bed"
    combined.write_text(
        # Adjacent rows mapping to different region ids but same chromosome
        # should merge in the chromosome BED but not in the region BED.
        "seqA\t0\t10\t1\n"  # chr1, rA
        "seqA\t10\t20\t2\n"  # chr1, rB
        "seqA\t20\t30\t3\n"  # chr2, rC
        # A 'novel' run (feature_id 0):
        "seqA\t30\t40\t0\n"  # novel
        # New sequence: starts a new run regardless of label
        "seqB\t0\t15\t1\n"  # chr1, rA
    )
    features = parse_features(unpacked_dummy_db / "features.tsv")

    out_paths = {
        "chromosome": tmp_path / "chrom.bed",
        "region": tmp_path / "region.bed",
    }
    _split_combined_bed(combined, ["chromosome", "region"], features, out_paths)

    chrom = _read_bed(out_paths["chromosome"])
    region = _read_bed(out_paths["region"])

    # Chromosome BED: rows 0-20 both map to chr1 and merge; 20-30 is chr2;
    # 30-40 is novel; seqB starts a fresh run at chr1.
    assert chrom == [
        ("seqA", 0, 20, "chr1"),
        ("seqA", 20, 30, "chr2"),
        ("seqA", 30, 40, "novel"),
        ("seqB", 0, 15, "chr1"),
    ]

    # Region BED: each row stays distinct (different region per row).
    assert region == [
        ("seqA", 0, 10, "rA"),
        ("seqA", 10, 20, "rB"),
        ("seqA", 20, 30, "rC"),
        ("seqA", 30, 40, "novel"),
        ("seqB", 0, 15, "rA"),
    ]


def test_split_combined_bed_raises_on_unknown_feature_id(
    tmp_path: Path, unpacked_dummy_db: Path
) -> None:
    """A featureID present in the BED but absent from features.tsv is an error.

    This usually signals a mismatch between the KMC index and
    features.tsv (different builds, or a stale features file). We
    refuse to silently emit a placeholder string because 'Unknown' can
    be a legitimate feature name in real databases (e.g., the repeats
    set).
    """
    from karyoscope.core.io.features import FeaturesError

    combined = tmp_path / "combined.bed"
    # featureID 999 has no row in the dummy db's features.tsv (which
    # only has rows for ids 1, 2, 3).
    combined.write_text("seqA\t0\t10\t999\n")
    features = parse_features(unpacked_dummy_db / "features.tsv")

    out_paths = {"chromosome": tmp_path / "chrom.bed"}
    with pytest.raises(FeaturesError, match="feature id 999 is not in features"):
        _split_combined_bed(combined, ["chromosome"], features, out_paths)


# --- end-to-end CLI ---------------------------------------------------


@pytest.fixture
def cli_with_populated_db(
    cli_runner: CliRunner,
    populated_db_root: Path,
    query_fasta: Path,
) -> tuple[CliRunner, Path, Path]:
    """A populated db root and query FASTA; returns (runner, db_root, fasta)."""
    return cli_runner, populated_db_root, query_fasta


def test_annotate_end_to_end_default_flags(
    cli_with_populated_db: tuple[CliRunner, Path, Path],
    tmp_path: Path,
) -> None:
    """`annotate` with all defaults produces BOTH presmoothed and smoothed BEDs.

    Default behaviour is ``--smooth --keep-presmoothed --bgzip``: write
    bgzipped versions of both outputs for every feature set.
    """
    if not _has_bgzip():
        pytest.skip("bgzip not on PATH; install htslib")

    runner, db_root, query_fasta = cli_with_populated_db
    outdir = tmp_path / "out"

    result = runner.invoke(
        main,
        [
            "annotate",
            "-i",
            str(query_fasta),
            "-o",
            str(outdir),
            "--db-root",
            str(db_root),
            "-t",
            "1",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    # All four files should exist (2 feature sets x 2 output types)
    for variant in ("presmoothed", "smoothed"):
        for fs in ("chromosome", "region"):
            p = outdir / f"my_query.KS_dummy_test_v1.{fs}.{variant}.bed.gz"
            assert p.is_file(), f"missing {p}"

    # The combined intermediate should be deleted by default.
    combined = outdir / "my_query.KS_dummy_test_v1.combined.presmoothed.featureIDs.bed"
    assert not combined.exists(), "intermediate should be deleted without --keep-intermediates"

    # Presmoothed content unchanged from previous stages.
    chrom_pre = outdir / "my_query.KS_dummy_test_v1.chromosome.presmoothed.bed.gz"
    assert _read_bed(chrom_pre) == [
        ("seq_with_features", 0, 2, "chr1"),
        ("seq_with_features", 2, 3, "chr2"),
        ("seq_novel", 0, 13, "novel"),
    ]
    region_pre = outdir / "my_query.KS_dummy_test_v1.region.presmoothed.bed.gz"
    assert _read_bed(region_pre) == [
        ("seq_with_features", 0, 1, "rA"),
        ("seq_with_features", 1, 2, "rB"),
        ("seq_with_features", 2, 3, "rC"),
        ("seq_novel", 0, 13, "novel"),
    ]

    # Smoothed content: for this query the dummy db has no novel runs
    # between known features (seq_with_features has no gaps; seq_novel
    # is entirely novel with no flanking features), so smoothing should
    # produce essentially the same output as presmoothed. The key thing
    # to verify is that the file exists and is well-formed.
    region_smo = outdir / "my_query.KS_dummy_test_v1.region.smoothed.bed.gz"
    smoothed_records = _read_bed(region_smo)
    assert smoothed_records, "smoothed BED should not be empty"
    # Every record should be a valid 4-column row (already enforced by _read_bed).
    # Each sequence should appear at least once.
    seqs = {r[0] for r in smoothed_records}
    assert seqs == {"seq_with_features", "seq_novel"}


def test_annotate_with_no_bgzip(
    cli_with_populated_db: tuple[CliRunner, Path, Path],
    tmp_path: Path,
) -> None:
    """--no-bgzip writes plain .bed files."""
    runner, db_root, query_fasta = cli_with_populated_db
    outdir = tmp_path / "out"

    result = runner.invoke(
        main,
        [
            "annotate",
            "-i",
            str(query_fasta),
            "-o",
            str(outdir),
            "--db-root",
            str(db_root),
            "-t",
            "1",
            "--no-bgzip",
            "--no-smooth",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    chrom_path = outdir / "my_query.KS_dummy_test_v1.chromosome.presmoothed.bed"
    region_path = outdir / "my_query.KS_dummy_test_v1.region.presmoothed.bed"
    assert chrom_path.is_file()
    assert region_path.is_file()
    # And NOT the .gz versions:
    assert not (outdir / f"{chrom_path.name}.gz").exists()


def test_annotate_with_keep_intermediates(
    cli_with_populated_db: tuple[CliRunner, Path, Path],
    tmp_path: Path,
) -> None:
    """--keep-intermediates retains the combined BED from the C++ step."""
    runner, db_root, query_fasta = cli_with_populated_db
    outdir = tmp_path / "out"

    result = runner.invoke(
        main,
        [
            "annotate",
            "-i",
            str(query_fasta),
            "-o",
            str(outdir),
            "--db-root",
            str(db_root),
            "-t",
            "1",
            "--no-bgzip",
            "--no-smooth",
            "--keep-intermediates",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    combined = outdir / "my_query.KS_dummy_test_v1.combined.presmoothed.featureIDs.bed"
    assert combined.is_file()


def test_annotate_filters_by_feature_set(
    cli_with_populated_db: tuple[CliRunner, Path, Path],
    tmp_path: Path,
) -> None:
    """--feature-set restricts which BEDs are produced (both variants)."""
    runner, db_root, query_fasta = cli_with_populated_db
    outdir = tmp_path / "out"

    result = runner.invoke(
        main,
        [
            "annotate",
            "-i",
            str(query_fasta),
            "-o",
            str(outdir),
            "--db-root",
            str(db_root),
            "-t",
            "1",
            "--no-bgzip",
            "--feature-set",
            "chromosome",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    assert (outdir / "my_query.KS_dummy_test_v1.chromosome.presmoothed.bed").is_file()
    assert (outdir / "my_query.KS_dummy_test_v1.chromosome.smoothed.bed").is_file()
    # The 'region' set was not requested.
    assert not (outdir / "my_query.KS_dummy_test_v1.region.presmoothed.bed").exists()
    assert not (outdir / "my_query.KS_dummy_test_v1.region.smoothed.bed").exists()


def test_annotate_unknown_feature_set_fails(
    cli_with_populated_db: tuple[CliRunner, Path, Path],
    tmp_path: Path,
) -> None:
    """Asking for a feature set that's not in the manifest is a clean error."""
    runner, db_root, query_fasta = cli_with_populated_db
    outdir = tmp_path / "out"

    result = runner.invoke(
        main,
        [
            "annotate",
            "-i",
            str(query_fasta),
            "-o",
            str(outdir),
            "--db-root",
            str(db_root),
            "--feature-set",
            "no_such_set",
        ],
    )
    assert result.exit_code != 0
    assert "no_such_set" in result.output


def test_annotate_no_database_installed_fails(
    cli_runner: CliRunner,
    isolated_db_root: Path,
    query_fasta: Path,
    tmp_path: Path,
) -> None:
    """A clear error if no database is installed."""
    result = cli_runner.invoke(
        main,
        [
            "annotate",
            "-i",
            str(query_fasta),
            "-o",
            str(tmp_path / "out"),
            "--db-root",
            str(isolated_db_root),
        ],
    )
    assert result.exit_code != 0
    assert (
        "no databases installed" in result.output.lower()
        or "not installed" in result.output.lower()
    )


# --- Stage 5c: smoothing-specific flag combinations ---------------------


def test_annotate_no_smooth_skips_smoothed_output(
    cli_with_populated_db: tuple[CliRunner, Path, Path],
    tmp_path: Path,
) -> None:
    """--no-smooth produces only the presmoothed BED, no smoothed one."""
    runner, db_root, query_fasta = cli_with_populated_db
    outdir = tmp_path / "out"

    result = runner.invoke(
        main,
        [
            "annotate",
            "-i",
            str(query_fasta),
            "-o",
            str(outdir),
            "--db-root",
            str(db_root),
            "-t",
            "1",
            "--no-bgzip",
            "--no-smooth",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    for fs in ("chromosome", "region"):
        assert (outdir / f"my_query.KS_dummy_test_v1.{fs}.presmoothed.bed").is_file()
        assert not (outdir / f"my_query.KS_dummy_test_v1.{fs}.smoothed.bed").exists()


def test_annotate_no_keep_presmoothed_writes_only_smoothed(
    cli_with_populated_db: tuple[CliRunner, Path, Path],
    tmp_path: Path,
) -> None:
    """--no-keep-presmoothed produces only the smoothed BED, no presmoothed."""
    runner, db_root, query_fasta = cli_with_populated_db
    outdir = tmp_path / "out"

    result = runner.invoke(
        main,
        [
            "annotate",
            "-i",
            str(query_fasta),
            "-o",
            str(outdir),
            "--db-root",
            str(db_root),
            "-t",
            "1",
            "--no-bgzip",
            "--no-keep-presmoothed",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    for fs in ("chromosome", "region"):
        assert (outdir / f"my_query.KS_dummy_test_v1.{fs}.smoothed.bed").is_file()
        assert not (outdir / f"my_query.KS_dummy_test_v1.{fs}.presmoothed.bed").exists()


def test_annotate_no_smooth_no_keep_presmoothed_is_error(
    cli_with_populated_db: tuple[CliRunner, Path, Path],
    tmp_path: Path,
) -> None:
    """Combining --no-smooth and --no-keep-presmoothed errors cleanly.

    There'd be no output to write, so the CLI should refuse rather than
    silently produce nothing.
    """
    runner, db_root, query_fasta = cli_with_populated_db
    outdir = tmp_path / "out"

    result = runner.invoke(
        main,
        [
            "annotate",
            "-i",
            str(query_fasta),
            "-o",
            str(outdir),
            "--db-root",
            str(db_root),
            "-t",
            "1",
            "--no-smooth",
            "--no-keep-presmoothed",
        ],
    )
    assert result.exit_code != 0
    assert "no output would be produced" in result.output.lower()


def test_annotate_smoothing_promotes_novel_gap_to_lca(
    cli_with_populated_db: tuple[CliRunner, Path, Path],
    tmp_path: Path,
) -> None:
    """End-to-end smoothing test that actually exercises LCA promotion.

    Build a query FASTA that starts with the dummy db's seed (so the
    first k-mers resolve to featureID 1 → region=rA), then a stretch
    of homopolymer (no k-mers in the index → featureID 0 → novel),
    then the seed shifted by 1 base (so the next k-mer resolves to
    featureID 2 → region=rB).

    The presmoothed region BED should show: rA, novel, rB.
    The smoothed region BED should promote the novel run to aSat,
    which is the LCA of rA and rB in the dummy db's hierarchy.

    We use a short novel gap (under the smoothing max_gap of 1000)
    to ensure the smoothing window actually spans it.
    """
    runner, db_root, _ = cli_with_populated_db
    outdir = tmp_path / "out"

    # Build a custom query FASTA for this scenario
    seed = "ACGTGCTAGCTAGGCTATCGTAC"  # 23 bp
    spacer = "TTTTTTTTTTTTTTTTTTTTTTTTTT"  # 26 bp homopolymer → 6 novel 21-mers
    # Layout: [k-mer m1][spacer][k-mer m2 starting position]
    # The seed's first 21-mer is m1 (featureID 1, region=rA)
    # The seed's second 21-mer is m2 (featureID 2, region=rB)
    # We want m1 to appear, then novel k-mers, then m2 to appear.
    # m1 ends at position 21 of the seed. m2 starts at position 1 of the
    # seed. Construct: seed[0:21] + spacer + seed[1:22] → m1, novels, m2
    constructed = seed[0:21] + spacer + seed[1:22]
    query_fa = tmp_path / "smoothing_query.fa"
    query_fa.write_text(f">smoothing_test\n{constructed}\n")

    result = runner.invoke(
        main,
        [
            "annotate",
            "-i",
            str(query_fa),
            "-o",
            str(outdir),
            "--db-root",
            str(db_root),
            "-t",
            "1",
            "--no-bgzip",
            "--feature-set",
            "region",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    pre = _read_bed(outdir / "smoothing_query.KS_dummy_test_v1.region.presmoothed.bed")
    smo = _read_bed(outdir / "smoothing_query.KS_dummy_test_v1.region.smoothed.bed")

    # Presmoothed: rA at the start, novel run in the middle, rB at the end.
    pre_features = [r[3] for r in pre]
    assert "rA" in pre_features
    assert "novel" in pre_features
    assert "rB" in pre_features

    # Smoothed: the novel run should be promoted to aSat (LCA of rA, rB).
    smo_features = [r[3] for r in smo]
    assert "rA" in smo_features
    assert "rB" in smo_features
    assert "aSat" in smo_features, (
        f"expected 'aSat' (LCA of rA, rB) in smoothed output, got: {smo_features}"
    )
    # The novel run should be gone after smoothing.
    assert "novel" not in smo_features


def test_annotate_no_preserve_order_matches_preserve_order(
    cli_with_populated_db: tuple[CliRunner, Path, Path],
    tmp_path: Path,
) -> None:
    """Reads-mode (``--no-preserve-order``) produces the same per-sequence
    content as assembly-mode (``--preserve-order``) -- just possibly with
    sequences in different order in the final BED.

    Verifies the two codepaths (per-sequence temp files vs straight
    ``imap_unordered``) emit byte-equivalent output when sorted by
    sequence + position. With a single sequence in the fixture, the
    sequence order question is moot and the files should be
    byte-identical.
    """
    runner, db_root, _ = cli_with_populated_db
    seed = "ACGTGCTAGCTAGGCTATCGTAC"
    query_fa = tmp_path / "single_seq.fa"
    query_fa.write_text(f">one_seq\n{seed}\n")

    out_preserve = tmp_path / "preserve"
    out_no_preserve = tmp_path / "no_preserve"

    for outdir, flag in (
        (out_preserve, "--preserve-order"),
        (out_no_preserve, "--no-preserve-order"),
    ):
        result = runner.invoke(
            main,
            [
                "annotate",
                "-i",
                str(query_fa),
                "-o",
                str(outdir),
                "--db-root",
                str(db_root),
                "-t",
                "1",
                "--no-bgzip",
                "--feature-set",
                "region",
                flag,
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, result.output

    pre_a = (out_preserve / "single_seq.KS_dummy_test_v1.region.presmoothed.bed").read_bytes()
    pre_b = (out_no_preserve / "single_seq.KS_dummy_test_v1.region.presmoothed.bed").read_bytes()
    assert pre_a == pre_b

    smo_a = (out_preserve / "single_seq.KS_dummy_test_v1.region.smoothed.bed").read_bytes()
    smo_b = (out_no_preserve / "single_seq.KS_dummy_test_v1.region.smoothed.bed").read_bytes()
    assert smo_a == smo_b


def test_annotate_preserve_order_with_multiple_sequences_keeps_input_order(
    cli_with_populated_db: tuple[CliRunner, Path, Path],
    tmp_path: Path,
) -> None:
    """With multiple input sequences, ``--preserve-order`` writes them
    in input order (the order they appeared in the source FASTA),
    independent of worker completion order.
    """
    runner, db_root, _ = cli_with_populated_db
    seed = "ACGTGCTAGCTAGGCTATCGTAC"
    # Three sequences, deliberately named so alphabetical order
    # (z_first, m_middle, a_last) differs from input order.
    query_fa = tmp_path / "ordered.fa"
    query_fa.write_text(f">z_first\n{seed}\n>m_middle\n{seed}\n>a_last\n{seed}\n")
    outdir = tmp_path / "out"

    result = runner.invoke(
        main,
        [
            "annotate",
            "-i",
            str(query_fa),
            "-o",
            str(outdir),
            "--db-root",
            str(db_root),
            "-t",
            "2",  # 2 workers; could finish in any order
            "--no-bgzip",
            "--feature-set",
            "region",
            "--preserve-order",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    bed = (outdir / "ordered.KS_dummy_test_v1.region.smoothed.bed").read_text()
    # First column of each line is the sequence name. Extract the
    # unique sequence names in first-occurrence order from the file
    # and verify it matches input order.
    seen: list[str] = []
    for line in bed.splitlines():
        if not line:
            continue
        seq = line.split("\t", 1)[0]
        if seq not in seen:
            seen.append(seq)
    assert seen == ["z_first", "m_middle", "a_last"]


# --- FASTQ input support --------------------------------------------------


def test_annotate_accepts_fastq_input(
    cli_with_populated_db: tuple[CliRunner, Path, Path],
    tmp_path: Path,
) -> None:
    """``annotate`` accepts ``.fastq`` input and produces equivalent
    presmoothed output to the FASTA case.

    The C++ binary reads FASTQ natively; our wrapper just needs to
    recognise the extension when deriving the output prefix.
    """
    runner, db_root, _ = cli_with_populated_db
    seed = "ACGTGCTAGCTAGGCTATCGTAC"  # same 21-mers as the FASTA fixture

    # Hand-write a tiny FASTQ. Quality string length must match sequence.
    fq_path = tmp_path / "reads.fastq"
    fq_path.write_text(
        f"@seq_with_features\n{seed}\n+\n" + ("I" * len(seed)) + "\n"
        f"@seq_novel\n{'A' * 33}\n+\n" + ("I" * 33) + "\n"
    )

    outdir = tmp_path / "out"
    result = runner.invoke(
        main,
        [
            "annotate",
            "-i",
            str(fq_path),
            "-o",
            str(outdir),
            "--db-root",
            str(db_root),
            "-t",
            "1",
            "--no-bgzip",
            "--no-smooth",
            "--feature-set",
            "region",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    # Prefix should be derived as "reads" (the FASTQ extension stripped).
    expected = outdir / "reads.KS_dummy_test_v1.region.presmoothed.bed"
    assert expected.is_file(), f"missing {expected}; got: {list(outdir.iterdir())}"

    # Same k-mer content as the FASTA fixture would produce; verifies
    # FASTQ parsing actually ran through the pipeline.
    bed = _read_bed(expected)
    features = [r[3] for r in bed]
    assert "rA" in features
    assert "rB" in features
    assert "rC" in features
    assert "novel" in features


def test_annotate_accepts_fastq_gz_input(
    cli_with_populated_db: tuple[CliRunner, Path, Path],
    tmp_path: Path,
) -> None:
    """Gzipped FASTQ (``.fastq.gz``) is also accepted; prefix derivation
    strips the full ``.fastq.gz`` suffix (not just ``.gz``)."""
    runner, db_root, _ = cli_with_populated_db
    seed = "ACGTGCTAGCTAGGCTATCGTAC"

    fq_text = (
        f"@seq_with_features\n{seed}\n+\n" + ("I" * len(seed)) + "\n"
        f"@seq_novel\n{'A' * 33}\n+\n" + ("I" * 33) + "\n"
    )
    fq_gz_path = tmp_path / "myreads.fastq.gz"
    fq_gz_path.write_bytes(gzip.compress(fq_text.encode()))

    outdir = tmp_path / "out"
    result = runner.invoke(
        main,
        [
            "annotate",
            "-i",
            str(fq_gz_path),
            "-o",
            str(outdir),
            "--db-root",
            str(db_root),
            "-t",
            "1",
            "--no-bgzip",
            "--no-smooth",
            "--feature-set",
            "region",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    # Prefix is "myreads" -- full ".fastq.gz" stripped, not just ".gz".
    expected = outdir / "myreads.KS_dummy_test_v1.region.presmoothed.bed"
    assert expected.is_file()


# --- scaffold rejects read-level inputs ----------------------------------


def test_scaffold_rejects_fastq_input(
    cli_runner: CliRunner,
    populated_db_root: Path,
    tmp_path: Path,
) -> None:
    """``scaffold`` requires a FASTA assembly; passing reads should
    fail cleanly with a message pointing the user at ``annotate``.

    Lives in test_command_annotate.py because the input-format
    validation lives in scaffold_run.py alongside the annotate
    cascade entry point; the scaffold command tests are organised
    around scaffold-specific behaviour rather than input rejection.
    """
    fq_path = tmp_path / "reads.fastq"
    fq_path.write_text("@r1\nACGT\n+\nIIII\n")

    result = cli_runner.invoke(
        main,
        [
            "scaffold",
            "-i",
            str(fq_path),
            "--db-root",
            str(populated_db_root),
        ],
    )
    assert result.exit_code != 0
    # The error message should mention FASTQ and point at annotate.
    output_lower = result.output.lower()
    assert "fasta" in output_lower
    assert "annotate" in output_lower
