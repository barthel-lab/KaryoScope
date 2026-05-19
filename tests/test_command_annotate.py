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
    """`annotate` with all defaults produces bgzipped BEDs for every feature set.

    Because the rebuilt dummy db's KMC counters are exactly 1, 2, 3
    (one per k-mer in the seed sequence), and ``query_fasta``'s
    seq_with_features is the seed verbatim, we can assert the exact
    output records here — not just structural validity.
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

    chrom_path = outdir / "my_query.KS_dummy_test_v1.chromosome.presmoothed.bed.gz"
    region_path = outdir / "my_query.KS_dummy_test_v1.region.presmoothed.bed.gz"
    assert chrom_path.is_file()
    assert region_path.is_file()

    # The combined intermediate should be deleted by default.
    combined = outdir / "my_query.KS_dummy_test_v1.combined.presmoothed.featureIDs.bed"
    assert not combined.exists(), "intermediate should be deleted without --keep-intermediates"

    # seq_with_features is the 23 bp seed → 3 k-mers with featureIDs 1, 2, 3.
    # In the chromosome BED: featureIDs 1 and 2 both translate to 'chr1' and
    # should merge across positions 0..2; featureID 3 translates to 'chr2'
    # at position 2..3. seq_novel is a long homopolymer → all k-mers featureID 0
    # → a single 'novel' run spanning the full sequence.
    assert _read_bed(chrom_path) == [
        ("seq_with_features", 0, 2, "chr1"),
        ("seq_with_features", 2, 3, "chr2"),
        ("seq_novel", 0, 13, "novel"),
    ]

    # In the region BED, each of the 3 features has a distinct name (rA,
    # rB, rC), so nothing merges.
    assert _read_bed(region_path) == [
        ("seq_with_features", 0, 1, "rA"),
        ("seq_with_features", 1, 2, "rB"),
        ("seq_with_features", 2, 3, "rC"),
        ("seq_novel", 0, 13, "novel"),
    ]


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
    """--feature-set restricts which BEDs are produced."""
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
    assert not (outdir / "my_query.KS_dummy_test_v1.region.presmoothed.bed").exists()


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
