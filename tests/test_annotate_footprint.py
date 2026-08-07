"""Tests for annotate's output-size estimate and dependency preflight.

The estimate exists so a user learns their disk is too small in a second
rather than twenty minutes. It only has to be roughly right, but "roughly"
is pinned here against a real measured run so a future change to the
constants can't silently drift by an order of magnitude.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.core.annotate import (
    GZIP_EXPANSION_FACTOR,
    _annotate_dependencies,
    estimate_input_bases,
    estimate_output_bytes,
)

#: The reference measurement, from a full HG002 v1.1 run against
#: HKS_human_CHM13_v2 at 16 threads (2026-07, --no-bgzip, six feature sets):
#: 5,999,424,718 bases in; 21.70 GB of presmoothed BED plus 7.26 GB of
#: smoothed BED out, with a largest-single-set lookup TSV of 5.42 GB live
#: at peak.
_HG002_BASES = 5_999_424_718
_HG002_MEASURED_PEAK = 34_400_000_000

#: The same run's peak once ``hks lookup`` began writing the presmoothed BED
#: itself: the 5.42 GB TSV was a second copy of a file we were keeping
#: anyway, and simply stopped existing. Outputs alone, 21.70 + 7.26 GB.
_HG002_HKS_PEAK = 28_960_000_000


# --- input size -------------------------------------------------------


def test_bases_come_from_a_fai_when_one_exists(tmp_path: Path) -> None:
    """An exact count beats any heuristic, and assemblies usually ship one."""
    fasta = tmp_path / "asm.fasta.gz"
    fasta.write_bytes(b"not actually gzip")
    (tmp_path / "asm.fasta.gz.fai").write_text(
        "chr1\t248956422\t112\t60\t61\nchr2\t242193529\t9000\t60\t61\n"
    )
    assert estimate_input_bases(fasta) == 248956422 + 242193529


def test_a_malformed_fai_falls_back_rather_than_raising(tmp_path: Path) -> None:
    fasta = tmp_path / "asm.fasta"
    fasta.write_bytes(b"A" * 1000)
    (tmp_path / "asm.fasta.fai").write_text("this is not a fai\n")
    assert estimate_input_bases(fasta) == pytest.approx(980, rel=0.01)


def test_plain_fasta_size_is_close_to_its_base_count(tmp_path: Path) -> None:
    fasta = tmp_path / "asm.fasta"
    fasta.write_bytes(b"A" * 1_000_000)
    assert estimate_input_bases(fasta) == pytest.approx(1_000_000, rel=0.05)


def test_fastq_is_discounted_for_quality_strings(tmp_path: Path) -> None:
    """A FASTQ carries one quality character per base, so ~half its bytes."""
    fq = tmp_path / "reads.fastq"
    fq.write_bytes(b"A" * 1_000_000)
    assert estimate_input_bases(fq) < 600_000


def test_gzipped_input_is_expanded_before_scaling(tmp_path: Path) -> None:
    plain = tmp_path / "asm.fasta"
    gz = tmp_path / "other.fasta.gz"
    plain.write_bytes(b"A" * 1_000_000)
    gz.write_bytes(b"A" * 1_000_000)
    assert estimate_input_bases(gz) == pytest.approx(
        estimate_input_bases(plain) * GZIP_EXPANSION_FACTOR, rel=0.01
    )


def test_missing_input_estimates_zero_rather_than_raising(tmp_path: Path) -> None:
    assert estimate_input_bases(tmp_path / "gone.fasta") == 0


def test_gz_heuristic_lands_near_the_real_hg002_file() -> None:
    """1.77 GB of gzip -> 6.0 Gbp. The constant must stay in that ballpark."""
    hg002_gz_bytes = 1_769_801_343
    estimated = hg002_gz_bytes * GZIP_EXPANSION_FACTOR * 0.98
    assert estimated == pytest.approx(_HG002_BASES, rel=0.05)


# --- output size ------------------------------------------------------


def test_estimate_matches_the_measured_hg002_run() -> None:
    estimate = estimate_output_bytes(
        input_bases=_HG002_BASES,
        n_feature_sets=6,
        keep_presmoothed=True,
        smooth=True,
    )
    assert estimate == pytest.approx(_HG002_MEASURED_PEAK, rel=0.10)


def test_estimate_scales_with_feature_set_count() -> None:
    """--feature-set is one of the escapes we suggest, so it must matter."""
    one = estimate_output_bytes(
        input_bases=_HG002_BASES, n_feature_sets=1, keep_presmoothed=True, smooth=True
    )
    six = estimate_output_bytes(
        input_bases=_HG002_BASES, n_feature_sets=6, keep_presmoothed=True, smooth=True
    )
    assert six > one * 3


def test_dropping_an_output_lowers_the_estimate() -> None:
    both = estimate_output_bytes(
        input_bases=_HG002_BASES, n_feature_sets=6, keep_presmoothed=True, smooth=True
    )
    smoothed_only = estimate_output_bytes(
        input_bases=_HG002_BASES, n_feature_sets=6, keep_presmoothed=False, smooth=True
    )
    presmoothed_only = estimate_output_bytes(
        input_bases=_HG002_BASES, n_feature_sets=6, keep_presmoothed=True, smooth=False
    )
    assert smoothed_only < both
    assert presmoothed_only < both
    # Presmoothed BEDs are the bulk of the output (0.60 vs 0.20 bytes/base).
    assert presmoothed_only > smoothed_only


def test_transient_intermediate_is_counted_when_there_is_one() -> None:
    """A single feature set still needs room for its intermediate alongside."""
    single = estimate_output_bytes(
        input_bases=1_000_000_000, n_feature_sets=1, keep_presmoothed=False, smooth=True
    )
    outputs_only = 1_000_000_000 * 0.20
    assert single > outputs_only * 2


def test_hks_keeping_presmoothed_is_charged_for_no_intermediate() -> None:
    """`hks lookup` writes the presmoothed BED itself, so nothing is duplicated.

    Charging for an intermediate that no longer exists would refuse runs
    that fit -- about 5 GB's worth on a diploid human assembly.
    """
    estimate = estimate_output_bytes(
        input_bases=_HG002_BASES,
        n_feature_sets=6,
        keep_presmoothed=True,
        smooth=True,
        index_type="hks",
    )
    assert estimate == pytest.approx(_HG002_HKS_PEAK, rel=0.10)
    assert estimate < _HG002_MEASURED_PEAK


def test_hks_discarding_presmoothed_still_needs_the_temp_file() -> None:
    """Without --keep-presmoothed the lookup output is a temp file smooth reads."""
    kept = estimate_output_bytes(
        input_bases=1_000_000_000,
        n_feature_sets=1,
        keep_presmoothed=False,
        smooth=True,
        index_type="hks",
    )
    outputs_only = 1_000_000_000 * 0.20
    assert kept > outputs_only * 2


def test_an_unknown_backend_errs_high() -> None:
    """The default must be the conservative one: over-estimating only warns."""
    assert estimate_output_bytes(
        input_bases=_HG002_BASES, n_feature_sets=6, keep_presmoothed=True, smooth=True
    ) == estimate_output_bytes(
        input_bases=_HG002_BASES,
        n_feature_sets=6,
        keep_presmoothed=True,
        smooth=True,
        index_type="kmc",
    )


# --- dependency selection ---------------------------------------------


def test_hks_database_asks_for_hks_not_kmc() -> None:
    deps = _annotate_dependencies(index_type="hks", input_path=Path("a.fasta"), bgzip=False)
    assert deps == ["hks"]


def test_kmc_database_asks_for_get_featureids() -> None:
    deps = _annotate_dependencies(index_type="kmc", input_path=Path("a.fasta"), bgzip=False)
    assert deps == ["get_featureIDs"]


def test_bam_input_adds_samtools_and_bgzip_adds_bgzip() -> None:
    deps = _annotate_dependencies(index_type="hks", input_path=Path("a.bam"), bgzip=True)
    assert set(deps) == {"hks", "samtools", "bgzip"}


def test_fasta_input_does_not_require_samtools() -> None:
    """Users shouldn't be told to install a tool their run will never call."""
    deps = _annotate_dependencies(index_type="hks", input_path=Path("a.fasta.gz"), bgzip=True)
    assert "samtools" not in deps
