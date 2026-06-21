"""Unit tests for :mod:`karyoscope.core.karyotype_run` helpers.

These cover the pure-Python staleness guard for binned-scaffolded BEDs
(no cairo / rendering needed). The end-to-end rendering tests live in
``test_karyotype.py``.
"""

from __future__ import annotations

from pathlib import Path

from karyoscope.core.io.scaffold_map import MapRow, map_signature
from karyoscope.core.karyotype_run import (
    _binned_bed_is_current,
    _binned_combined_scaffolded_bed_path,
    _binned_mapsig_path,
    _combined_centromeres_bed_path,
    _combined_scaffolded_bed_path,
    _common_base,
    _write_binned_mapsig,
)


class TestCommonBase:
    def test_single_input_returns_its_stem(self) -> None:
        assert _common_base(["GM00392.assembly"]) == "GM00392.assembly"

    def test_haplotype_pair_collapses_to_sample(self) -> None:
        assert _common_base(["GM04890.haplotype1", "GM04890.haplotype2"]) == "GM04890"

    def test_hap_pair_collapses_to_sample(self) -> None:
        assert _common_base(["BJ.hap1", "BJ.hap2"]) == "BJ"

    def test_maternal_paternal_collapses_to_sample(self) -> None:
        assert _common_base(["HG002.maternal", "HG002.paternal"]) == "HG002"

    def test_three_inputs_with_unassigned(self) -> None:
        assert _common_base(["S.hap1", "S.hap2", "S.unassigned"]) == "S"

    def test_deeper_shared_prefix_kept(self) -> None:
        assert _common_base(["HG002.hifiasm.hap1", "HG002.hifiasm.hap2"]) == "HG002.hifiasm"

    def test_underscore_separator(self) -> None:
        assert _common_base(["GM04890_hap1", "GM04890_hap2"]) == "GM04890"

    def test_no_common_separator_prefix_falls_back_to_first(self) -> None:
        # Shared chars but no separator boundary -> don't emit a partial
        # token; fall back to the first stem.
        assert _common_base(["abc1", "abc2"]) == "abc1"

    def test_no_overlap_falls_back_to_first(self) -> None:
        assert _common_base(["alpha", "beta"]) == "alpha"

    def test_identical_stems_kept_whole(self) -> None:
        assert _common_base(["sampleX", "sampleX"]) == "sampleX"

    def test_empty_returns_empty(self) -> None:
        assert _common_base([]) == ""


def _row(new_name: str, *, hap: str, flipped: bool = False) -> MapRow:
    return MapRow(
        new_name=new_name,
        original_name=new_name.split("_")[-1],
        input_file="x.fa",
        hap=hap,
        chromosome=new_name.split("_")[0],
        flipped=flipped,
        length=1000,
        stats="PCQ",
    )


def _binned(tmp_path: Path) -> Path:
    # The helpers only touch the sidecar, not the binned BED's bytes, so
    # an empty placeholder file is enough to stand in for a built BED.
    p = tmp_path / "s.db.region.smoothed.scaffolded.binned100000.bed.gz"
    p.write_bytes(b"")
    return p


def test_mapsig_path_is_sidecar(tmp_path: Path) -> None:
    binned = _binned(tmp_path)
    assert _binned_mapsig_path(binned) == binned.with_name(binned.name + ".mapsig")


def test_fresh_sidecar_is_current(tmp_path: Path) -> None:
    binned = _binned(tmp_path)
    rows = [_row("chr1_hap1_a", hap="hap1"), _row("chr1_hap2_b", hap="hap2")]
    _write_binned_mapsig(binned, rows)
    assert _binned_bed_is_current(binned, rows)


def test_changed_map_is_not_current(tmp_path: Path) -> None:
    binned = _binned(tmp_path)
    # Built when both contigs were hap1 (the pre-fix collapse) ...
    _write_binned_mapsig(binned, [_row("chr1_hap1_a", hap="hap1"), _row("chr1_hap1_b", hap="hap1")])
    # ... now hap inference is corrected: the second contig is hap2.
    corrected = [_row("chr1_hap1_a", hap="hap1"), _row("chr1_hap2_b", hap="hap2")]
    assert not _binned_bed_is_current(binned, corrected)


def test_missing_sidecar_is_not_current(tmp_path: Path) -> None:
    # A binned BED from before the guard existed has no sidecar; it must
    # be treated as stale (rebuilt once) rather than trusted.
    binned = _binned(tmp_path)
    assert not _binned_bed_is_current(binned, [_row("chr1_hap1_a", hap="hap1")])


def test_empty_sidecar_is_not_current(tmp_path: Path) -> None:
    binned = _binned(tmp_path)
    _binned_mapsig_path(binned).write_text("")
    assert not _binned_bed_is_current(binned, [_row("chr1_hap1_a", hap="hap1")])


def test_none_map_rows_preserves_legacy_reuse(tmp_path: Path) -> None:
    # When the caller can't supply a map to check against, fall back to
    # the historical reuse-if-present behaviour.
    binned = _binned(tmp_path)
    assert _binned_bed_is_current(binned, None)
    _write_binned_mapsig(binned, None)  # no-op, must not raise or create a file
    assert not _binned_mapsig_path(binned).exists()


def test_sidecar_contents_match_signature(tmp_path: Path) -> None:
    binned = _binned(tmp_path)
    rows = [_row("chr1_hap1_a", hap="hap1")]
    _write_binned_mapsig(binned, rows)
    assert _binned_mapsig_path(binned).read_text().strip() == map_signature(rows)


# --- combine-chromosome path naming ---------------------------------


def test_combined_scaffolded_bed_path_defaults_to_gz(tmp_path: Path) -> None:
    p = _combined_scaffolded_bed_path(tmp_path, "samp", "DB", "region")
    assert p.name == "samp.DB.region.smoothed.scaffolded.combined_chromosomes.bed.gz"


def test_combined_scaffolded_bed_path_prefers_existing_plain(tmp_path: Path) -> None:
    plain = tmp_path / "samp.DB.region.smoothed.scaffolded.combined_chromosomes.bed"
    plain.write_text("")
    assert _combined_scaffolded_bed_path(tmp_path, "samp", "DB", "region") == plain


def test_binned_combined_scaffolded_bed_path(tmp_path: Path) -> None:
    p = _binned_combined_scaffolded_bed_path(tmp_path, "samp", "DB", "region", 1_000_000)
    assert p.name == (
        "samp.DB.region.smoothed.scaffolded.combined_chromosomes.binned1000000.bed.gz"
    )


def test_combined_centromeres_bed_path(tmp_path: Path) -> None:
    p = _combined_centromeres_bed_path(tmp_path, "samp", "DB")
    assert p.name == "samp.DB.centromeres.combined_chromosomes.bed.gz"
