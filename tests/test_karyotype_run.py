"""Unit tests for :mod:`karyoscope.core.karyotype_run` helpers.

These cover the pure-Python staleness guard for binned-scaffolded BEDs
(no cairo / rendering needed). The end-to-end rendering tests live in
``test_karyotype.py``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import karyoscope.core.karyotype_run as kr
from karyoscope.core.io.scaffold_map import MapRow, map_signature
from karyoscope.core.karyotype_run import (
    _binned_bed_is_current,
    _binned_combined_scaffolded_bed_path,
    _binned_mapsig_path,
    _binned_scaffolded_bed_path,
    _combined_centromeres_bed_path,
    _combined_scaffolded_bed_path,
    _common_base,
    _ensure_binned_scaffolded,
    _write_binned_mapsig,
)
from karyoscope.exceptions import KaryotypeError


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


# --- _ensure_binned_scaffolded: the reuse/rebuild/fail decisions ------
#
# The mapsig primitives above say whether a binned BED is current; this
# function decides what to DO about it. These tests pin those decisions
# with the binning layer stubbed, so each branch is observable.


class TestEnsureBinnedScaffolded:
    @property
    def ROWS(self) -> list[MapRow]:
        return [_row("chr1_hap1_a", hap="hap1"), _row("chr1_hap2_b", hap="hap2")]

    def _call(self, out_dir: Path, **overrides):
        kwargs = dict(
            out_dir=out_dir,
            stem="s",
            db_id="db",
            fs="region",
            bin_size=100_000,
            leaf_set=set(),
            auto=True,
            input_name="s.fa",
            threads=1,
            map_rows=self.ROWS,
        )
        kwargs.update(overrides)
        return _ensure_binned_scaffolded(**kwargs)

    def _out_path(self, out_dir: Path) -> Path:
        return _binned_scaffolded_bed_path(out_dir, "s", "db", "region", 100_000)

    def test_reuses_a_current_binned_bed(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        out = self._out_path(tmp_path)
        out.write_bytes(b"")
        _write_binned_mapsig(out, self.ROWS)
        monkeypatch.setattr(
            kr, "bin_features", lambda *a, **kw: pytest.fail("must not rebin a current BED")
        )
        assert self._call(tmp_path) == out

    def test_stale_map_without_auto_fails_loudly(self, tmp_path: Path) -> None:
        out = self._out_path(tmp_path)
        out.write_bytes(b"")
        _write_binned_mapsig(
            out, [_row("chr1_hap1_a", hap="hap1"), _row("chr1_hap1_b", hap="hap1")]
        )
        with pytest.raises(KaryotypeError, match="stale binned scaffolded BED"):
            self._call(tmp_path, auto=False)

    def test_missing_without_auto_fails_loudly(self, tmp_path: Path) -> None:
        with pytest.raises(KaryotypeError, match="missing binned scaffolded BED"):
            self._call(tmp_path, auto=False)

    def test_stale_map_rebuilds_and_rerecords_the_signature(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        out = self._out_path(tmp_path)
        out.write_bytes(b"")
        _write_binned_mapsig(
            out, [_row("chr1_hap1_a", hap="hap1"), _row("chr1_hap1_b", hap="hap1")]
        )
        src = tmp_path / "s.db.region.smoothed.scaffolded.bed.gz"
        src.write_bytes(b"")
        binned_from: list[Path] = []

        def fake_bin(source, dest, **kw):
            binned_from.append(source)
            dest.write_bytes(b"rebuilt")

        monkeypatch.setattr(kr, "bin_features", fake_bin)
        result = self._call(tmp_path)
        assert result == out
        assert binned_from == [src]
        assert _binned_bed_is_current(out, self.ROWS)

    def test_no_scaffolded_bed_falls_back_to_annotation_plus_map(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        annotation = tmp_path / "s.db.region.smoothed.bed.gz"
        annotation.write_bytes(b"")
        binned_from: list[Path] = []
        rewrites: list[tuple[Path, Path]] = []

        def fake_bin(source, dest, **kw):
            binned_from.append(source)
            dest.write_bytes(b"binned")

        def fake_rewrite(source, dest, *, map_rows):
            assert map_rows == self.ROWS
            rewrites.append((source, dest))
            dest.write_bytes(b"renamed")

        monkeypatch.setattr(kr, "bin_features", fake_bin)
        monkeypatch.setattr(kr, "rewrite_bed", fake_rewrite)
        out = self._call(tmp_path)
        assert binned_from == [annotation]
        assert rewrites and rewrites[0][1] == out
        assert _binned_bed_is_current(out, self.ROWS)

    def test_no_scaffolded_bed_and_no_map_is_an_error(self, tmp_path: Path) -> None:
        with pytest.raises(KaryotypeError, match="no scaffold map provided"):
            self._call(tmp_path, map_rows=None)

    def test_nothing_to_bin_from_is_an_error(self, tmp_path: Path) -> None:
        with pytest.raises(KaryotypeError, match="smoothed BED missing"):
            self._call(tmp_path)
