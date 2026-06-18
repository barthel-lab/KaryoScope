"""Tests for the standalone scaffold-map BED remap (`karyoscope remap-bed`).

Covers the compatibility guards in
:func:`karyoscope.core.scaffold.remap_bed_with_map` and the CLI wrapper.
The algorithmic coordinate rewrite itself (flip mirroring, renaming) is owned
by :func:`rewrite_bed` and exercised in ``tests/test_scaffold.py``; here we
confirm the wrapper validates correctly and still delegates the rewrite.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from karyoscope.cli import main
from karyoscope.core.io.scaffold_map import MapRow, write_map
from karyoscope.core.scaffold import remap_bed_with_map
from karyoscope.exceptions import ScaffoldError


def _row(**kw: object) -> MapRow:
    defaults: dict[str, object] = {
        "new_name": "chr1_hap1_ctg1",
        "original_name": "ctg1",
        "input_file": "hap1.fa.gz",
        "hap": "hap1",
        "chromosome": "chr1",
        "flipped": False,
        "length": 1_000_000,
        "stats": "TPCQT",
    }
    defaults.update(kw)
    return MapRow(**defaults)  # type: ignore[arg-type]


def _write_bed(path: Path, rows: list[tuple[str, int, int, str]]) -> None:
    path.write_text("".join(f"{c}\t{s}\t{e}\t{n}\n" for c, s, e, n in rows))


def _write_map(path: Path, rows: list[MapRow]) -> None:
    write_map(rows, path)


def test_remap_renames_and_drops_unscaffolded(tmp_path: Path) -> None:
    # ctg1 is in the map (-> renamed); ctg_short is not (-> dropped).
    bed = tmp_path / "hap1.KS_cytoband.cytoband.smoothed.bed"
    _write_bed(bed, [("ctg1", 0, 100, "1p36.1"), ("ctg1", 100, 200, "1p36.2"), ("ctg_short", 0, 50, "x")])
    mp = tmp_path / "hap1.KS_v2.scaffold_map.tsv"
    _write_map(mp, [_row()])

    out = tmp_path / "hap1.cytoband.scaffolded.bed"
    stats = remap_bed_with_map(bed, out, mp)

    lines = [ln for ln in out.read_text().splitlines() if ln]
    assert all(ln.startswith("chr1_hap1_ctg1\t") for ln in lines)
    assert len(lines) == 2  # ctg_short dropped
    assert stats.mapped_contigs == 1
    assert stats.dropped_contigs == 1
    assert stats.bed_contigs == 2


def test_remap_mirrors_flipped_contig(tmp_path: Path) -> None:
    bed = tmp_path / "hap1.x.bed"
    _write_bed(bed, [("ctg1", 0, 100, "a"), ("ctg1", 900_000, 1_000_000, "b")])
    mp = tmp_path / "hap1.scaffold_map.tsv"
    _write_map(mp, [_row(flipped=True, new_name="chr1_hap1_ctg1_rc")])

    out = tmp_path / "out.bed"
    remap_bed_with_map(bed, out, mp)
    rows = [ln.split("\t") for ln in out.read_text().splitlines() if ln]
    # length=1_000_000: [0,100) -> [999_900, 1_000_000); [900_000,1_000_000) -> [0,100_000)
    # and reversed so output stays coordinate-sorted.
    assert rows[0][:3] == ["chr1_hap1_ctg1_rc", "0", "100000"]
    assert rows[1][:3] == ["chr1_hap1_ctg1_rc", "999900", "1000000"]


def test_zero_overlap_is_hard_error(tmp_path: Path) -> None:
    bed = tmp_path / "hap1.bed"
    _write_bed(bed, [("other_ctg", 0, 100, "a")])
    mp = tmp_path / "hap1.scaffold_map.tsv"
    _write_map(mp, [_row()])
    with pytest.raises(ScaffoldError, match="do not appear to describe the same assembly"):
        remap_bed_with_map(bed, tmp_path / "out.bed", mp)


def test_interval_beyond_contig_length_is_hard_error(tmp_path: Path) -> None:
    bed = tmp_path / "hap1.bed"
    _write_bed(bed, [("ctg1", 0, 2_000_000, "a")])  # end > map length 1_000_000
    mp = tmp_path / "hap1.scaffold_map.tsv"
    _write_map(mp, [_row()])
    with pytest.raises(ScaffoldError, match="beyond the contig lengths"):
        remap_bed_with_map(bed, tmp_path / "out.bed", mp)


def test_stem_mismatch_warns_but_succeeds(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    # BED filename stem 'sampleB' != map source FASTA stem 'hap1'.
    bed = tmp_path / "sampleB.cytoband.smoothed.bed"
    _write_bed(bed, [("ctg1", 0, 100, "a")])
    mp = tmp_path / "hap1.scaffold_map.tsv"
    _write_map(mp, [_row()])
    with caplog.at_level("WARNING"):
        remap_bed_with_map(bed, tmp_path / "out.bed", mp)
    assert any("same assembly" in r.message for r in caplog.records)


def test_strict_promotes_stem_mismatch_to_error(tmp_path: Path) -> None:
    bed = tmp_path / "sampleB.bed"
    _write_bed(bed, [("ctg1", 0, 100, "a")])
    mp = tmp_path / "hap1.scaffold_map.tsv"
    _write_map(mp, [_row()])
    with pytest.raises(ScaffoldError):
        remap_bed_with_map(bed, tmp_path / "out.bed", mp, strict=True)


def test_strict_errors_on_map_contig_absent_from_bed(tmp_path: Path) -> None:
    bed = tmp_path / "hap1.bed"
    _write_bed(bed, [("ctg1", 0, 100, "a")])  # ctg2 (in map) has no BED records
    mp = tmp_path / "hap1.scaffold_map.tsv"
    _write_map(mp, [_row(), _row(original_name="ctg2", new_name="chr2_hap1_ctg2", chromosome="chr2")])
    with pytest.raises(ScaffoldError, match="no records"):
        remap_bed_with_map(bed, tmp_path / "out.bed", mp, strict=True)


def test_empty_map_is_error(tmp_path: Path) -> None:
    bed = tmp_path / "hap1.bed"
    _write_bed(bed, [("ctg1", 0, 100, "a")])
    mp = tmp_path / "hap1.scaffold_map.tsv"
    _write_map(mp, [])
    with pytest.raises(ScaffoldError, match="empty"):
        remap_bed_with_map(bed, tmp_path / "out.bed", mp)


def test_gzip_output_inferred_from_suffix(tmp_path: Path) -> None:
    import gzip

    bed = tmp_path / "hap1.bed"
    _write_bed(bed, [("ctg1", 0, 100, "a")])
    mp = tmp_path / "hap1.scaffold_map.tsv"
    _write_map(mp, [_row()])
    out = tmp_path / "out.bed.gz"
    remap_bed_with_map(bed, out, mp)
    with gzip.open(out, "rt") as h:
        assert h.read().startswith("chr1_hap1_ctg1\t")


# --- CLI ------------------------------------------------------------


def test_cli_remap_bed(tmp_path: Path) -> None:
    bed = tmp_path / "hap1.cytoband.smoothed.bed"
    _write_bed(bed, [("ctg1", 0, 100, "1p36.1"), ("drop", 0, 10, "x")])
    mp = tmp_path / "hap1.scaffold_map.tsv"
    _write_map(mp, [_row()])
    out = tmp_path / "hap1.cytoband.scaffolded.bed"

    result = CliRunner().invoke(
        main, ["remap-bed", "-b", str(bed), "-m", str(mp), "-o", str(out)]
    )
    assert result.exit_code == 0, result.output
    assert out.is_file()
    assert "1/2 contigs placed, 1 dropped" in result.output


def test_cli_remap_bed_mismatch_errors(tmp_path: Path) -> None:
    bed = tmp_path / "hap1.bed"
    _write_bed(bed, [("nope", 0, 100, "a")])
    mp = tmp_path / "hap1.scaffold_map.tsv"
    _write_map(mp, [_row()])
    result = CliRunner().invoke(
        main, ["remap-bed", "-b", str(bed), "-m", str(mp), "-o", str(tmp_path / "o.bed")]
    )
    assert result.exit_code != 0
    assert "same assembly" in result.output
