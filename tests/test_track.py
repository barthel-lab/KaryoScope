"""Tests for the shared annotation-track reader and the karyozip wrapper.

:mod:`karyoscope.core.io.track` is the one reader every command uses for
annotate output, whether BED, bgzipped BED or ``.kst``. The ``.kst`` stream
semantics (errors surface, early stops do not) are tested against a fake
``karyozip`` script, so they run in the default pass; the round trips through
the real binary are ``integration``.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from karyoscope.core.external import ExternalToolError
from karyoscope.core.io import karyozip as kmz
from karyoscope.core.io.track import (
    Track,
    annotation_track,
    as_track,
    is_plain_bed,
    open_track,
    open_track_text,
)
from karyoscope.exceptions import KaryoscopeError

requires_karyozip = pytest.mark.skipif(
    shutil.which("karyozip") is None, reason="karyozip not on PATH"
)

ROWS = "a\t0\t5\tchr1\na\t5\t9\tnovel\nb\t0\t3\tchr2\n"


# --- Track and lookup --------------------------------------------------


def test_track_is_path_like(tmp_path: Path) -> None:
    p = tmp_path / "x.kst"
    p.write_bytes(b"")
    t = Track(p, smoothed=True, feature_set="chromosome")
    assert os.fspath(t) == str(p)
    assert t.name == "x.kst" and t.suffix == ".kst" and t.parent == tmp_path
    assert t.is_file()
    assert "smoothed chromosome" in str(t)
    assert str(Track(tmp_path / "x.bed")) == str(tmp_path / "x.bed")


def test_as_track_overrides_only_what_is_given(tmp_path: Path) -> None:
    base = Track(tmp_path / "x.kst", smoothed=True, feature_set="a", db_dir=tmp_path)
    assert as_track(base) == base
    assert as_track(base, smoothed=False).smoothed is False
    assert as_track(base, feature_set="b").feature_set == "b"
    assert as_track(str(tmp_path / "y.bed")).path == tmp_path / "y.bed"


def test_is_plain_bed() -> None:
    assert is_plain_bed(Path("a.bed"))
    assert not is_plain_bed(Path("a.bed.gz"))
    assert not is_plain_bed(Path("a.kst"))
    assert not is_plain_bed(Track(Path("a.kst")))
    assert not is_plain_bed("-")


def test_annotation_track_prefers_bed_then_per_set_then_merged_kst(tmp_path: Path) -> None:
    def find(variant: str = "smoothed"):
        return annotation_track(tmp_path, "s", "DB", "chromosome", variant, db_dir=tmp_path)

    assert find() == tmp_path / "s.DB.chromosome.smoothed.bed.gz"  # nothing on disk yet
    merged = tmp_path / "s.DB.kst"
    merged.write_bytes(b"")
    assert find() == Track(merged, smoothed=True, feature_set="chromosome", db_dir=tmp_path)
    assert find("presmoothed") == Track(
        merged, smoothed=False, feature_set="chromosome", db_dir=tmp_path
    )
    per_set = tmp_path / "s.DB.chromosome.kst"
    per_set.write_bytes(b"")
    assert find().path == per_set
    plain = tmp_path / "s.DB.chromosome.smoothed.bed"
    plain.write_text(ROWS)
    assert find() == plain
    gz = tmp_path / "s.DB.chromosome.smoothed.bed.gz"
    gz.write_bytes(gzip.compress(ROWS.encode()))
    assert find() == gz


# --- BED reading -------------------------------------------------------


@pytest.mark.parametrize("suffix", [".bed", ".bed.gz"])
def test_open_track_reads_bed_and_gz(tmp_path: Path, suffix: str) -> None:
    p = tmp_path / f"x{suffix}"
    data = "# header\n" + ROWS
    p.write_bytes(gzip.compress(data.encode()) if suffix.endswith(".gz") else data.encode())
    # The kst-only options are ignored for a BED.
    rows = list(open_track(p, smoothed=True, feature_set="whatever"))
    assert rows == [("a", 0, 5, "chr1"), ("a", 5, 9, "novel"), ("b", 0, 3, "chr2")]
    with open_track_text(p) as h:
        assert h.read() == data


def test_open_track_reports_malformed_rows(tmp_path: Path) -> None:
    p = tmp_path / "x.bed"
    p.write_text("a\t0\t5\n")
    with pytest.raises(KaryoscopeError, match="expected 4"):
        list(open_track(p))
    p.write_text("a\tx\t5\tchr1\n")
    with pytest.raises(KaryoscopeError, match="non-integer"):
        list(open_track(p))


# --- .kst header and command lines --------------------------------------


def _kst_header(path: Path, *, merged: bool, magic: bytes = b"KST00001") -> Path:
    flags = (1 << 58) if merged else 0
    path.write_bytes(magic + (0).to_bytes(8, "little") + flags.to_bytes(8, "little") + b"\0" * 64)
    return path


def test_header_says_whether_a_track_holds_feature_ids(tmp_path: Path) -> None:
    assert kmz.kst_holds_feature_ids(_kst_header(tmp_path / "m.kst", merged=True))
    assert not kmz.kst_holds_feature_ids(_kst_header(tmp_path / "s.kst", merged=False))
    (tmp_path / "short.kst").write_bytes(b"KST")
    with pytest.raises(kmz.KstError, match="too short"):
        kmz.kst_holds_feature_ids(tmp_path / "short.kst")


@pytest.fixture
def fake_binary(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(kmz, "karyozip_binary", lambda: "karyozip")


def test_decompress_args(tmp_path: Path, fake_binary: None) -> None:
    f, h = tmp_path / "features.tsv", tmp_path / "hierarchy.tsv"
    p = tmp_path / "x.kst"
    assert kmz.decompress_args(p, merged=False) == ["karyozip", "decompress", "-i", str(p)]
    assert kmz.decompress_args(p, merged=False, smoothed=True, hierarchy=h, feature_set="a") == [
        "karyozip", "decompress", "-i", str(p), "--hier", str(h), "--set", "a", "--smooth",
    ]  # fmt: skip
    assert kmz.decompress_args(
        p, merged=True, smoothed=True, feature_set="a", features=f, hierarchy=h
    ) == [
        "karyozip", "decompress", "-i", str(p), "-d", str(f), "--hier", str(h), "--set", "a",
        "--smooth",
    ]  # fmt: skip
    with pytest.raises(kmz.KstError, match="say which one"):
        kmz.decompress_args(p, merged=True, features=f)
    with pytest.raises(kmz.KstError, match=r"features\.tsv"):
        kmz.decompress_args(p, merged=True, feature_set="a")
    with pytest.raises(kmz.KstError, match=r"hierarchy\.tsv"):
        kmz.decompress_args(p, merged=False, smoothed=True)


def test_parse_karyozip_version() -> None:
    assert kmz.parse_karyozip_version("karyozip 0.2.0, format KST00001") == "0.2.0"
    assert kmz.parse_karyozip_version(None) is None
    assert kmz.parse_karyozip_version("something else") is None


# --- .kst streaming, against a fake karyozip ----------------------------


@pytest.fixture
def fake_karyozip(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Put a scripted ``karyozip`` first on PATH; return a setter for its body."""
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    script = bindir / "karyozip"
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")

    def set_body(body: str) -> None:
        script.write_text("#!/bin/sh\n" + body + "\n")
        script.chmod(0o755)

    return set_body


def test_kst_stream_yields_karyozip_output(tmp_path: Path, fake_karyozip) -> None:
    fake_karyozip(f"printf '{ROWS.encode('unicode_escape').decode()}'")
    kst = _kst_header(tmp_path / "x.kst", merged=False)
    assert list(open_track(kst)) == [("a", 0, 5, "chr1"), ("a", 5, 9, "novel"), ("b", 0, 3, "chr2")]


def test_kst_stream_failure_is_an_error_not_a_short_track(tmp_path: Path, fake_karyozip) -> None:
    fake_karyozip("printf 'a\\t0\\t5\\tchr1\\n'; echo 'karyozip: corrupt block' >&2; exit 1")
    kst = _kst_header(tmp_path / "x.kst", merged=False)
    with pytest.raises(ExternalToolError, match="corrupt block"):
        list(open_track(kst))


def test_kst_stream_early_stop_is_not_an_error(tmp_path: Path, fake_karyozip) -> None:
    # Far more output than a pipe buffer: the child is still writing when the
    # reader stops, and must be stopped rather than reported as failed.
    fake_karyozip("yes 'a\t0\t5\tchr1' | head -n 2000000; exit 1")
    kst = _kst_header(tmp_path / "x.kst", merged=False)
    with open_track_text(kst) as h:
        assert h.readline() == "a\t0\t5\tchr1\n"


def test_kst_stream_needs_the_database_for_a_merged_track(tmp_path: Path, fake_karyozip) -> None:
    fake_karyozip("exit 0")
    kst = _kst_header(tmp_path / "m.kst", merged=True)
    with pytest.raises(kmz.KstError, match="say which one"):
        list(open_track(kst))
    with pytest.raises(kmz.KstError, match=r"features\.tsv"):
        list(open_track(kst, feature_set="chromosome"))


def test_missing_kst_is_reported(tmp_path: Path) -> None:
    with pytest.raises(KaryoscopeError, match="not found"):
        list(open_track(tmp_path / "absent.kst"))


# --- the real karyozip ---------------------------------------------------


def _karyozip(*args: object) -> str:
    return subprocess.run(
        ["karyozip", *map(str, args)], check=True, capture_output=True, text=True
    ).stdout


@pytest.mark.integration
@requires_karyozip
@pytest.mark.parametrize(
    "names",
    [
        pytest.param(["read_a", "read_b", "read_c"], id="names"),
        pytest.param(["0", "1", "2"], id="ranks"),
        pytest.param(["5", "6", "7"], id="ranks-offset"),
    ],
)
def test_names_and_ranks_round_trip(tmp_path: Path, names: list[str]) -> None:
    """``--query-names`` writes names, reads write ranks: both come back exactly."""
    bed = tmp_path / "x.bed"
    bed.write_text(
        "".join(f"{n}\t0\t10\tchr1\n{n}\t10\t12\tnovel\n{n}\t12\t40\tchr2\n" for n in names)
    )
    kst = tmp_path / "x.kst"
    kmz.compress(inputs=[bed], output=kst)
    assert not kmz.kst_holds_feature_ids(kst)
    with open_track_text(kst) as h:
        assert h.read() == bed.read_text()
    assert kmz.index(kst).is_file()
    assert _karyozip("query", "-i", kst, "-r", names[1]) == "".join(
        line + "\n" for line in bed.read_text().splitlines()[3:6]
    )


@pytest.mark.integration
@requires_karyozip
def test_real_merged_track_reads_back_per_set(unpacked_dummy_db: Path, tmp_path: Path) -> None:
    """``compress -d`` of KMC-style featureIDs; each set decodes, raw and smoothed."""
    combined = tmp_path / "combined.bed"
    # dummy features.tsv: 1 -> chr1/rA, 2 -> chr1/rB, 3 -> chr2/rC; 0 is novel.
    combined.write_text(
        "seqA\t0\t10\t1\nseqA\t10\t14\t0\nseqA\t14\t30\t1\nseqA\t30\t40\t2\nseqB\t0\t20\t3\n"
    )
    kst = tmp_path / "m.kst"
    kmz.compress(
        inputs=[combined],
        output=kst,
        features=unpacked_dummy_db / "features.tsv",
        hierarchy=unpacked_dummy_db / "hierarchy.tsv",
    )
    assert kmz.kst_holds_feature_ids(kst)
    chrom = list(open_track(kst, feature_set="chromosome", db_dir=unpacked_dummy_db))
    assert chrom == [
        ("seqA", 0, 10, "chr1"),
        ("seqA", 10, 14, "novel"),
        ("seqA", 14, 40, "chr1"),
        ("seqB", 0, 20, "chr2"),
    ]
    region = list(open_track(kst, feature_set="region", db_dir=unpacked_dummy_db))
    assert [r[3] for r in region] == ["rA", "novel", "rA", "rB", "rC"]
    smoothed = list(
        open_track(kst, feature_set="chromosome", db_dir=unpacked_dummy_db, smoothed=True)
    )
    assert smoothed == [("seqA", 0, 40, "chr1"), ("seqB", 0, 20, "chr2")]
