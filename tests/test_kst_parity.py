"""Parity of the four smoothers on fixed real tracks (plan 5.1).

Smoothing exists in four places: ``hks smooth``, KaryoScope's Python
smoother (:mod:`karyoscope.core.smooth`, the KMC backend's), and karyozip's
``decompress --smooth`` in its C and Rust builds. ``annotate --format kst``
stores only the presmoothed track and relies on karyozip to derive the
smoothed one on read, so they must agree byte for byte.

The fixture (``tests/data/kst_parity``) is 122 HG002 reads of a real
``annotate`` run against ``HKS_human_CHM13_v2``: every feature set's
presmoothed BED and the smoothed BED ``annotate`` wrote from it with
``hks smooth``, read names shortened to ``r001``.. (see the README there).
Each smoother is run on the presmoothed BED and compared with that smoothed
BED.

Which karyozip builds are checked: ``karyozip`` on ``$PATH``, plus
``karyozip-c`` if that is on ``$PATH`` too, plus any listed (colon
separated) in ``$KARYOSCOPE_TEST_KARYOZIP``.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from karyoscope.core.io.features import NOVEL_NAME
from karyoscope.core.io.hierarchy import parse_hierarchy
from karyoscope.core.smooth import (
    HierarchyIndex,
    Interval,
    _render_for_output,
    merge_adjacent,
    smooth_intervals,
)

DATA = Path(__file__).resolve().parent / "data" / "kst_parity"
SETS = ("acrocentric", "chromosome", "gene", "region", "repeat", "subtelomeric")


def _karyozip_builds() -> list[str]:
    found: list[str] = []
    for name in ("karyozip", "karyozip-c"):
        path = shutil.which(name)
        if path:
            found.append(path)
    for extra in os.environ.get("KARYOSCOPE_TEST_KARYOZIP", "").split(":"):
        if extra and extra not in found:
            found.append(extra)
    return found


def _text(path: Path) -> str:
    with gzip.open(path, "rt") as h:
        return h.read()


@pytest.fixture(scope="module")
def plain_presmoothed(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """Each set's presmoothed BED, decompressed once (hks and karyozip read files)."""
    d = tmp_path_factory.mktemp("kst_parity")
    out: dict[str, Path] = {}
    for fs in SETS:
        p = d / f"{fs}.presmoothed.bed"
        p.write_text(_text(DATA / f"{fs}.presmoothed.bed.gz"))
        out[fs] = p
    return out


def _python_smooth(presmoothed: str, fs: str) -> str:
    """KaryoScope's Python smoother over a per-set presmoothed BED.

    The KMC path feeds it the unmerged id runs of the combined BED; here the
    input is the presmoothed (already merged) track, which is what every
    other smoother reads, so the comparison is like for like.
    """
    index = HierarchyIndex.from_hierarchy(parse_hierarchy(DATA / "hierarchy.tsv"), fs)
    root = index.root
    out: list[str] = []
    seq: str | None = None
    ivs: list[Interval] = []

    def flush() -> None:
        if ivs:
            for iv in merge_adjacent(smooth_intervals(ivs, index), root):
                out.append(_render_for_output(iv, root))

    for line in presmoothed.splitlines():
        name, start, end, label = line.split("\t")
        if name != seq:
            flush()
            seq, ivs = name, []
        novel = label == NOVEL_NAME
        ivs.append(Interval(name, int(start), int(end), root if novel else label, novel))
    flush()
    return "".join(out)


def _first_difference(got: str, want: str) -> str:
    g, w = got.splitlines(), want.splitlines()
    for i, (a, b) in enumerate(zip(g, w, strict=False), start=1):
        if a != b:
            return f"line {i}: got {a!r}, expected {b!r}"
    return f"got {len(g)} lines, expected {len(w)}"


@pytest.mark.parametrize("fs", SETS)
def test_python_smoother_matches_annotate(fs: str) -> None:
    want = _text(DATA / f"{fs}.smoothed.bed.gz")
    got = _python_smooth(_text(DATA / f"{fs}.presmoothed.bed.gz"), fs)
    assert got == want, f"Python smoother differs for {fs}: {_first_difference(got, want)}"


@pytest.mark.integration
@pytest.mark.parametrize("fs", SETS)
def test_hks_smooth_matches_annotate(
    fs: str, plain_presmoothed: dict[str, Path], tmp_path: Path
) -> None:
    if shutil.which("hks") is None and not os.environ.get("KARYOSCOPE_HKS"):
        pytest.skip("hks is not on PATH")
    from karyoscope.core.io.hks import run_hks_smooth

    out = tmp_path / f"{fs}.smoothed.bed"
    run_hks_smooth(
        hierarchy_file=DATA / "hks" / f"{fs}.hierarchy.txt",
        input_path=plain_presmoothed[fs],
        output_path=out,
        threads=2,
        capture=True,
    )
    got, want = out.read_text(), _text(DATA / f"{fs}.smoothed.bed.gz")
    assert got == want, f"hks smooth differs for {fs}: {_first_difference(got, want)}"


@pytest.mark.integration
@pytest.mark.parametrize("binary", _karyozip_builds() or ["karyozip"])
@pytest.mark.parametrize("fs", SETS)
def test_karyozip_smooth_matches_annotate(
    fs: str, binary: str, plain_presmoothed: dict[str, Path], tmp_path: Path
) -> None:
    """Per-set ``.kst``: round trip, ``decompress --smooth`` and ``query --smooth``."""
    if shutil.which(binary) is None:
        pytest.skip(f"{binary} is not on PATH")
    kst = tmp_path / f"{fs}.kst"
    subprocess.run(
        [binary, "compress", "-i", plain_presmoothed[fs], "-o", kst, "-e", "low"],
        check=True,
        capture_output=True,
    )
    pre = subprocess.run(
        [binary, "decompress", "-i", kst], check=True, capture_output=True, text=True
    ).stdout
    assert pre == plain_presmoothed[fs].read_text(), f"{binary}: {fs} round trip differs"

    hier = DATA / "hierarchy.tsv"
    smo = subprocess.run(
        [binary, "decompress", "-i", kst, "--smooth", "--hier", hier, "--set", fs],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    want = _text(DATA / f"{fs}.smoothed.bed.gz")
    assert smo == want, f"{binary} --smooth differs for {fs}: {_first_difference(smo, want)}"


@pytest.mark.integration
@pytest.mark.parametrize("binary", _karyozip_builds() or ["karyozip"])
def test_karyozip_merged_track_matches_per_set_beds(
    binary: str, plain_presmoothed: dict[str, Path], tmp_path: Path
) -> None:
    """One merged ``.kst`` (``-d``) decodes to every set's presmoothed and smoothed BED.

    ``features.tsv`` here is built from the fixture itself: one featureID per
    distinct label tuple across the six sets, which is what a merged track
    stores. The real database's table is 690 KB; the tuples are what matter.
    """
    if shutil.which(binary) is None:
        pytest.skip(f"{binary} is not on PATH")
    from karyoscope.core.io.track import open_track_text

    features = _fixture_features(plain_presmoothed, tmp_path / "features.tsv")
    columns = features.read_text().splitlines()[0].split("\t")[1:]
    kst = tmp_path / "merged.kst"
    cmd = [binary, "compress", "-d", features, "--hier", DATA / "hierarchy.tsv"]
    for fs in columns:
        cmd += ["-i", plain_presmoothed[fs]]
    subprocess.run([*cmd, "-o", kst, "-e", "low"], check=True, capture_output=True)

    db = _fixture_db(tmp_path / "db", features)
    for fs in SETS:
        with open_track_text(kst, feature_set=fs, db_dir=db) as h:
            assert h.read() == plain_presmoothed[fs].read_text(), f"{fs} presmoothed"
        with open_track_text(kst, feature_set=fs, db_dir=db, smoothed=True) as h:
            got = h.read()
        want = _text(DATA / f"{fs}.smoothed.bed.gz")
        assert got == want, f"{fs} smoothed: {_first_difference(got, want)}"


def _fixture_features(beds: dict[str, Path], out: Path) -> Path:
    """A ``features.tsv`` with one row per label tuple the six tracks hold.

    Walks all six BEDs in step, cutting each read at every boundary of any
    set, the way a k-mer's featureID covers all sets at once.
    """
    import itertools

    tuples: dict[tuple[str, ...], int] = {}
    handles = {fs: beds[fs].open() for fs in SETS}
    try:
        streams = {fs: (ln.rstrip("\n").split("\t") for ln in handles[fs]) for fs in SETS}
        rows = {fs: next(streams[fs], None) for fs in SETS}
        while rows[SETS[0]] is not None:
            name = rows[SETS[0]][0]
            pos = int(rows[SETS[0]][1])
            while all(r is not None and r[0] == name for r in rows.values()):
                labels = tuple(rows[fs][3] for fs in SETS)
                if NOVEL_NAME not in labels:
                    tuples.setdefault(labels, len(tuples) + 1)
                elif any(lab != NOVEL_NAME for lab in labels):
                    raise AssertionError(f"{name}:{pos} is novel in some sets only")
                end = min(int(rows[fs][2]) for fs in SETS)
                pos = end
                for fs in SETS:
                    if int(rows[fs][2]) == end:
                        rows[fs] = next(streams[fs], None)
    finally:
        for h in handles.values():
            h.close()
    lines = ["featureID\t" + "\t".join(SETS)]
    lines += [f"{i}\t" + "\t".join(t) for t, i in itertools.islice(tuples.items(), None)]
    out.write_text("\n".join(lines) + "\n")
    return out


def _fixture_db(db: Path, features: Path) -> Path:
    """The three files :func:`open_track_text` reads from a database directory."""
    db.mkdir()
    shutil.copy(features, db / "features.tsv")
    shutil.copy(DATA / "hierarchy.tsv", db / "hierarchy.tsv")
    sets = "".join(f"  - {fs}\n" for fs in SETS)
    (db / "manifest.yaml").write_text(
        "id: kst_parity\n"
        'version: "1.0.0"\n'
        'karyoscope_min_version: "2.0.0"\n'
        "index:\n  type: hks\n  basename: index/features\n"
        "hierarchy: hierarchy.tsv\nfeatures: features.tsv\ncolors: colors.tsv\n"
        "kmer:\n  size: 31\n  type: fixed\n  max_size: 31\n"
        f"feature_sets:\n{sets}"
    )
    return db
