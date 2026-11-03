"""``annotate --format kst``: planning, validation, and runs through real karyozip.

The flag checks, the merged-versus-per-set decision, the dependency list and
the disk estimate are pure and run in the default pass. The end-to-end runs
stub only the k-mer backends (``get_featureIDs``, ``hks lookup``) and drive
the real ``karyozip``, so they are ``integration``: they check that every
``.kst`` decodes to exactly the BED ``--format bed`` writes, and that the
downstream readers (``bin``, ``remap-bed``, the scaffold cascade's lookup)
take a ``.kst`` where they took a BED.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

import karyoscope.core.annotate as ann
from karyoscope.cli import main
from karyoscope.core.annotate import annotate, annotate_batch
from karyoscope.core.io.kmc import combined_bed_path, write_combined_marker
from karyoscope.core.io.track import Track, open_track, open_track_text
from karyoscope.exceptions import KaryoscopeError

DUMMY_DB_ID = "KS_dummy_test_v1"
HKS_DB_ID = "HKS_kst_test_v1"

requires_karyozip = pytest.mark.skipif(
    shutil.which("karyozip") is None, reason="karyozip not on PATH"
)
requires_bgzip = pytest.mark.skipif(shutil.which("bgzip") is None, reason="bgzip not on PATH")

#: Combined-BED content for the dummy db (features.tsv: 1 -> chr1/rA,
#: 2 -> chr1/rB, 3 -> chr2/rC; id 0 is the novel sentinel).
_COMBINED_LINES = (
    "seqA\t0\t10\t1\nseqA\t10\t14\t0\nseqA\t14\t30\t1\nseqA\t30\t40\t2\nseqB\t0\t20\t3\n"
)


# --- pure checks -------------------------------------------------------


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"output_format": "bed", "smoothed_bed": True}, "--smoothed-bed only applies"),
        ({"output_format": "kst", "keep_presmoothed": False}, "--no-keep-presmoothed"),
        ({"output_format": "kst", "smoothed_bed": True, "smooth": False}, "--no-smooth"),
        ({"output_format": "kst", "kst_effort": "high"}, "--kst-effort"),
        ({"output_format": "zip"}, "unknown output format"),
    ],
)
def test_format_flag_combinations_are_checked(kwargs: dict, match: str) -> None:
    args = {
        "output_format": "bed",
        "kst_effort": "low",
        "smooth": True,
        "keep_presmoothed": True,
        "smoothed_bed": False,
    }
    args.update(kwargs)
    with pytest.raises(KaryoscopeError, match=match):
        ann._check_format_args(**args)


def test_valid_format_flags_pass() -> None:
    for fmt in ("bed", "kst"):
        ann._check_format_args(
            output_format=fmt, kst_effort="low", smooth=True, keep_presmoothed=True,
            smoothed_bed=fmt == "kst",
        )  # fmt: skip


def test_karyozip_is_a_dependency_only_for_kst(tmp_path: Path) -> None:
    fa = tmp_path / "q.fa"
    deps = ann._annotate_dependencies(index_type="hks", input_path=fa, bgzip=True)
    assert "karyozip" not in deps and "bgzip" in deps
    deps = ann._annotate_dependencies(
        index_type="hks", input_path=fa, bgzip=True, output_format="kst"
    )
    assert "karyozip" in deps and "bgzip" not in deps  # no BED to compress
    deps = ann._annotate_dependencies(
        index_type="hks", input_path=fa, bgzip=True, output_format="kst", smoothed_bed=True
    )
    assert "bgzip" in deps


def test_kst_disk_estimate_is_the_presmoothed_peak() -> None:
    bases = 1_000_000_000
    pre = int(bases * ann.PRESMOOTHED_BYTES_PER_BASE)
    est = ann.estimate_output_bytes
    common = {"input_bases": bases, "n_feature_sets": 6, "keep_presmoothed": True, "smooth": True}
    # HKS merged keeps every set's BED until the end; per-set compresses as it goes.
    assert est(**common, index_type="hks", output_format="kst", kst_merged=True) == 6 * pre
    assert est(**common, index_type="hks", output_format="kst") == pre
    with_smoothed = est(
        **common, index_type="hks", output_format="kst", kst_merged=True, smoothed_bed=True
    )
    assert with_smoothed == int(6 * pre + bases * ann.SMOOTHED_BYTES_PER_BASE * 6)
    # KMC merged compresses its combined BED directly: only the transient.
    kmc = est(**common, index_type="kmc", output_format="kst", kst_merged=True)
    assert kmc == int(bases * ann.PRESMOOTHED_BYTES_PER_BASE * ann.TRANSIENT_INTERMEDIATE_FACTOR)
    # And BED is unchanged and larger.
    assert est(**common, index_type="hks") > est(**common, index_type="hks", output_format="kst")


def test_plan_kst_merges_only_when_every_set_and_features_tsv(unpacked_dummy_db: Path) -> None:
    from karyoscope.manifest import validate_database_layout

    manifest = validate_database_layout(unpacked_dummy_db)
    plan = ann._plan_kst(
        manifest=manifest,
        db_dir=unpacked_dummy_db,
        requested=["region", "chromosome"],
        effort="low",
        smoothed_bed=False,
        bgzip=True,
    )
    assert plan.merged and plan.columns == ("chromosome", "region")
    assert plan.features == unpacked_dummy_db / "features.tsv"
    sub = ann._plan_kst(
        manifest=manifest,
        db_dir=unpacked_dummy_db,
        requested=["chromosome"],
        effort="low",
        smoothed_bed=False,
        bgzip=True,
    )
    assert not sub.merged and sub.features is None
    manifest.features = None
    assert not ann._plan_kst(
        manifest=manifest,
        db_dir=unpacked_dummy_db,
        requested=["chromosome", "region"],
        effort="low",
        smoothed_bed=False,
        bgzip=True,
    ).merged


def test_cli_rejects_smoothed_bed_without_kst(cli_runner: CliRunner, tmp_path: Path) -> None:
    fa = tmp_path / "q.fa"
    fa.write_text(">a\nACGT\n")
    db_root = tmp_path / "db"
    db_root.mkdir()
    result = cli_runner.invoke(
        main, ["annotate", "-i", str(fa), "--db-root", str(db_root), "--smoothed-bed"]
    )
    assert result.exit_code != 0


def test_cli_help_lists_the_kst_options(cli_runner: CliRunner) -> None:
    result = cli_runner.invoke(main, ["annotate", "--help"])
    assert result.exit_code == 0
    for opt in ("--format [bed|kst]", "--kst-effort [low|medium]", "--smoothed-bed"):
        assert opt in result.output


# --- KMC backend, real karyozip ------------------------------------------


@pytest.fixture
def query_fasta(tmp_path: Path) -> Path:
    fa = tmp_path / "q.fasta"
    fa.write_text(">seqA\n" + "A" * 40 + "\n>seqB\n" + "C" * 20 + "\n")
    return fa


@pytest.fixture
def stub_kmc(monkeypatch: pytest.MonkeyPatch) -> dict:
    """Stub get_featureIDs (writes the fixed combined BED) and the preflight."""
    calls: dict = {"get_featureids": 0}
    monkeypatch.setattr(ann.preflight, "require", lambda *a, **kw: None)

    def fake_get_featureids(*, db_path, input_path, output_dir, threads, prefix, **kw):
        calls["get_featureids"] += 1
        bed = combined_bed_path(output_dir, prefix)
        bed.write_text(_COMBINED_LINES)
        write_combined_marker(bed, prefix=prefix, db_path=db_path, input_path=input_path)
        return bed

    monkeypatch.setattr(ann, "run_get_featureids", fake_get_featureids)
    return calls


def _text(path: Path) -> str:
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as h:
            return h.read()
    return path.read_text()


@pytest.mark.integration
@requires_karyozip
@requires_bgzip
def test_kmc_merged_kst_decodes_to_the_bed_run(
    populated_db_root: Path, query_fasta: Path, tmp_path: Path, stub_kmc: dict
) -> None:
    bed_out, kst_out = tmp_path / "bed", tmp_path / "kst"
    beds = annotate(
        input_path=query_fasta, output_dir=bed_out, db_root=populated_db_root, threads=1
    )
    result = annotate(
        input_path=query_fasta,
        output_dir=kst_out,
        db_root=populated_db_root,
        threads=1,
        output_format="kst",
        smoothed_bed=True,
    )

    merged = kst_out / f"q.{DUMMY_DB_ID}.kst"
    assert result.kst_merged
    assert result.kst_paths == {"chromosome": merged, "region": merged}
    assert result.presmoothed_paths == {}
    assert (kst_out / f"q.{DUMMY_DB_ID}.kst.kzi").is_file()
    db_dir = populated_db_root / DUMMY_DB_ID
    for fs in ("chromosome", "region"):
        with open_track_text(merged, feature_set=fs, db_dir=db_dir) as h:
            assert h.read() == _text(beds.presmoothed_paths[fs]), fs
        with open_track_text(merged, feature_set=fs, db_dir=db_dir, smoothed=True) as h:
            assert h.read() == _text(beds.smoothed_paths[fs]), fs
        # --smoothed-bed wrote the same smoothed BED, bgzipped.
        assert result.smoothed_paths[fs].name == f"q.{DUMMY_DB_ID}.{fs}.smoothed.bed.gz"
        assert _text(result.smoothed_paths[fs]) == _text(beds.smoothed_paths[fs])

    # Nothing else is left behind: no presmoothed BED, no combined intermediate.
    assert sorted(p.name for p in kst_out.iterdir()) == sorted(
        [
            merged.name,
            merged.name + ".kzi",
            f"q.{DUMMY_DB_ID}.annotate_manifest.json",
            f"q.{DUMMY_DB_ID}.chromosome.smoothed.bed.gz",
            f"q.{DUMMY_DB_ID}.region.smoothed.bed.gz",
        ]
    )
    doc = json.loads(result.manifest_path.read_text())
    assert doc["format"] == "kst" and doc["kst"]["merged"] and doc["kst"]["effort"] == "low"
    assert doc["kst"]["karyozip_version"].startswith("karyozip ")
    by_name = {o["path"]: o for o in doc["outputs"]}
    assert by_name[merged.name]["sha256"] == hashlib.sha256(merged.read_bytes()).hexdigest()
    assert by_name[merged.name]["feature_sets"] == ["chromosome", "region"]
    assert by_name[merged.name + ".kzi"]["kind"] == "kst_index"


@pytest.mark.integration
@requires_karyozip
def test_kmc_feature_set_subset_writes_per_set_kst(
    populated_db_root: Path, query_fasta: Path, tmp_path: Path, stub_kmc: dict
) -> None:
    out = tmp_path / "out"
    result = annotate(
        input_path=query_fasta,
        output_dir=out,
        db_root=populated_db_root,
        feature_sets=["chromosome"],
        threads=1,
        output_format="kst",
        kst_effort="medium",
    )
    kst = out / f"q.{DUMMY_DB_ID}.chromosome.kst"
    assert result.kst_paths == {"chromosome": kst} and not result.kst_merged
    assert result.smoothed_paths == {}
    assert list(open_track(kst)) == [
        ("seqA", 0, 10, "chr1"),
        ("seqA", 10, 14, "novel"),
        ("seqA", 14, 40, "chr1"),
        ("seqB", 0, 20, "chr2"),
    ]
    assert not list(out.glob("*.bed*"))


# --- HKS backend, real karyozip ------------------------------------------

_FEATURES = "featureID\tchromosome\tregion\n1\tchr1\trA\n2\tchr1\trB\n3\tchr2\trC\n"
_HIERARCHY = (
    "feature_set\tchild\tparent\n"
    "chromosome\tautosome\tcategorized\n"
    "chromosome\tchr1\tautosome\n"
    "chromosome\tchr2\tautosome\n"
    "region\tcentromeric\tcategorized\n"
    "region\taSat\tcentromeric\n"
    "region\tHSat\tcentromeric\n"
    "region\trA\taSat\n"
    "region\trB\taSat\n"
    "region\trC\tHSat\n"
)

#: What the stubbed `hks lookup` writes per set (name column filled in).
_LOOKUP = {
    "chromosome": [(0, 10, "chr1"), (10, 14, "novel"), (14, 40, "chr1")],
    "region": [(0, 10, "rA"), (10, 14, "novel"), (14, 30, "rA"), (30, 40, "rB")],
}


@pytest.fixture
def hks_db_root(tmp_path: Path) -> Path:
    """A layout-valid two-set HKS database, with a features.tsv for merged tracks."""
    from karyoscope.installed import InstalledRecord, now_iso, record_install

    root = tmp_path / "hks_root"
    db = root / HKS_DB_ID
    (db / "index").mkdir(parents=True)
    (db / "manifest.yaml").write_text(
        f"id: {HKS_DB_ID}\n"
        'version: "1.0.0"\n'
        'karyoscope_min_version: "2.0.0"\n'
        "index:\n  type: hks\n  basename: index/features\n"
        "hierarchy: hierarchy.tsv\nfeatures: features.tsv\ncolors: colors.tsv\n"
        "kmer:\n  size: 31\n  type: fixed\n  max_size: 31\n"
        "feature_sets:\n  - chromosome\n  - region\n"
        "roles:\n  chromosome_assignment: chromosome\n  region_assignment: region\n"
    )
    (db / "hierarchy.tsv").write_text(_HIERARCHY)
    (db / "features.tsv").write_text(_FEATURES)
    (db / "colors.tsv").write_text("chromosome\tchr1\t#112233\n")
    (db / "index" / "features.hksb").write_bytes(b"hksb")
    for fs in ("chromosome", "region"):
        (db / "index" / f"features.{fs}.hksf").write_bytes(b"hksf")
        (db / "index" / f"features.{fs}.hierarchy.txt").write_text("x\ty\n")
    record_install(
        root,
        HKS_DB_ID,
        InstalledRecord(
            version="1.0.0",
            installed_at=now_iso(),
            source_url="file://test",
            source_sha256="0" * 64,
            directory=HKS_DB_ID,
        ),
    )
    return root


@pytest.fixture
def stub_hks(monkeypatch: pytest.MonkeyPatch) -> list[dict]:
    """Stub `hks lookup` to write `_LOOKUP` (names or ranks); forbid `hks smooth`."""
    calls: list[dict] = []
    monkeypatch.setattr(ann.preflight, "require", lambda *a, **kw: None)

    def fake_lookup(*, feature_set_file, io_pairs, report_query_names, **kw):
        fs = feature_set_file.name.split(".")[1]
        fmt = kw.get("output_format", "bed")
        calls.append({"fs": fs, "names": report_query_names, "n": len(io_pairs), "format": fmt})
        for i, (query, out) in enumerate(io_pairs):
            name = query.stem if report_query_names else str(i)
            rows = "".join(f"{name}\t{s}\t{e}\t{lab}\n" for s, e, lab in _LOOKUP[fs])
            if fmt == "bed":
                out.write_text(rows)
                continue
            # Stand-in for `hks lookup --format kst`: the same rows, as a .kst,
            # without a BED ever landing in the output directory.
            assert out.suffix == ".kst"
            from karyoscope.core.io import karyozip as kmz

            out.parent.mkdir(parents=True, exist_ok=True)
            src = query.parent / f".{out.name}.rows.bed"
            src.write_text(rows)
            kmz.compress(inputs=[src], output=out, effort=kw["kst_effort"])
            src.unlink()
            if kw.get("kst_index"):
                kmz.index(out)

    def no_smooth(**kw):
        raise AssertionError("--format kst must not run hks smooth")

    monkeypatch.setattr(ann, "run_hks_lookup_batch", fake_lookup)
    monkeypatch.setattr(ann, "run_hks_smooth", no_smooth)
    # Step 1 (BED, then karyozip) unless a test opts into hks writing .kst.
    monkeypatch.setattr(ann, "_hks_writes_kst", lambda: False)
    return calls


def _expected(fs: str, name: str) -> str:
    return "".join(f"{name}\t{s}\t{e}\t{lab}\n" for s, e, lab in _LOOKUP[fs])


@pytest.mark.integration
@requires_karyozip
def test_hks_merged_kst_and_no_bed_left(
    hks_db_root: Path, query_fasta: Path, tmp_path: Path, stub_hks: list
) -> None:
    out = tmp_path / "out"
    result = annotate(
        input_path=query_fasta,
        output_dir=out,
        db_root=hks_db_root,
        threads=1,
        output_format="kst",
        smoothed_bed=True,
        bgzip=False,
    )
    merged = out / f"q.{HKS_DB_ID}.kst"
    assert result.kst_merged and set(result.kst_paths.values()) == {merged}
    db_dir = hks_db_root / HKS_DB_ID
    for fs in ("chromosome", "region"):
        with open_track_text(merged, feature_set=fs, db_dir=db_dir) as h:
            assert h.read() == _expected(fs, "q")
    assert _text(result.smoothed_paths["chromosome"]) == "q\t0\t40\tchr1\n"
    assert _text(result.smoothed_paths["region"]) == ("q\t0\t30\trA\nq\t30\t40\trB\n")
    assert not list(out.glob("*.presmoothed*"))


@pytest.mark.integration
@requires_karyozip
@pytest.mark.parametrize("query_names", [True, False], ids=["names", "ranks"])
def test_hks_per_set_kst_round_trips_names_and_ranks(
    hks_db_root: Path, query_fasta: Path, tmp_path: Path, stub_hks: list, query_names: bool
) -> None:
    out = tmp_path / "out"
    result = annotate(
        input_path=query_fasta,
        output_dir=out,
        db_root=hks_db_root,
        feature_sets=["region"],
        threads=1,
        output_format="kst",
        query_names=query_names,
    )
    kst = out / f"q.{HKS_DB_ID}.region.kst"
    assert result.kst_paths == {"region": kst}
    with open_track_text(kst) as h:
        assert h.read() == _expected("region", "q" if query_names else "0")
    assert stub_hks[-1]["names"] is query_names
    assert not list(out.glob("*.bed"))


@pytest.mark.integration
@requires_karyozip
def test_hks_batch_writes_one_merged_kst_per_input(
    hks_db_root: Path, tmp_path: Path, stub_hks: list
) -> None:
    inputs = []
    for stem in ("a", "b"):
        fa = tmp_path / f"{stem}.fa"
        fa.write_text(">seqA\n" + "A" * 40 + "\n")
        inputs.append(fa)
    out = tmp_path / "out"
    results = annotate_batch(
        input_paths=inputs,
        output_dir=out,
        db_root=hks_db_root,
        threads=1,
        output_format="kst",
    )
    assert [c["n"] for c in stub_hks] == [2, 2]  # one lookup per set for both inputs
    for p in inputs:
        r = results[p]
        merged = out / f"{p.stem}.{HKS_DB_ID}.kst"
        assert r.kst_merged and r.kst_paths["region"] == merged
        rows = list(open_track(merged, feature_set="region", db_dir=hks_db_root / HKS_DB_ID))
        assert [row[0] for row in rows] == [p.stem] * 4
        assert r.manifest_path == out / f"{p.stem}.{HKS_DB_ID}.annotate_manifest.json"
    assert not list(out.glob("*.bed"))


@pytest.mark.integration
@requires_karyozip
def test_failed_compress_leaves_no_temp_beds(
    hks_db_root: Path, query_fasta: Path, tmp_path: Path, stub_hks: list, monkeypatch
) -> None:
    def boom(**kw):
        raise KaryoscopeError("compress failed")

    monkeypatch.setattr(ann.kmz, "compress", boom)
    out = tmp_path / "out"
    with pytest.raises(KaryoscopeError, match="compress failed"):
        annotate(
            input_path=query_fasta,
            output_dir=out,
            db_root=hks_db_root,
            threads=1,
            output_format="kst",
        )
    assert not list(out.glob("*.bed"))


# --- downstream readers take a .kst --------------------------------------


@pytest.fixture
def kst_and_beds(
    populated_db_root: Path, query_fasta: Path, tmp_path: Path, stub_kmc: dict
) -> tuple[Path, dict[str, Path], Path]:
    """A merged .kst and the BED run it must equal, from the same combined BED."""
    beds = annotate(
        input_path=query_fasta,
        output_dir=tmp_path / "bed",
        db_root=populated_db_root,
        threads=1,
        bgzip=False,
    )
    kst = annotate(
        input_path=query_fasta,
        output_dir=tmp_path / "kst",
        db_root=populated_db_root,
        threads=1,
        output_format="kst",
    )
    return kst.kst_paths["chromosome"], beds.smoothed_paths, populated_db_root / DUMMY_DB_ID


@pytest.mark.integration
@requires_karyozip
@pytest.mark.parametrize("threads", [1, 2])
def test_bin_features_reads_a_kst(kst_and_beds, tmp_path: Path, threads: int) -> None:
    from karyoscope.core.bin import bin_features

    kst, smoothed_beds, db_dir = kst_and_beds
    from_bed, from_kst = tmp_path / "a.bed", tmp_path / "b.bed"
    bin_features(smoothed_beds["region"], from_bed, bin_size=8, threads=threads, chunk_size=1)
    src = Track(kst, smoothed=True, feature_set="region", db_dir=db_dir)
    bin_features(src, from_kst, bin_size=8, threads=threads, chunk_size=1)
    assert from_kst.read_text() == from_bed.read_text()


@pytest.mark.integration
@requires_karyozip
def test_bin_cli_reads_a_kst(kst_and_beds, tmp_path: Path, populated_db_root: Path) -> None:
    kst, smoothed_beds, _ = kst_and_beds
    runner = CliRunner()
    want, got = tmp_path / "want.bed", tmp_path / "got.bed"
    base = ["bin", "-b", "8", "--feature-set", "chromosome", "--db", DUMMY_DB_ID]
    base += ["--db-root", str(populated_db_root)]
    r1 = runner.invoke(main, [*base, "-i", str(smoothed_beds["chromosome"]), "-o", str(want)])
    r2 = runner.invoke(main, [*base, "-i", str(kst), "-o", str(got)])
    assert r1.exit_code == 0, r1.output
    assert r2.exit_code == 0, r2.output
    assert got.read_text() == want.read_text()
    # Smoothing a .kst needs the database; without --feature-set there is none.
    r3 = runner.invoke(main, ["bin", "-b", "8", "-i", str(kst), "-o", str(got)])
    assert r3.exit_code != 0 and "needs the database" in r3.output
    r4 = runner.invoke(
        main, ["bin", "-b", "8", "-i", str(want), "-o", str(got), "--kst-track", "smoothed"]
    )
    assert r4.exit_code != 0 and "only applies to a .kst" in r4.output


@pytest.mark.integration
@requires_karyozip
def test_remap_bed_reads_a_kst(kst_and_beds, tmp_path: Path, populated_db_root: Path) -> None:
    from karyoscope.core.io.scaffold_map import MapRow, write_map

    kst, smoothed_beds, _ = kst_and_beds
    mp = tmp_path / "q.KS.scaffold_map.tsv"
    write_map(
        [
            MapRow(
                new_name="chr1_hap1_seqA",
                original_name="seqA",
                input_file="q.fasta",
                hap="hap1",
                chromosome="chr1",
                flipped=True,
                length=40,
                stats="TPCQT",
            )
        ],
        mp,
    )
    runner = CliRunner()
    want, got = tmp_path / "want.bed", tmp_path / "got.bed"
    r1 = runner.invoke(
        main, ["remap-bed", "-b", str(smoothed_beds["region"]), "-m", str(mp), "-o", str(want)]
    )
    r2 = runner.invoke(
        main,
        [
            "remap-bed", "-b", str(kst), "-m", str(mp), "-o", str(got), "--feature-set", "region",
            "--db", DUMMY_DB_ID, "--db-root", str(populated_db_root),
        ],
    )  # fmt: skip
    assert r1.exit_code == 0, r1.output
    assert r2.exit_code == 0, r2.output
    assert got.read_text() == want.read_text() != ""
    r3 = runner.invoke(main, ["remap-bed", "-b", str(kst), "-m", str(mp), "-o", str(got),
                              "--db-root", str(populated_db_root)])  # fmt: skip
    assert r3.exit_code != 0 and "--feature-set" in r3.output


@pytest.mark.integration
@requires_karyozip
def test_scaffold_cascade_finds_and_reads_a_kst(kst_and_beds, populated_db_root: Path) -> None:
    from karyoscope.core.scaffold_run import _annotation_bed_path, _load_label_mass

    kst, smoothed_beds, db_dir = kst_and_beds
    found = _annotation_bed_path(kst.parent, "q", DUMMY_DB_ID, "chromosome", db_dir=db_dir)
    assert found == Track(kst, smoothed=True, feature_set="chromosome", db_dir=db_dir)
    assert found.is_file()
    assert _load_label_mass(found) == _load_label_mass(smoothed_beds["chromosome"])


@pytest.mark.integration
@requires_karyozip
def test_scaffold_run_on_kst_matches_bed(populated_db_root: Path, tmp_path: Path) -> None:
    """``scaffold --mode bed`` from per-set ``.kst`` writes what it writes from BEDs.

    The binned BEDs are derived (auto) in both runs, so the ``.kst`` run bins
    its smoothed track through karyozip, tallies label mass from it, and
    rewrites it into scaffold coordinates.
    """
    from karyoscope.core.io import karyozip as kmz
    from karyoscope.core.scaffold_run import InputSpec, scaffold_run

    mb = 1_000_000
    tracks = {
        "chromosome": f"ctgA\t0\t{6 * mb}\tchr1\nctgB\t0\t{6 * mb}\tchr2\n",
        "region": f"ctgA\t0\t{3 * mb}\trA\nctgA\t{3 * mb}\t{6 * mb}\trB\nctgB\t0\t{6 * mb}\trC\n",
    }
    results = {}
    for fmt in ("bed", "kst"):
        d = tmp_path / fmt
        d.mkdir()
        fasta = d / "samp.fasta"
        fasta.write_text(">ctgA\nACGT\n>ctgB\nACGT\n")
        (d / "samp.telo").write_text("")
        for fs, rows in tracks.items():
            bed = d / f"samp.{DUMMY_DB_ID}.{fs}.smoothed.bed"
            bed.write_text(rows)
            if fmt == "kst":
                kmz.compress(inputs=[bed], output=d / f"samp.{DUMMY_DB_ID}.{fs}.kst")
                bed.unlink()
        out = scaffold_run(
            inputs=[InputSpec(name=None, path=fasta, telo_path=None)],
            db_root=populated_db_root,
            mode="bed",
            bgzip=False,
            auto=True,
            output_dir=d,
        )
        r = out["samp.fasta"]
        results[fmt] = (
            r.map_path.read_text(),
            {fs: p.read_text() for fs, p in r.scaffolded_beds.items()},
        )
    assert results["kst"] == results["bed"]
    assert set(results["kst"][1]) == {"chromosome", "region"}


# --- hks lookup --format kst (step 2) -------------------------------------


@pytest.mark.integration
@requires_karyozip
@pytest.mark.parametrize("subset", [None, ["region"]], ids=["merged", "per-set"])
def test_hks_writing_kst_matches_bed_then_compress(
    hks_db_root: Path,
    query_fasta: Path,
    tmp_path: Path,
    stub_hks: list,
    monkeypatch: pytest.MonkeyPatch,
    subset: list[str] | None,
) -> None:
    """With an hks that writes .kst, no BED is written and the bytes are unchanged."""
    outs = {}
    for via_hks in (False, True):
        monkeypatch.setattr(ann, "_hks_writes_kst", lambda v=via_hks: v)
        out = tmp_path / ("hks" if via_hks else "bed")
        result = annotate(
            input_path=query_fasta,
            output_dir=out,
            db_root=hks_db_root,
            feature_sets=subset,
            threads=1,
            output_format="kst",
            smoothed_bed=True,
            bgzip=False,
        )
        outs[via_hks] = (out, result)
    assert {c["format"] for c in stub_hks[: len(stub_hks) // 2]} == {"bed"}
    assert {c["format"] for c in stub_hks[len(stub_hks) // 2 :]} == {"kst"}
    (bed_dir, _), (hks_dir, hks_res) = outs[False], outs[True]
    assert hks_res.kst_merged is (subset is None)
    names = sorted(p.name for p in hks_dir.iterdir())
    assert names == sorted(p.name for p in bed_dir.iterdir())  # no temp left behind
    for name in names:
        if name.endswith(".json"):
            continue
        assert (hks_dir / name).read_bytes() == (bed_dir / name).read_bytes(), name
    doc = json.loads(hks_res.manifest_path.read_text())
    assert doc["kst"]["written_by"] == "hks lookup --format kst"


def test_direct_kst_disk_estimate_has_no_presmoothed_bed() -> None:
    bases = 1_000_000_000
    est = ann.estimate_output_bytes
    common = {"input_bases": bases, "n_feature_sets": 6, "keep_presmoothed": True, "smooth": True}
    direct = est(**common, index_type="hks", output_format="kst", kst_merged=True, kst_via_hks=True)
    assert direct == int(bases * ann.KST_BYTES_PER_BASE * 12)
    assert direct < est(**common, index_type="hks", output_format="kst", kst_merged=True) / 5


def test_hks_writes_kst_needs_hks_support_and_karyozip_0_3(monkeypatch) -> None:
    monkeypatch.setattr(ann, "hks_supports_kst", lambda: False)
    assert not ann._hks_writes_kst()
    monkeypatch.setattr(ann, "hks_supports_kst", lambda: True)
    for line, want in (
        ("karyozip 0.3.0, format KST00001", True),
        ("karyozip 0.2.0, format KST00001", False),
        (None, False),
    ):
        monkeypatch.setattr(ann.kmz, "karyozip_version", lambda line=line: line)
        assert ann._hks_writes_kst() is want


# --- centromeres and karyotype from a .kst, end to end ---------------------


def _get_featureids_available() -> bool:
    from karyoscope.core.external import ToolNotFoundError
    from karyoscope.core.io.kmc import get_featureids_binary

    try:
        get_featureids_binary()
    except ToolNotFoundError:
        return False
    return True


@pytest.mark.integration
@requires_karyozip
@requires_bgzip
@pytest.mark.skipif(shutil.which("seqtk") is None, reason="needs seqtk on PATH")
@pytest.mark.skipif(not _get_featureids_available(), reason="needs get_featureIDs")
@pytest.mark.parametrize("subset", [None, ["chromosome", "region"]], ids=["merged", "per-set"])
def test_centromeres_and_karyotype_from_kst_match_bed(
    populated_db_root: Path, tmp_path: Path, subset: list[str] | None
) -> None:
    """The real KMC annotate on the dummy db, as BED and as .kst, then the cascade.

    ``centromeres`` and ``karyotype`` (genome and centromere modes) run on each
    directory; every file they write must be byte-identical (``.gz`` compared
    decompressed). The SVGs embed no path or time, so nothing is normalised.
    """
    seed = "ACGTGCTAGCTAGGCTATCGTAC"
    runner = CliRunner()
    produced: dict[str, dict[str, bytes]] = {}
    for fmt in ("bed", "kst"):
        d = tmp_path / fmt
        d.mkdir()
        fa = d / "assembly.fa"
        fa.write_text(f">seq_for_chr1\n{seed[:21]}\n>seq_for_chr2\n{seed[2:23]}\n")
        annotate(
            input_path=fa,
            output_dir=d,
            db_root=populated_db_root,
            feature_sets=subset,
            threads=1,
            output_format=fmt,
        )
        before = {p.name for p in d.iterdir()}
        sets = [a for fs in (subset or []) for a in ("--feature-set", fs)]
        for args in (
            ["centromeres", "--coarse-bin-size", "10", "--fine-bin-size", "0"],
            ["karyotype", "--mode", "genome", "--bin-size", "10", *sets],
            ["karyotype", "--mode", "centromere", "--bin-size", "10", *sets],
        ):
            result = runner.invoke(
                main,
                [*args, "-i", f"hap1={fa}", "--outdir", str(d), "--min-scaffold-length", "1"],
            )
            assert result.exit_code == 0, (fmt, args, result.output)
        produced[fmt] = {
            p.name: (gzip.decompress(p.read_bytes()) if p.suffix == ".gz" else p.read_bytes())
            for p in d.iterdir()
            if p.name not in before
        }
    assert any(n.endswith(".svg") for n in produced["kst"])
    assert any("centromeres" in n for n in produced["kst"])
    assert produced["kst"] == produced["bed"]
