"""Orchestration tests for :func:`karyoscope.core.annotate.annotate`.

The audit found the three big orchestrators covered almost entirely by
integration tests, with only their pure helpers unit-tested. These tests
run ``annotate()`` / ``annotate_batch()`` end to end with ONLY the
external binaries stubbed (``run_get_featureids``, ``bgzip_file``, the
HKS backend, the tool preflight) — database resolution, manifest and
hierarchy validation, feature translation, the real smoothing pass, path
naming, intermediate reuse, and result assembly all execute for real.
No external tool is needed, so these run in the default unit pass.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

import karyoscope.core.annotate as ann
from karyoscope.core.annotate import AnnotateResult, annotate, annotate_batch
from karyoscope.core.io.kmc import combined_bed_path, combined_marker_path, write_combined_marker
from karyoscope.exceptions import DatabaseNotFoundError, KaryoscopeError

DUMMY_DB_ID = "KS_dummy_test_v1"

#: Combined-BED content for the dummy db (features.tsv: 1 -> chr1/rA,
#: 2 -> chr1/rB, 3 -> chr2/rC; id 0 is the novel sentinel).
_COMBINED_LINES = (
    "seqA\t0\t10\t1\nseqA\t10\t14\t0\nseqA\t14\t30\t1\nseqA\t30\t40\t2\nseqB\t0\t20\t3\n"
)


@pytest.fixture
def query_fasta(tmp_path: Path) -> Path:
    fa = tmp_path / "q.fasta"
    fa.write_text(">seqA\n" + "A" * 40 + "\n>seqB\n" + "C" * 20 + "\n")
    return fa


@pytest.fixture
def stub_externals(monkeypatch: pytest.MonkeyPatch) -> dict:
    """Stub every external binary; record what the orchestration asked for.

    ``bgzip_file`` renames in place like the real one; the fake
    ``run_get_featureids`` writes a real combined BED so everything
    downstream of the binary runs unstubbed.
    """
    calls: dict = {"get_featureids": 0, "bgzipped": []}

    monkeypatch.setattr(ann.preflight, "require", lambda *a, **kw: None)

    def fake_bgzip(path: Path, threads: int = 1) -> Path:
        out = Path(str(path) + ".gz")
        shutil.move(path, out)
        calls["bgzipped"].append(out.name)
        return out

    monkeypatch.setattr(ann, "bgzip_file", fake_bgzip)

    def fake_get_featureids(*, db_path, input_path, output_dir, threads, prefix, **kw):
        calls["get_featureids"] += 1
        bed = combined_bed_path(output_dir, prefix)
        bed.write_text(_COMBINED_LINES)
        write_combined_marker(bed, prefix=prefix, db_path=db_path, input_path=input_path)
        return bed

    monkeypatch.setattr(ann, "run_get_featureids", fake_get_featureids)
    return calls


def _read_bed(path: Path) -> list[tuple[str, int, int, str]]:
    out = []
    for line in path.read_text().splitlines():
        seq, start, end, name = line.split("\t")
        out.append((seq, int(start), int(end), name))
    return out


# --- KMC path ---------------------------------------------------------


def test_kmc_full_run_writes_both_outputs_per_set(
    populated_db_root: Path, query_fasta: Path, tmp_path: Path, stub_externals: dict
) -> None:
    """The default flag set produces presmoothed + smoothed bgzipped BEDs
    for every manifest feature set, named by input and database, and
    removes the combined intermediate and its marker."""
    out = tmp_path / "out"
    result = annotate(
        input_path=query_fasta,
        output_dir=out,
        db_root=populated_db_root,
        threads=1,
    )

    assert isinstance(result, AnnotateResult)
    for fs in ("chromosome", "region"):
        pre = out / f"q.{DUMMY_DB_ID}.{fs}.presmoothed.bed.gz"
        smo = out / f"q.{DUMMY_DB_ID}.{fs}.smoothed.bed.gz"
        assert result.presmoothed_paths[fs] == pre and pre.is_file()
        assert result.smoothed_paths[fs] == smo and smo.is_file()
    assert sorted(result.all_output_paths) == sorted(
        list(result.presmoothed_paths.values()) + list(result.smoothed_paths.values())
    )
    assert stub_externals["get_featureids"] == 1
    assert len(stub_externals["bgzipped"]) == 4

    combined = combined_bed_path(out, f"q.{DUMMY_DB_ID}")
    assert result.combined_intermediate is None
    assert not combined.exists()
    assert not combined_marker_path(combined).exists()


def test_kmc_no_smooth_split_translates_ids(
    populated_db_root: Path, query_fasta: Path, tmp_path: Path, stub_externals: dict
) -> None:
    """--no-smooth takes the in-process split path: presmoothed only,
    ids translated per feature set, adjacent same-name runs merged."""
    out = tmp_path / "out"
    result = annotate(
        input_path=query_fasta,
        output_dir=out,
        db_root=populated_db_root,
        threads=1,
        smooth=False,
        bgzip=False,
    )

    assert result.smoothed_paths == {}
    chrom = _read_bed(result.presmoothed_paths["chromosome"])
    # ids 1 and 2 are both chr1, so 14-30 and 30-40 merge; the novel gap
    # at 10-14 keeps 0-10 separate.
    assert chrom == [
        ("seqA", 0, 10, "chr1"),
        ("seqA", 10, 14, "novel"),
        ("seqA", 14, 40, "chr1"),
        ("seqB", 0, 20, "chr2"),
    ]
    region = _read_bed(result.presmoothed_paths["region"])
    assert region == [
        ("seqA", 0, 10, "rA"),
        ("seqA", 10, 14, "novel"),
        ("seqA", 14, 30, "rA"),
        ("seqA", 30, 40, "rB"),
        ("seqB", 0, 20, "rC"),
    ]
    assert stub_externals["bgzipped"] == []


def test_kmc_reuses_verified_combined_bed(
    populated_db_root: Path,
    query_fasta: Path,
    tmp_path: Path,
    stub_externals: dict,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A complete combined BED with a matching marker skips the k-mer
    query; --force regenerates it; a marker-less file is not trusted."""
    out = tmp_path / "out"
    out.mkdir()
    prefix = f"q.{DUMMY_DB_ID}"
    combined = combined_bed_path(out, prefix)
    combined.write_text(_COMBINED_LINES)
    write_combined_marker(combined, prefix=prefix, db_path=Path("/db"), input_path=query_fasta)

    common = dict(
        input_path=query_fasta,
        output_dir=out,
        db_root=populated_db_root,
        threads=1,
        smooth=False,
        bgzip=False,
        keep_intermediates=True,
    )
    annotate(**common)
    assert stub_externals["get_featureids"] == 0  # reused

    annotate(**common, force=True)
    assert stub_externals["get_featureids"] == 1  # regenerated

    combined.write_text(_COMBINED_LINES + "seqB\t20\t25\t1\n")  # marker now stale
    annotate(**common)
    assert stub_externals["get_featureids"] == 2  # not trusted


def test_kmc_keep_intermediates_reports_the_kept_path(
    populated_db_root: Path, query_fasta: Path, tmp_path: Path, stub_externals: dict
) -> None:
    out = tmp_path / "out"
    result = annotate(
        input_path=query_fasta,
        output_dir=out,
        db_root=populated_db_root,
        threads=1,
        smooth=False,
        bgzip=False,
        keep_intermediates=True,
    )
    assert result.combined_intermediate == combined_bed_path(out, f"q.{DUMMY_DB_ID}")
    assert result.combined_intermediate.is_file()


def test_feature_set_subset_and_unknown(
    populated_db_root: Path, query_fasta: Path, tmp_path: Path, stub_externals: dict
) -> None:
    out = tmp_path / "out"
    result = annotate(
        input_path=query_fasta,
        output_dir=out,
        db_root=populated_db_root,
        threads=1,
        smooth=False,
        bgzip=False,
        feature_sets=["region"],
    )
    assert list(result.presmoothed_paths) == ["region"]

    with pytest.raises(KaryoscopeError, match="not declared"):
        annotate(
            input_path=query_fasta,
            output_dir=out,
            db_root=populated_db_root,
            threads=1,
            feature_sets=["nope"],
        )


def test_rejections_before_any_work(
    populated_db_root: Path, query_fasta: Path, tmp_path: Path, stub_externals: dict
) -> None:
    """The no-output combination and a missing input fail before the
    database is even resolved / the query step is reached."""
    with pytest.raises(KaryoscopeError, match="no output would be produced"):
        annotate(
            input_path=query_fasta,
            output_dir=tmp_path,
            db_root=tmp_path / "nonexistent-root",
            smooth=False,
            keep_presmoothed=False,
        )
    with pytest.raises(KaryoscopeError, match="input file not found"):
        annotate(
            input_path=tmp_path / "missing.fasta",
            output_dir=tmp_path,
            db_root=populated_db_root,
        )
    with pytest.raises(DatabaseNotFoundError):
        annotate(
            input_path=query_fasta,
            output_dir=tmp_path,
            db_root=tmp_path / "empty-root",
        )
    assert stub_externals["get_featureids"] == 0


# --- HKS path ---------------------------------------------------------


HKS_DB_ID = "HKS_fake_test_v1"

_HKS_HIERARCHY = (
    "feature_set\tchild\tparent\n"
    "chromosome\tautosome\tcategorized\n"
    "chromosome\tchr1\tautosome\n"
    "chromosome\tchr2\tautosome\n"
)


@pytest.fixture
def hks_db_root(tmp_path: Path) -> Path:
    """A db root holding a minimal but layout-valid HKS-backend database.

    The index files are stand-ins (the backend is stubbed in these
    tests); everything the manifest references exists, so
    ``validate_database_layout`` passes for real.
    """
    from karyoscope.installed import InstalledRecord, now_iso, record_install

    root = tmp_path / "hks_root"
    db = root / HKS_DB_ID
    (db / "index").mkdir(parents=True)
    (db / "manifest.yaml").write_text(
        f"id: {HKS_DB_ID}\n"
        'version: "1.0.0"\n'
        'karyoscope_min_version: "2.0.0"\n'
        "index:\n"
        "  type: hks\n"
        "  basename: index/features\n"
        "hierarchy: hierarchy.tsv\n"
        "colors: colors.tsv\n"
        "kmer:\n"
        "  size: 31\n"
        "  type: variable\n"
        "  max_size: 31\n"
        "feature_sets:\n"
        "  - chromosome\n"
        "roles:\n"
        "  chromosome_assignment: chromosome\n"
    )
    (db / "hierarchy.tsv").write_text(_HKS_HIERARCHY)
    (db / "colors.tsv").write_text("chromosome\tchr1\t#112233\n")
    (db / "index" / "features.hksb").write_bytes(b"hksb")
    (db / "index" / "features.chromosome.hksf").write_bytes(b"hksf")
    (db / "index" / "features.chromosome.hierarchy.txt").write_text("chr1\tautosome\n")
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
def stub_hks_backend(monkeypatch: pytest.MonkeyPatch) -> list[dict]:
    """Replace the HKS backend with a recorder that writes what it is told."""
    calls: list[dict] = []

    def fake_backend(**kwargs):
        calls.append(kwargs)
        for p, by_fs in kwargs["presmoothed_by_input"].items():
            for path in by_fs.values():
                path.write_text(f"{p.name}\t0\t1\tchr1\n")
        for p, by_fs in kwargs["smoothed_by_input"].items():
            for path in by_fs.values():
                path.write_text(f"{p.name}\t0\t1\tchr1\n")

    monkeypatch.setattr(ann, "_run_hks_backend", fake_backend)
    monkeypatch.setattr(ann.preflight, "require", lambda *a, **kw: None)
    return calls


def test_hks_dispatch_paths_and_result(
    hks_db_root: Path, query_fasta: Path, tmp_path: Path, stub_hks_backend: list
) -> None:
    """The HKS branch hands the backend db-derived paths and produces a
    result with no combined intermediate."""
    out = tmp_path / "out"
    result = annotate(
        input_path=query_fasta,
        output_dir=out,
        db_root=hks_db_root,
        threads=1,
        bgzip=False,
    )
    assert len(stub_hks_backend) == 1
    call = stub_hks_backend[0]
    assert call["input_paths"] == [query_fasta]
    assert call["prefixes"] == {query_fasta: f"q.{HKS_DB_ID}"}
    assert call["requested"] == ["chromosome"]

    assert result.combined_intermediate is None
    assert (
        result.presmoothed_paths["chromosome"].name == f"q.{HKS_DB_ID}.chromosome.presmoothed.bed"
    )
    assert result.smoothed_paths["chromosome"].is_file()


def test_hks_k_override_tags_outputs_and_validates(
    hks_db_root: Path, query_fasta: Path, tmp_path: Path, stub_hks_backend: list
) -> None:
    """On a variable-k index an explicit k tags the output names; a k
    beyond max_size is refused."""
    out = tmp_path / "out"
    result = annotate(
        input_path=query_fasta,
        output_dir=out,
        db_root=hks_db_root,
        threads=1,
        bgzip=False,
        k=15,
    )
    assert result.presmoothed_paths["chromosome"].name == (
        f"q.{HKS_DB_ID}.k15.chromosome.presmoothed.bed"
    )
    assert stub_hks_backend[-1]["k"] == 15

    with pytest.raises(KaryoscopeError, match="max_size"):
        annotate(
            input_path=query_fasta,
            output_dir=out,
            db_root=hks_db_root,
            threads=1,
            k=40,
        )


def test_batch_single_input_delegates_to_annotate(
    hks_db_root: Path, query_fasta: Path, tmp_path: Path, stub_hks_backend: list
) -> None:
    result = annotate_batch(
        input_paths=[query_fasta],
        output_dir=tmp_path / "out",
        db_root=hks_db_root,
        threads=1,
        bgzip=False,
    )
    assert set(result) == {query_fasta}
    assert isinstance(result[query_fasta], AnnotateResult)
    assert len(stub_hks_backend) == 1


def test_batch_multi_input_shares_one_backend_invocation(
    hks_db_root: Path, tmp_path: Path, stub_hks_backend: list
) -> None:
    """Several inputs go to the backend together (one index load per
    feature set for the cohort), with per-input prefixes and results."""
    a = tmp_path / "a.fasta"
    b = tmp_path / "b.fasta"
    a.write_text(">s\nACGT\n")
    b.write_text(">s\nACGT\n")

    results = annotate_batch(
        input_paths=[a, b],
        output_dir=tmp_path / "out",
        db_root=hks_db_root,
        threads=1,
        bgzip=False,
    )
    assert len(stub_hks_backend) == 1
    call = stub_hks_backend[0]
    assert call["input_paths"] == [a, b]
    assert call["prefixes"] == {a: f"a.{HKS_DB_ID}", b: f"b.{HKS_DB_ID}"}

    assert set(results) == {a, b}
    for p in (a, b):
        assert results[p].presmoothed_paths["chromosome"].name.startswith(f"{p.stem}.{HKS_DB_ID}.")
        assert results[p].presmoothed_paths["chromosome"].is_file()
