"""Orchestration tests for :func:`karyoscope.core.scaffold_run.scaffold_run`.

Companion to ``test_annotate_orchestration.py``: the pipeline runs for
real from prepared on-disk intermediates (``auto=False``), so database
resolution, role resolution, contig classification and orientation, map
and stats writing, and the scaffolded-BED rewrite all execute unstubbed
— no external tool is involved in ``mode="bed"`` with ``bgzip=False``.
The auto-derive cascade is tested separately with the three derive
layers stubbed, pinning *when* each is invoked.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

import karyoscope.core.scaffold_run as sr
from karyoscope.core.io.scaffold_map import read_map
from karyoscope.core.scaffold_run import InputSpec, ScaffoldResult, scaffold_run
from karyoscope.exceptions import ScaffoldError

DUMMY_DB_ID = "KS_dummy_test_v1"
MB = 1_000_000


@pytest.fixture
def prepared_input(tmp_path: Path) -> Path:
    """A FASTA plus every intermediate the ``auto=False`` path requires.

    Two 6 Mb contigs (over the 5 Mb no-telomere floor): ``ctgA`` binned
    entirely to chr1 and ``ctgB`` to chr2, with matching full-resolution
    smoothed BEDs for both dummy-db feature sets, an empty telo file
    (no telomeres), and per-role binned BEDs at the default 1 Mb.
    """
    fasta = tmp_path / "samp.fasta"
    fasta.write_text(">ctgA\nACGT\n>ctgB\nACGT\n")

    def bed(name: str, text: str) -> None:
        (tmp_path / f"samp.{DUMMY_DB_ID}.{name}").write_text(text)

    def binned_gz(fs: str, rows: str) -> None:
        with gzip.open(tmp_path / f"samp.{DUMMY_DB_ID}.{fs}.smoothed.binned{MB}.bed.gz", "wt") as h:
            h.write(rows)

    chrom_bins = "".join(f"ctgA\t{i * MB}\t{(i + 1) * MB}\tchr1\n" for i in range(6)) + "".join(
        f"ctgB\t{i * MB}\t{(i + 1) * MB}\tchr2\n" for i in range(6)
    )
    region_bins = "".join(f"ctgA\t{i * MB}\t{(i + 1) * MB}\trA\n" for i in range(6)) + "".join(
        f"ctgB\t{i * MB}\t{(i + 1) * MB}\trC\n" for i in range(6)
    )
    binned_gz("chromosome", chrom_bins)
    binned_gz("region", region_bins)

    bed(
        "chromosome.smoothed.bed",
        f"ctgA\t0\t{6 * MB}\tchr1\nctgB\t0\t{6 * MB}\tchr2\n",
    )
    bed(
        "region.smoothed.bed",
        f"ctgA\t0\t{3 * MB}\trA\nctgA\t{3 * MB}\t{6 * MB}\trB\nctgB\t0\t{6 * MB}\trC\n",
    )
    (tmp_path / "samp.telo").write_text("")
    return fasta


def _run(fasta: Path, db_root: Path, **overrides) -> dict[str, ScaffoldResult]:
    kwargs = dict(
        inputs=[InputSpec(name=None, path=fasta, telo_path=None)],
        db_root=db_root,
        mode="bed",
        bgzip=False,
        auto=False,
        output_dir=fasta.parent,
    )
    kwargs.update(overrides)
    return scaffold_run(**kwargs)


# --- input validation guards (no database work needed for the first two)


def test_guards_reject_bad_invocations(populated_db_root: Path, tmp_path: Path) -> None:
    with pytest.raises(ScaffoldError, match="at least one"):
        scaffold_run(inputs=[], db_root=populated_db_root)
    fa = tmp_path / "x.fasta"
    fa.write_text(">c\nA\n")
    with pytest.raises(ScaffoldError, match="unknown mode"):
        _run(fa, populated_db_root, mode="nope")
    with pytest.raises(ScaffoldError, match="combine_chromosomes requires a FASTA"):
        _run(fa, populated_db_root, combine_chromosomes=True)
    with pytest.raises(ScaffoldError, match="not declared in manifest"):
        _run(fa, populated_db_root, feature_sets=["nope"])


# --- auto=False: every missing prerequisite fails loudly ---------------


def test_auto_off_reports_each_missing_prerequisite(
    populated_db_root: Path, prepared_input: Path
) -> None:
    out_dir = prepared_input.parent

    telo = out_dir / "samp.telo"
    telo.unlink()
    with pytest.raises(ScaffoldError, match="auto-derive"):
        _run(prepared_input, populated_db_root)
    telo.write_text("")

    binned = out_dir / f"samp.{DUMMY_DB_ID}.chromosome.smoothed.binned{MB}.bed.gz"
    moved = binned.with_suffix(".bak")
    binned.rename(moved)
    with pytest.raises(ScaffoldError, match="auto-derive"):
        _run(prepared_input, populated_db_root)
    moved.rename(binned)

    ann = out_dir / f"samp.{DUMMY_DB_ID}.chromosome.smoothed.bed"
    ann.unlink()
    with pytest.raises(ScaffoldError, match="auto-derive"):
        _run(prepared_input, populated_db_root)


# --- the real pipeline from prepared intermediates ---------------------


def test_bed_mode_end_to_end_from_prepared_intermediates(
    populated_db_root: Path, prepared_input: Path
) -> None:
    """classify + orient, map/stats writing and the scaffolded-BED
    rewrite run for real; the outputs agree with each other."""
    out_dir = prepared_input.parent
    results = _run(prepared_input, populated_db_root)

    assert set(results) == {"samp.fasta"}
    res = results["samp.fasta"]
    assert res.combined is False
    assert res.scaffolded_fasta is None
    assert res.map_path == out_dir / f"samp.{DUMMY_DB_ID}.scaffold_map.tsv"
    assert res.stats_path == out_dir / f"samp.{DUMMY_DB_ID}.scaffold_stats.tsv"
    assert res.map_path.is_file() and res.stats_path.is_file()

    rows = read_map(res.map_path)
    assert {r.original_name for r in rows} == {"ctgA", "ctgB"}
    by_orig = {r.original_name: r for r in rows}
    assert by_orig["ctgA"].chromosome == "chr1"
    assert by_orig["ctgB"].chromosome == "chr2"
    assert all(r.length == 6 * MB for r in rows)

    # The rewritten BEDs use the map's encoded names, with every record
    # from the source BED accounted for.
    assert set(res.scaffolded_beds) == {"chromosome", "region"}
    for fs, n_src_rows in (("chromosome", 2), ("region", 3)):
        out = res.scaffolded_beds[fs]
        assert out == out_dir / f"samp.{DUMMY_DB_ID}.{fs}.smoothed.scaffolded.bed"
        lines = [ln.split("\t") for ln in out.read_text().splitlines()]
        assert len(lines) == n_src_rows
        assert {ln[0] for ln in lines} == {r.new_name for r in rows}


def test_combined_mode_reuses_the_discovery_pass_lengths(
    populated_db_root: Path, prepared_input: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """combine_chromosomes derives contig names and true lengths from a
    single read of the input FASTA: ``read_fasta_lengths`` runs exactly
    once and ``read_fasta_contig_names`` is never called."""
    calls = {"lengths": 0}
    real_lengths = sr.read_fasta_lengths

    def counting_lengths(path: Path):
        calls["lengths"] += 1
        return real_lengths(path)

    def forbidden_names(path: Path):
        raise AssertionError("combined mode must not re-scan the FASTA for names")

    monkeypatch.setattr(sr, "read_fasta_lengths", counting_lengths)
    monkeypatch.setattr(sr, "read_fasta_contig_names", forbidden_names)

    results = _run(prepared_input, populated_db_root, mode="fasta", combine_chromosomes=True)
    res = results["samp.fasta"]
    assert res.combined is True
    assert calls["lengths"] == 1

    assert res.scaffolded_fasta is not None and res.scaffolded_fasta.is_file()
    heads = sorted(
        ln[1:].split()[0]
        for ln in res.scaffolded_fasta.read_text().splitlines()
        if ln.startswith(">")
    )
    assert len(heads) == 2
    assert heads[0].startswith("chr1") and heads[1].startswith("chr2")


def test_bed_mode_can_skip_the_bed_rewrite(populated_db_root: Path, prepared_input: Path) -> None:
    """write_scaffolded_beds=False still writes the map but no BEDs
    (the karyotype cascade's arrangement)."""
    results = _run(prepared_input, populated_db_root, write_scaffolded_beds=False)
    res = results["samp.fasta"]
    assert res.map_path.is_file()
    assert res.scaffolded_beds == {}


def test_feature_set_restriction_limits_the_rewrite(
    populated_db_root: Path, prepared_input: Path
) -> None:
    results = _run(prepared_input, populated_db_root, feature_sets=["region"])
    res = results["samp.fasta"]
    assert set(res.scaffolded_beds) == {"region"}


# --- the auto-derive cascade (derive layers stubbed) -------------------


def test_auto_derives_only_whats_missing(
    populated_db_root: Path, prepared_input: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """With every intermediate present, auto=True derives nothing; with
    the telo and binned files gone, exactly those two layers run."""
    calls: dict[str, int] = {"annotate": 0, "telo": 0, "bin": 0}

    def fake_annotate_batch(**kwargs):
        calls["annotate"] += 1
        return {}

    def fake_telo(path: Path, out: Path, *, motif=None) -> Path:
        calls["telo"] += 1
        out.write_text("")
        return out

    def fake_bin(src: Path, out: Path, *, bin_size, leaf_set, threads) -> None:
        calls["bin"] += 1
        # Rebuild the same binned content the prepared fixture provides.
        label = {"chromosome": ("chr1", "chr2"), "region": ("rA", "rC")}[
            "chromosome" if ".chromosome." in out.name else "region"
        ]
        with gzip.open(out, "wt") as h:
            for ctg, lab in zip(("ctgA", "ctgB"), label, strict=True):
                for i in range(6):
                    h.write(f"{ctg}\t{i * MB}\t{(i + 1) * MB}\t{lab}\n")

    monkeypatch.setattr(sr, "annotate_batch", fake_annotate_batch)
    monkeypatch.setattr(sr, "run_seqtk_telo", fake_telo)
    monkeypatch.setattr(sr, "bin_features", fake_bin)

    out_dir = prepared_input.parent
    _run(prepared_input, populated_db_root, auto=True)
    assert calls == {"annotate": 0, "telo": 0, "bin": 0}

    (out_dir / "samp.telo").unlink()
    for fs in ("chromosome", "region"):
        (out_dir / f"samp.{DUMMY_DB_ID}.{fs}.smoothed.binned{MB}.bed.gz").unlink()
    _run(prepared_input, populated_db_root, auto=True)
    assert calls == {"annotate": 0, "telo": 1, "bin": 2}
