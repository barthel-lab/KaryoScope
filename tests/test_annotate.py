"""Unit tests for :mod:`karyoscope.core.annotate` helpers.

The CLI-level end-to-end tests live in ``test_command_annotate.py`` and
need the ``get_featureIDs`` binary; everything here is pure-Python and
runs on the default unit-test pass.
"""

from __future__ import annotations

import json
import os
import threading
from pathlib import Path
from types import SimpleNamespace

import pytest

from karyoscope.core.annotate import (
    _derive_input_basename,
    _detect_worker_death,
    _peak_child_rss_bytes,
    _quiet_worker_pipe_errors,
    _rate,
    _split_combined_bed,
)
from karyoscope.core.io.features import parse_features
from karyoscope.core.io.kmc import (
    clear_combined_marker,
    combined_bed_is_complete,
    combined_marker_path,
    write_combined_marker,
)

# --- _detect_worker_death (pure function, no threads, no os._exit) -----

# Background: ``multiprocessing.Pool`` does not reassign tasks whose
# worker dies from an external signal -- ``imap_unordered`` just blocks
# forever. The watchdog in :mod:`karyoscope.core.annotate` polls this
# helper to detect the death within a few seconds and fail loudly.
# These tests cover the detection logic without involving real Pool
# objects or threads (the ``os._exit`` path can only be verified by
# manual / subprocess testing -- documented in CHANGELOG).


def _fake_pool(*workers: tuple[int | None, int | None]) -> SimpleNamespace:
    """Build a Pool-shaped stand-in. Each worker is ``(pid, exitcode)``.

    The real :func:`_detect_worker_death` only reads ``pool._pool`` and
    each worker's ``.pid`` / ``.exitcode``, so ``SimpleNamespace`` is
    sufficient.
    """
    return SimpleNamespace(_pool=[SimpleNamespace(pid=pid, exitcode=ec) for pid, ec in workers])


def test_detect_worker_death_returns_none_when_all_alive() -> None:
    pool = _fake_pool((100, None), (101, None), (102, None))
    assert _detect_worker_death(pool, {100, 101, 102}) is None


def test_detect_worker_death_catches_fresh_dead_worker() -> None:
    """A worker with non-None exitcode is reported even before the pool
    replaces it."""
    pool = _fake_pool((100, None), (101, -9), (102, None))
    result = _detect_worker_death(pool, {100, 101, 102})
    assert result == ([101], [-9])


def test_detect_worker_death_catches_pool_replacement() -> None:
    """When the pool's _handle_workers thread has already replaced a
    dead worker, the new PID is evidence of the original's death."""
    # 101 died; pool spawned 200 to replace it.
    pool = _fake_pool((100, None), (200, None), (102, None))
    result = _detect_worker_death(pool, {100, 101, 102})
    # No fresh-dead worker (exitcodes empty); died inferred from
    # missing PID.
    assert result == ([101], [])


def test_detect_worker_death_reports_multiple_deaths() -> None:
    pool = _fake_pool((100, None), (101, -9), (102, -15))
    result = _detect_worker_death(pool, {100, 101, 102})
    assert result == ([101, 102], [-15, -9])


def test_detect_worker_death_handles_simultaneous_death_and_replacement() -> None:
    """Mixed: one dead-not-yet-replaced AND one already-replaced.

    Prefers the precise (fresh-dead) info in that case, since we have
    its exitcode -- the inference fallback only kicks in when there are
    zero fresh-dead workers to report on.
    """
    # 101 died and was replaced by 200; 102 died and not yet replaced.
    pool = _fake_pool((100, None), (200, None), (102, -9))
    result = _detect_worker_death(pool, {100, 101, 102})
    # Reports the worker we actually still have evidence for.
    assert result == ([102], [-9])


def test_detect_worker_death_ignores_zero_exitcode() -> None:
    """A worker with exitcode 0 ought to be unusual during normal
    operation, but is still a failure mode worth surfacing if it
    happens -- the pool hangs the same way regardless of the signal."""
    pool = _fake_pool((100, None), (101, 0))
    result = _detect_worker_death(pool, {100, 101})
    # Reports the pid; exitcode 0 is included in the exitcodes list.
    assert result == ([101], [0])


# --- combined-BED completion marker (pure-Python; no binary) -----------
#
# The marker is what makes "reuse an existing combined BED" safe: it is
# written only after get_featureIDs exits 0, and records the BED's size
# and mtime so a file that was truncated by a killed run (or otherwise
# modified after the fact) is never silently trusted. These tests pin
# the fail-closed behaviour.


def _make_combined(tmp_path: Path) -> Path:
    bed = tmp_path / "x.KS_test.combined.presmoothed.featureIDs.bed"
    bed.write_text("seqA\t0\t10\t1\n")
    return bed


def _write_marker(bed: Path) -> None:
    write_combined_marker(
        bed, prefix="x.KS_test", db_path=Path("/db/features"), input_path=Path("in.fa")
    )


def test_marker_roundtrip_marks_complete(tmp_path: Path) -> None:
    bed = _make_combined(tmp_path)
    assert not combined_bed_is_complete(bed)  # no marker yet
    _write_marker(bed)
    assert combined_bed_is_complete(bed)


def test_missing_bed_is_incomplete(tmp_path: Path) -> None:
    bed = tmp_path / "nope.combined.presmoothed.featureIDs.bed"
    assert not combined_bed_is_complete(bed)


def test_size_change_invalidates_marker(tmp_path: Path) -> None:
    """A combined BED appended to / truncated after the marker was
    written (e.g. a killed run that left a partial file) must not be
    trusted."""
    bed = _make_combined(tmp_path)
    _write_marker(bed)
    assert combined_bed_is_complete(bed)
    with bed.open("a") as fh:
        fh.write("seqB\t0\t5\t2\n")
    assert not combined_bed_is_complete(bed)


def test_mtime_change_invalidates_marker(tmp_path: Path) -> None:
    bed = _make_combined(tmp_path)
    _write_marker(bed)
    st = bed.stat()
    os.utime(bed, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))
    assert not combined_bed_is_complete(bed)


def test_corrupt_marker_is_incomplete(tmp_path: Path) -> None:
    bed = _make_combined(tmp_path)
    _write_marker(bed)
    combined_marker_path(bed).write_text("{ not valid json")
    assert not combined_bed_is_complete(bed)


def test_wrong_schema_is_incomplete(tmp_path: Path) -> None:
    bed = _make_combined(tmp_path)
    _write_marker(bed)
    marker = combined_marker_path(bed)
    data = json.loads(marker.read_text())
    data["schema"] = 999
    marker.write_text(json.dumps(data))
    assert not combined_bed_is_complete(bed)


def test_clear_marker_is_idempotent(tmp_path: Path) -> None:
    bed = _make_combined(tmp_path)
    _write_marker(bed)
    assert combined_marker_path(bed).is_file()
    clear_combined_marker(bed)
    assert not combined_marker_path(bed).is_file()
    clear_combined_marker(bed)  # second call must not raise


# --- _quiet_worker_pipe_errors -----------------------------------------


def _raise_in_thread(exc: BaseException) -> None:
    t = threading.Thread(target=lambda: (_ for _ in ()).throw(exc))
    t.start()
    t.join()


def test_quiet_worker_pipe_errors_suppresses_and_restores() -> None:
    """BrokenPipe/EOF from threads are swallowed inside the context;
    everything else is delegated to the previous hook, which is restored
    on exit."""
    original = threading.excepthook
    seen: list[type] = []

    def record(args: threading.ExceptHookArgs) -> None:
        if args.exc_type is not None:
            seen.append(args.exc_type)

    threading.excepthook = record
    try:
        with _quiet_worker_pipe_errors():
            _raise_in_thread(BrokenPipeError())
            _raise_in_thread(EOFError())
            _raise_in_thread(ValueError("real failure"))
        # Hook restored to whatever it was on entry (our `record`).
        assert threading.excepthook is record
    finally:
        threading.excepthook = original

    # Benign pipe errors suppressed; the real error still reached the
    # delegate.
    assert BrokenPipeError not in seen
    assert EOFError not in seen
    assert ValueError in seen


# --- query-k resolution (variable-k support) -------------------------


from karyoscope.core.annotate import _resolve_query_k  # noqa: E402
from karyoscope.exceptions import KaryoscopeError  # noqa: E402


def _manifest(size: int, type_: str, max_size: int) -> SimpleNamespace:
    return SimpleNamespace(kmer=SimpleNamespace(size=size, type=type_, max_size=max_size))


def test_resolve_query_k_defaults_to_manifest_size() -> None:
    assert _resolve_query_k(_manifest(31, "fixed", 31), None, "db") == 31


def test_resolve_query_k_same_as_size_always_ok() -> None:
    assert _resolve_query_k(_manifest(31, "fixed", 31), 31, "db") == 31


def test_resolve_query_k_variable_allows_smaller() -> None:
    assert _resolve_query_k(_manifest(31, "variable", 31), 21, "db") == 21


def test_resolve_query_k_fixed_rejects_override() -> None:
    with pytest.raises(KaryoscopeError, match="fixed-k index"):
        _resolve_query_k(_manifest(31, "fixed", 31), 21, "db")


def test_resolve_query_k_variable_rejects_above_max() -> None:
    with pytest.raises(KaryoscopeError, match="exceeds"):
        _resolve_query_k(_manifest(31, "variable", 31), 40, "db")


def test_resolve_query_k_rejects_below_one() -> None:
    with pytest.raises(KaryoscopeError, match=">= 1"):
        _resolve_query_k(_manifest(31, "variable", 31), 0, "db")


def _read_bed(path: Path) -> list[tuple[str, int, int, str]]:
    """Read a plain BED and return parsed records."""
    out: list[tuple[str, int, int, str]] = []
    for line in path.read_text().splitlines():
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


# --- run-resource reporting -------------------------------------------


def test_rate_reports_gb_per_second() -> None:
    assert _rate(2 * 1024**3, 2.0) == "1.00 GB/s"


def test_rate_declines_to_divide_by_zero() -> None:
    """A phase that measured as instantaneous has no meaningful throughput."""
    assert _rate(1024, 0.0) == "-"
    assert _rate(1024, -1.0) == "-"


def test_peak_child_rss_is_reported_in_bytes() -> None:
    """The units differ by platform; the caller must not have to know that.

    Linux's ru_maxrss is kilobytes and macOS's is bytes, so a helper that
    passed the raw value through would be wrong by 1024x on one of them.
    """
    peak = _peak_child_rss_bytes()
    if peak is None:  # pragma: no cover — platform without `resource`
        pytest.skip("resource module unavailable")
    assert peak >= 0
    # pytest has already reaped children, so this is a plausible RSS in
    # bytes and an implausible one in kilobytes.
    assert peak < 1024**4


# --- CRAM input ------------------------------------------------------------


def test_cram_is_recognised_as_a_read_level_input() -> None:
    """CRAM holds reads, so it takes the read dispatch path, not the assembly one."""
    from karyoscope.core.annotate import _is_reads_input

    assert _is_reads_input(Path("sample.cram")) is True
    assert _is_reads_input(Path("sample.bam")) is True
    assert _is_reads_input(Path("assembly.fa.gz")) is False


def test_cram_extension_is_stripped_from_the_output_basename() -> None:
    """Outputs are named after the input stem, so .cram must be a known suffix."""
    from karyoscope.core.annotate import _derive_input_basename

    assert _derive_input_basename(Path("ALTseq_00001A_tumor.cram")) == "ALTseq_00001A_tumor"


def test_cram_input_declares_a_samtools_dependency() -> None:
    """The preflight resolves tools from the input format rather than assuming."""
    from karyoscope.core.annotate import _annotate_dependencies

    needed = _annotate_dependencies(
        index_type="hks", input_path=Path("s.cram"), bgzip=False
    )
    assert "samtools" in needed
    assert "hks" in needed


def test_cram_bases_estimate_exceeds_bam_per_byte() -> None:
    """CRAM packs far harder than BAM, so a byte of it means far more sequence.

    Sharing BAM's ~1 base/byte would under-estimate a CRAM run's output by
    most of an order of magnitude -- and this figure feeds the disk check
    whose job is to refuse a run that would fill the filesystem.
    """
    from karyoscope.core.annotate import _BASES_PER_BAM_BYTE, _BASES_PER_CRAM_BYTE

    assert _BASES_PER_CRAM_BYTE > _BASES_PER_BAM_BYTE * 5


def test_estimate_input_bases_uses_the_cram_factor(tmp_path: Path) -> None:
    """A .cram file's size is scaled by the CRAM factor, not the FASTA fraction."""
    from karyoscope.core.annotate import _BASES_PER_CRAM_BYTE, estimate_input_bases

    cram = tmp_path / "s.cram"
    cram.write_bytes(b"x" * 1000)
    assert estimate_input_bases(cram) == int(1000 * _BASES_PER_CRAM_BYTE)


def test_query_names_is_refused_for_read_level_input() -> None:
    """Reads cannot carry names: hks buffers all of them and the run is OOM-killed.

    Measured on a 64x human WGS CRAM -- 1.315 G names at ~68 bytes is ~90 GB on
    top of the index. The refusal is deliberate rather than a warning, because
    the failure lands ~20 minutes in, after the decode has been paid for, with
    nothing written and nothing to resume from.
    """
    from karyoscope.core.annotate import _reject_query_names_for_reads

    for name in ("s.cram", "s.bam", "s.fastq.gz", "s.fq"):
        with pytest.raises(KaryoscopeError, match="not supported for read-level input"):
            _reject_query_names_for_reads([Path(name)], True)


def test_query_names_is_allowed_for_assemblies() -> None:
    """An assembly has a few thousand contig names; buffering those is free."""
    from karyoscope.core.annotate import _reject_query_names_for_reads

    _reject_query_names_for_reads([Path("asm.fa.gz")], True)
    _reject_query_names_for_reads([Path("asm.fasta")], True)


def test_reads_may_still_opt_out_of_names_explicitly() -> None:
    """--no-query-names and the default are both fine for reads; only True is refused."""
    from karyoscope.core.annotate import _reject_query_names_for_reads

    _reject_query_names_for_reads([Path("s.cram")], False)
    _reject_query_names_for_reads([Path("s.cram")], None)


def test_query_names_refusal_names_the_offending_inputs() -> None:
    """A cohort run should say WHICH inputs are the problem, not just that one is."""
    from karyoscope.core.annotate import _reject_query_names_for_reads

    with pytest.raises(KaryoscopeError) as excinfo:
        _reject_query_names_for_reads([Path("asm.fa"), Path("tumor.cram")], True)
    assert "tumor.cram" in str(excinfo.value)
    assert "asm.fa" not in str(excinfo.value)


def _hks_cmd_for(tmp_path: Path, monkeypatch, *, input_name: str, query_names):
    """Run the HKS backend with hks stubbed, returning the cmd it would invoke."""
    from karyoscope.core import annotate as ann_mod

    captured: dict = {}

    def _fake_batch(**kwargs):
        captured["report_query_names"] = kwargs["report_query_names"]
        # the backend checks the lookup output exists before smoothing
        for _inp, out in kwargs["io_pairs"]:
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text("")

    monkeypatch.setattr(ann_mod, "run_hks_lookup_batch", _fake_batch)
    monkeypatch.setattr(ann_mod, "run_hks_smooth", lambda **kw: None)

    class _Manifest:
        class index:
            basename = "features"

    inp = tmp_path / input_name
    inp.write_text("")
    ann_mod._run_hks_backend(
        manifest=_Manifest(),
        db_dir=tmp_path,
        input_paths=[inp],
        prefixes={inp: "p"},
        output_dir=tmp_path / "out",
        requested=["cytoband"],
        smooth=False,
        keep_presmoothed=True,
        presmoothed_by_input={inp: {"cytoband": tmp_path / "out" / "p.cytoband.bed"}},
        smoothed_by_input={inp: {}},
        threads=1,
        k=31,
        query_names=query_names,
    )
    return captured["report_query_names"]


def test_reads_default_to_ranks_with_no_flag_at_all(tmp_path: Path, monkeypatch) -> None:
    """Passing neither flag must not silently hand hks the OOM-inducing option.

    This is the path almost every user takes, so the safe behaviour has to be
    the DEFAULT rather than something opted into.
    """
    for name in ("s.cram", "s.bam", "s.fastq"):
        assert _hks_cmd_for(tmp_path, monkeypatch, input_name=name, query_names=None) is False


def test_assemblies_still_default_to_names(tmp_path: Path, monkeypatch) -> None:
    """Contig names are the useful identifier for an assembly and cost nothing."""
    assert _hks_cmd_for(tmp_path, monkeypatch, input_name="asm.fa", query_names=None) is True


def test_explicit_no_query_names_forces_ranks_on_an_assembly(tmp_path: Path, monkeypatch) -> None:
    """The override works downward too, for databases with long feature names."""
    assert _hks_cmd_for(tmp_path, monkeypatch, input_name="asm.fa", query_names=False) is False
