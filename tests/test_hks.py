"""Tests for :mod:`karyoscope.core.io.hks`."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.core.external import ToolNotFoundError
from karyoscope.core.io.features import NOVEL_NAME
from karyoscope.core.io.hks import (
    _HKS_MISS_LABEL,
    ENV_OVERRIDE,
    _infer_prefix,
    get_hks_binary,
    run_hks_lookup,
    run_hks_smooth,
)
from karyoscope.exceptions import KaryoscopeError


def _capture_lookup_cmd(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, report_query_names: bool
) -> list[str]:
    """Run ``run_hks_lookup`` with the binary + subprocess stubbed, return the cmd."""
    captured: dict[str, list[str]] = {}
    monkeypatch.setattr("karyoscope.core.io.hks.get_hks_binary", lambda: "hks")
    monkeypatch.setattr(
        "karyoscope.core.io.hks.run_tool",
        lambda cmd, capture=False, **kw: captured.__setitem__("cmd", cmd),
    )
    run_hks_lookup(
        base_path=tmp_path / "features.hksb",
        feature_set_file=tmp_path / "features.chromosome.hksf",
        k=31,
        input_path=tmp_path / "input.fa",
        output_path=tmp_path / "out.tsv",
        report_query_names=report_query_names,
    )
    return captured["cmd"]


def test_run_hks_lookup_reports_names_when_requested(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Assemblies pass --report-query-names so contig names reach the output."""
    cmd = _capture_lookup_cmd(tmp_path, monkeypatch, report_query_names=True)
    assert "--report-query-names" in cmd
    assert "--report-misses" in cmd


def test_run_hks_lookup_omits_names_for_reads(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reads drop --report-query-names (integer ranks) but still report misses."""
    cmd = _capture_lookup_cmd(tmp_path, monkeypatch, report_query_names=False)
    assert "--report-query-names" not in cmd
    assert "--report-misses" in cmd


def _capture_smooth_cmd(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, out: Path | None = None
) -> list[str]:
    """Run ``run_hks_smooth`` with the binary + subprocess stubbed, return the cmd."""
    captured: dict[str, list[str]] = {}
    monkeypatch.setattr("karyoscope.core.io.hks.get_hks_binary", lambda: "hks")
    monkeypatch.setattr(
        "karyoscope.core.io.hks.run_tool",
        lambda cmd, capture=False, **kw: captured.__setitem__("cmd", cmd),
    )
    run_hks_smooth(
        hierarchy_file=tmp_path / "features.chromosome.hierarchy.txt",
        input_path=tmp_path / "raw.tsv",
        output_path=out if out is not None else tmp_path / "smoothed.bed",
    )
    return captured["cmd"]


def _flag_value(cmd: list[str], flag: str) -> str:
    return cmd[cmd.index(flag) + 1]


def test_run_hks_smooth_writes_the_bed_directly(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """smooth is pointed at the final BED, not a temp TSV needing a rewrite."""
    out = tmp_path / "smoothed.bed"
    cmd = _capture_smooth_cmd(tmp_path, monkeypatch, out=out)
    assert _flag_value(cmd, "-o") == str(out)
    assert "--no-header" in cmd
    assert _flag_value(cmd, "--miss-label") == NOVEL_NAME


def test_run_hks_lookup_writes_the_presmoothed_bed_directly(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The lookup output is the presmoothed BED, so it carries no header."""
    cmd = _capture_lookup_cmd(tmp_path, monkeypatch, report_query_names=True)
    assert "--no-header" in cmd


def test_lookup_and_smooth_agree_on_the_label_vocabulary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Neither may ask for numeric label ids.

    With no header in between to declare which vocabulary the fourth column
    uses, this is the only thing keeping the two in step: smooth resolves
    those tokens against the hierarchy, and would reject every name if it
    expected ids.
    """
    lookup = _capture_lookup_cmd(tmp_path, monkeypatch, report_query_names=True)
    smooth = _capture_smooth_cmd(tmp_path, monkeypatch)
    assert "--report-label-ids" not in lookup
    assert "--report-label-ids" not in smooth


def test_lookup_and_smooth_agree_on_the_miss_label(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The invariant the whole scheme rests on.

    ``hks smooth`` parses the miss token out of its input as well as writing
    it, so if lookup wrote a different one every miss run would be read as an
    unknown feature name. Neither may quietly fall back to HKS's ``none``.
    """
    lookup = _capture_lookup_cmd(tmp_path, monkeypatch, report_query_names=True)
    smooth = _capture_smooth_cmd(tmp_path, monkeypatch)
    assert _flag_value(lookup, "--miss-label") == _flag_value(smooth, "--miss-label")
    assert _flag_value(lookup, "--miss-label") == NOVEL_NAME
    assert _flag_value(lookup, "--miss-label") != _HKS_MISS_LABEL


def test_run_hks_lookup_emits_the_karyoscope_miss_label(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Misses come out of lookup already labelled ``novel``, not ``none``."""
    cmd = _capture_lookup_cmd(tmp_path, monkeypatch, report_query_names=True)
    assert _flag_value(cmd, "--miss-label") == NOVEL_NAME


def test_infer_prefix_strips_extension() -> None:
    assert _infer_prefix(Path("sample.fasta.gz"), "HKS_db") == "sample.HKS_db"
    assert _infer_prefix(Path("reads.fq"), "HKS_db") == "reads.HKS_db"


def test_get_hks_binary_env_override(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fake = tmp_path / "hks"
    fake.write_text("#!/bin/sh\n")
    monkeypatch.setenv(ENV_OVERRIDE, str(fake))
    assert get_hks_binary() == str(fake)


def test_get_hks_binary_env_override_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv(ENV_OVERRIDE, str(tmp_path / "does_not_exist"))
    with pytest.raises(ToolNotFoundError, match="no file exists"):
        get_hks_binary()


def test_get_hks_binary_not_found(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(ENV_OVERRIDE, raising=False)
    monkeypatch.setattr("karyoscope.core.io.hks.shutil.which", lambda _name: None)
    with pytest.raises(ToolNotFoundError, match="was not found on PATH"):
        get_hks_binary()


# --- index construction wrappers -------------------------------------


def _capture_cmd(monkeypatch: pytest.MonkeyPatch) -> dict:
    captured: dict = {}
    monkeypatch.setattr("karyoscope.core.io.hks.get_hks_binary", lambda: "hks")
    monkeypatch.setattr(
        "karyoscope.core.io.hks.run_tool",
        lambda cmd, capture=False, **kw: captured.__setitem__("cmd", [str(c) for c in cmd]),
    )
    return captured


def test_build_base_input_and_external_memory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from karyoscope.core.io.hks import run_hks_build_base

    captured = _capture_cmd(monkeypatch)
    run_hks_build_base(
        output_path=tmp_path / "features.hksb",
        s=31,
        input_path=tmp_path / "g.fa",
        threads=8,
        mem_gigas=16,
        external_memory=tmp_path / "scratch",
    )
    cmd = captured["cmd"]
    assert cmd[:3] == ["hks", "build-base", "-s"]
    assert "--input" in cmd and "--input-file-list" not in cmd
    assert "--external-memory" in cmd
    assert cmd[cmd.index("--mem-gigas") + 1] == "16"


def test_build_base_requires_exactly_one_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from karyoscope.core.io.hks import run_hks_build_base

    _capture_cmd(monkeypatch)
    with pytest.raises(ValueError, match="exactly one"):
        run_hks_build_base(output_path=tmp_path / "b.hksb", s=31)


def test_add_feature_set_priority_mode(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from karyoscope.core.io.hks import run_hks_add_feature_set

    captured = _capture_cmd(monkeypatch)
    run_hks_add_feature_set(
        base_path=tmp_path / "b.hksb",
        output_path=tmp_path / "b.repeat.hksf",
        feature_set_name="repeat",
        feature_names=tmp_path / "names.txt",
        feature_hierarchy=tmp_path / "h.txt",
        feature_file_list=tmp_path / "fof.txt",
        feature_priorities=tmp_path / "prio.txt",
    )
    cmd = captured["cmd"]
    assert cmd[:2] == ["hks", "add-feature-set"]
    assert "--feature-priorities" in cmd
    assert "--variable-k-support" not in cmd
    assert cmd[cmd.index("--feature-set-name") + 1] == "repeat"


def test_add_feature_set_variable_k_and_priority_conflict(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from karyoscope.core.io.hks import run_hks_add_feature_set

    _capture_cmd(monkeypatch)
    with pytest.raises(ValueError, match="mutually exclusive"):
        run_hks_add_feature_set(
            base_path=tmp_path / "b.hksb",
            output_path=tmp_path / "o.hksf",
            feature_set_name="x",
            feature_names=tmp_path / "n.txt",
            feature_hierarchy=tmp_path / "h.txt",
            feature_file_list=tmp_path / "fof.txt",
            feature_priorities=tmp_path / "prio.txt",
            variable_k_support=True,
        )


def test_validate_sibling_priorities() -> None:
    from karyoscope.core.io.hks import validate_sibling_priorities

    parent_of = {"a": "root", "b": "root", "c": "root"}
    # all distinct -> ok
    assert validate_sibling_priorities(parent_of, {"a": 1, "b": 2, "c": 3}) == []
    # all equal -> ok
    assert validate_sibling_priorities(parent_of, {"a": 0, "b": 0, "c": 0}) == []
    # mixed (two share, one distinct) -> flagged
    issues = validate_sibling_priorities(parent_of, {"a": 1, "b": 1, "c": 2})
    assert len(issues) == 1 and "mixed priorities" in issues[0]


# --- batch lookup -----------------------------------------------------
#
# The batch path is a second copy of the same pipeline, so it can drift
# from the single-input one silently. These pin the properties that
# would break correctness rather than merely differ.


def _capture_batch_cmd(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, n_inputs: int = 2
) -> tuple[list[str], list[tuple[Path, Path]]]:
    from karyoscope.core.io.hks import run_hks_lookup_batch

    captured: dict[str, list[str]] = {}
    monkeypatch.setattr("karyoscope.core.io.hks.get_hks_binary", lambda: "hks")
    monkeypatch.setattr(
        "karyoscope.core.io.hks.run_tool",
        lambda cmd, capture=False, **kw: captured.__setitem__("cmd", cmd),
    )
    io_pairs = [(tmp_path / f"in{i}.fa", tmp_path / f"out{i}.bed") for i in range(n_inputs)]
    run_hks_lookup_batch(
        base_path=tmp_path / "features.hksb",
        feature_set_file=tmp_path / "features.chromosome.hksf",
        k=31,
        io_pairs=io_pairs,
        report_query_names=True,
    )
    return captured["cmd"], io_pairs


def test_batch_lookup_writes_the_presmoothed_bed_format(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The batch path must produce the same bytes as the single-input one.

    Its outputs are read back by the same `hks smooth` invocation, so a
    batch lookup that kept HKS's defaults would hand smooth a header and a
    `none` token it is not expecting.
    """
    cmd, _ = _capture_batch_cmd(tmp_path, monkeypatch)
    assert _flag_value(cmd, "--miss-label") == NOVEL_NAME
    assert "--no-header" in cmd
    assert "--report-label-ids" not in cmd


def test_batch_lookup_pairs_every_query_with_its_own_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """-q and -o are positional pairs; a mismatch would cross-write outputs."""
    cmd, io_pairs = _capture_batch_cmd(tmp_path, monkeypatch, n_inputs=3)
    pairs = [
        (cmd[i + 1], cmd[i + 3])
        for i, tok in enumerate(cmd)
        if tok == "-q" and i + 3 < len(cmd) and cmd[i + 2] == "-o"
    ]
    assert pairs == [(str(inp), str(out)) for inp, out in io_pairs]


def test_batch_lookup_matches_the_single_input_output_flags(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Whatever governs output format must be identical across both paths."""
    single = _capture_lookup_cmd(tmp_path, monkeypatch, report_query_names=True)
    batch, _ = _capture_batch_cmd(tmp_path, monkeypatch)
    fmt = ("--miss-label", "--no-header", "--report-misses", "--report-query-names")
    for flag in fmt:
        in_single = flag in single
        in_batch = flag in batch
        assert in_single == in_batch, flag
        if in_single and flag == "--miss-label":
            assert _flag_value(single, flag) == _flag_value(batch, flag)


def test_batch_lookup_with_no_inputs_does_nothing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from karyoscope.core.io.hks import run_hks_lookup_batch

    calls: list[object] = []
    monkeypatch.setattr("karyoscope.core.io.hks.get_hks_binary", lambda: "hks")
    monkeypatch.setattr(
        "karyoscope.core.io.hks.run_tool", lambda cmd, capture=False, **kw: calls.append(cmd)
    )
    run_hks_lookup_batch(
        base_path=tmp_path / "b.hksb",
        feature_set_file=tmp_path / "b.chromosome.hksf",
        k=31,
        io_pairs=[],
    )
    assert calls == []


# --- CRAM / BAM conversion -------------------------------------------------


def _capture_samtools_cmd(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    suffix: str,
    reference: Path | None,
) -> list[str]:
    """Run a lookup over a BAM/CRAM input, returning the samtools cmd it built.

    ``hks`` itself is stubbed out; the assertion target is the conversion step
    that runs before it.
    """
    captured: dict[str, list[str]] = {}

    class _Result:
        returncode = 0
        stderr = b""

    def _fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        # The caller sizes and then hands the temp FASTA to hks, so it must exist.
        stdout = kwargs.get("stdout")
        if stdout is not None:
            stdout.close()
        return _Result()

    monkeypatch.setattr("karyoscope.core.io.hks.get_hks_binary", lambda: "hks")
    monkeypatch.setattr("karyoscope.core.io.hks.require_tool", lambda *a, **kw: "samtools")
    monkeypatch.setattr("karyoscope.core.io.hks.subprocess.run", _fake_run)
    monkeypatch.setattr("karyoscope.core.io.hks.run_tool", lambda cmd, **kw: None)

    aln = tmp_path / f"input{suffix}"
    aln.write_bytes(b"")
    run_hks_lookup(
        base_path=tmp_path / "features.hksb",
        feature_set_file=tmp_path / "features.cytoband.hksf",
        k=31,
        input_path=aln,
        output_path=tmp_path / "out.bed",
        reference=reference,
    )
    return captured["cmd"]


def test_cram_conversion_passes_the_reference(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """CRAM decode supplies --reference, not -T (which is samtools' taglist)."""
    ref = tmp_path / "genome.fasta"
    ref.write_text(">chr1\nACGT\n")
    cmd = _capture_samtools_cmd(tmp_path, monkeypatch, suffix=".cram", reference=ref)
    assert "--reference" in cmd
    assert cmd[cmd.index("--reference") + 1] == str(ref)
    # -T here would silently mean "copy these tags to the header" and leave the
    # decode with no reference at all.
    assert "-T" not in cmd


def test_cram_without_reference_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A CRAM cannot be decoded without its reference, so the run stops early."""
    with pytest.raises(KaryoscopeError, match="CRAM"):
        _capture_samtools_cmd(tmp_path, monkeypatch, suffix=".cram", reference=None)


def test_bam_conversion_needs_no_reference(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """BAM is self-contained; no --reference is added and no error is raised."""
    cmd = _capture_samtools_cmd(tmp_path, monkeypatch, suffix=".bam", reference=None)
    assert "--reference" not in cmd


def test_conversion_states_flag_filter_and_mate_suffix_explicitly(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """-F 0x900 -N are pinned rather than left to samtools' defaults.

    ``-F 0x900`` keeps exactly one full-length record per read (primaries),
    while ``-N`` forces the /1,/2 suffix that distinguishes mates whose QNAME
    the aligner stripped. Without -N both mates share one name and the pairing
    is unrecoverable downstream.
    """
    cmd = _capture_samtools_cmd(tmp_path, monkeypatch, suffix=".bam", reference=None)
    assert cmd[cmd.index("-F") + 1] == "0x900"
    assert "-N" in cmd


def test_materialised_queries_decodes_once_and_passes_through(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Already-seekable inputs cost nothing; alignments are decoded exactly once.

    The decode used to live inside run_hks_lookup_batch, which the annotate
    backend calls once per feature set -- so a six-set database decoded the same
    CRAM six times, ~25 minutes each on a 56 GB input.
    """
    from karyoscope.core.io.hks import materialised_queries

    calls: list[list[str]] = []

    class _Result:
        returncode = 0
        stderr = b""

    def _fake_run(cmd, **kwargs):
        calls.append(cmd)
        stdout = kwargs.get("stdout")
        if stdout is not None:
            stdout.close()
        return _Result()

    monkeypatch.setattr("karyoscope.core.io.hks.require_tool", lambda *a, **kw: "samtools")
    monkeypatch.setattr("karyoscope.core.io.hks.subprocess.run", _fake_run)

    plain = tmp_path / "reads.fastq"
    plain.write_text("")
    bam = tmp_path / "aln.bam"
    bam.write_bytes(b"")
    with materialised_queries([plain, bam], threads=4) as resolved:
        # a seekable input is handed straight back, unconverted
        assert resolved[plain] == plain
        # the alignment became a temp FASTA
        assert resolved[bam] != bam
        tmp_fasta = resolved[bam]
    assert len(calls) == 1, "exactly one decode for one alignment"
    assert not tmp_fasta.exists(), "temp FASTA removed on exit"


def test_materialised_queries_tees_the_name_sidecar(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The sidecar comes off the SAME decode, not a second pass over the input."""
    from karyoscope.core.io.hks import materialised_queries

    captured: dict[str, object] = {}

    class _Result:
        returncode = 0
        stderr = b""

    def _fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        return _Result()

    monkeypatch.setattr("karyoscope.core.io.hks.require_tool", lambda *a, **kw: "samtools")
    monkeypatch.setattr("karyoscope.core.io.hks.subprocess.run", _fake_run)

    bam = tmp_path / "aln.bam"
    bam.write_bytes(b"")
    sidecar = tmp_path / "names.txt.gz"
    with materialised_queries([bam], threads=2, query_names_sidecar={bam: sidecar}):
        pass
    joined = " ".join(captured["cmd"])
    # one process: samtools teed into the sidecar writer AND the FASTA
    assert "tee" in joined
    assert str(sidecar) in joined
    assert "samtools fasta" in joined
