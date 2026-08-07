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


def _capture_lookup_cmd(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, report_query_names: bool
) -> list[str]:
    """Run ``run_hks_lookup`` with the binary + subprocess stubbed, return the cmd."""
    captured: dict[str, list[str]] = {}
    monkeypatch.setattr("karyoscope.core.io.hks.get_hks_binary", lambda: "hks")
    monkeypatch.setattr(
        "karyoscope.core.io.hks.run_tool",
        lambda cmd, capture=False: captured.__setitem__("cmd", cmd),
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
        lambda cmd, capture=False: captured.__setitem__("cmd", cmd),
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
        lambda cmd, capture=False: captured.__setitem__("cmd", [str(c) for c in cmd]),
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


def _reference_convert(tsv_text: str) -> str:
    """The original line-by-line implementation, as an oracle."""
    miss = f"\t{_HKS_MISS_LABEL}\n"
    novel = f"\t{NOVEL_NAME}\n"
    lines = tsv_text.splitlines(keepends=True)
    return "".join(line.replace(miss, novel) for line in lines[1:])
