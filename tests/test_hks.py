"""Tests for :mod:`karyoscope.core.io.hks`."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.core.external import ToolNotFoundError
from karyoscope.core.io.features import NOVEL_NAME
from karyoscope.core.io.hks import (
    ENV_OVERRIDE,
    _infer_prefix,
    convert_hks_tsv_to_bed,
    get_hks_binary,
    run_hks_lookup,
)


def test_convert_hks_tsv_to_bed_strips_header_and_maps_none(tmp_path: Path) -> None:
    tsv = tmp_path / "raw.tsv"
    tsv.write_text(
        "query_name\tfrom_kmer\tto_kmer\tlabel_name\n"
        "chr1\t0\t10\tchr1\n"
        "chr1\t10\t20\tnone\n"
        "chr1\t20\t30\tcentromere\n"
    )
    bed = tmp_path / "out.bed"
    convert_hks_tsv_to_bed(tsv, bed)

    lines = bed.read_text().splitlines()
    # Header dropped.
    assert lines[0] == "chr1\t0\t10\tchr1"
    # `none` miss label rewritten to the KaryoScope `novel` sentinel.
    assert lines[1] == f"chr1\t10\t20\t{NOVEL_NAME}"
    assert lines[2] == "chr1\t20\t30\tcentromere"
    assert len(lines) == 3


def test_convert_hks_tsv_only_none_in_label_column(tmp_path: Path) -> None:
    """A feature literally containing 'none' as a substring is not corrupted."""
    tsv = tmp_path / "raw.tsv"
    tsv.write_text(
        "query_name\tfrom_kmer\tto_kmer\tlabel_name\nchr1\t0\t10\tnonesuch\nchr1\t10\t20\tnone\n"
    )
    bed = tmp_path / "out.bed"
    convert_hks_tsv_to_bed(tsv, bed)
    lines = bed.read_text().splitlines()
    assert lines[0] == "chr1\t0\t10\tnonesuch"
    assert lines[1] == f"chr1\t10\t20\t{NOVEL_NAME}"


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


def test_convert_hks_tsv_header_only_yields_empty_bed(tmp_path: Path) -> None:
    """A TSV with only a header (no records) produces an empty BED."""
    tsv = tmp_path / "raw.tsv"
    tsv.write_text("query_name\tfrom_kmer\tto_kmer\tlabel_name\n")
    bed = tmp_path / "out.bed"
    convert_hks_tsv_to_bed(tsv, bed)
    assert bed.read_text() == ""


def test_convert_hks_tsv_last_line_without_newline(tmp_path: Path) -> None:
    """A final row lacking a trailing newline is passed through unchanged."""
    tsv = tmp_path / "raw.tsv"
    # No trailing newline on the last line: the miss-label pattern requires a
    # trailing newline, so it is (correctly) not rewritten -- matching the
    # previous splitlines(keepends=True) behavior.
    tsv.write_text("query_name\tfrom_kmer\tto_kmer\tlabel_name\n0\t0\t10\tchr1\n0\t10\t20\tnone")
    bed = tmp_path / "out.bed"
    convert_hks_tsv_to_bed(tsv, bed)
    lines = bed.read_text().splitlines()
    assert lines == ["0\t0\t10\tchr1", "0\t10\t20\tnone"]


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
