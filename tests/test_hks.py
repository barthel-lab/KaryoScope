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
