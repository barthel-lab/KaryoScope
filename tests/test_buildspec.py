"""Tests for :mod:`karyoscope.core.buildspec`."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.core.buildspec import DEFAULT_BACKGROUND, BuildSpec, FeatureSetSpec
from karyoscope.exceptions import BuildError


@pytest.fixture
def genome(tmp_path: Path) -> Path:
    p = tmp_path / "g.fa"
    p.write_text(">chr1\nACGT\n")
    return p


@pytest.fixture
def bed(tmp_path: Path) -> Path:
    p = tmp_path / "r.bed"
    p.write_text("chr1\t0\t4\tLINE\n")
    return p


def test_from_flags_defaults_background(genome: Path, bed: Path) -> None:
    spec = BuildSpec.from_flags(
        db_id="HKS_x", version="1.0.0", sequence=genome, feature_beds={"repeat": bed}
    )
    (fs,) = spec.feature_sets
    assert fs.mode == "bed"
    assert fs.background == DEFAULT_BACKGROUND
    assert spec.s == 31


def test_from_yaml_parses_all_fields(tmp_path: Path, genome: Path, bed: Path) -> None:
    spec_file = tmp_path / "spec.yaml"
    spec_file.write_text(
        f"""
id: HKS_y
version: "2.0.0"
sequence: {genome}
kmer: {{ s: 21 }}
build: {{ threads: 8, mem_gigas: 4 }}
feature_sets:
  - name: repeat
    bed: {bed}
    background: nonrepeat
  - name: chrom
    bed: {bed}
    background: null
roles: {{ chromosome_assignment: chrom }}
exclude: [ChrM, ChrC]
"""
    )
    spec = BuildSpec.from_yaml(spec_file)
    assert spec.id == "HKS_y"
    assert spec.version == "2.0.0"
    assert spec.s == 21
    assert spec.threads == 8
    assert spec.roles == {"chromosome_assignment": "chrom"}
    assert spec.exclude == ["ChrM", "ChrC"]
    by_name = {fs.name: fs for fs in spec.feature_sets}
    assert by_name["repeat"].background == "nonrepeat"
    assert by_name["chrom"].background is None  # explicit null disables


def test_from_yaml_relative_paths_resolved(tmp_path: Path) -> None:
    (tmp_path / "g.fa").write_text(">c\nACGT\n")
    (tmp_path / "r.bed").write_text("c\t0\t4\tL\n")
    spec_file = tmp_path / "spec.yaml"
    spec_file.write_text(
        'id: X\nversion: "1"\nsequence: g.fa\nfeature_sets:\n  - name: r\n    bed: r.bed\n'
    )
    spec = BuildSpec.from_yaml(spec_file)
    assert spec.sequence == tmp_path / "g.fa"
    assert spec.feature_sets[0].bed == tmp_path / "r.bed"


def test_validate_rejects_duplicate_names(genome: Path, bed: Path) -> None:
    spec = BuildSpec(
        id="x",
        version="1",
        sequence=genome,
        feature_sets=[FeatureSetSpec("r", bed=bed), FeatureSetSpec("r", bed=bed)],
    )
    with pytest.raises(BuildError, match="duplicate feature-set name"):
        spec.validate()


def test_validate_rejects_background_in_mode_b(genome: Path) -> None:
    spec = BuildSpec(
        id="x",
        version="1",
        feature_sets=[FeatureSetSpec("r", fastas=[genome], background="bg")],
    )
    with pytest.raises(BuildError, match="only applies to BED"):
        spec.validate()


def test_validate_requires_sequence_for_bed(bed: Path) -> None:
    spec = BuildSpec(
        id="x", version="1", sequence=None, feature_sets=[FeatureSetSpec("r", bed=bed)]
    )
    with pytest.raises(BuildError, match="no 'sequence'"):
        spec.validate()


def test_validate_rejects_variable_k_with_priority(tmp_path: Path, genome: Path, bed: Path) -> None:
    prio = tmp_path / "p.txt"
    prio.write_text("LINE 1\n")
    spec = BuildSpec(
        id="x",
        version="1",
        sequence=genome,
        feature_sets=[FeatureSetSpec("r", bed=bed, priority=prio, variable_k=True)],
    )
    with pytest.raises(BuildError, match="mutually exclusive"):
        spec.validate()


def test_validate_rejects_reserved_background_none(genome: Path, bed: Path) -> None:
    spec = BuildSpec(
        id="x",
        version="1",
        sequence=genome,
        feature_sets=[FeatureSetSpec("r", bed=bed, background="none")],
    )
    with pytest.raises(BuildError, match="reserved"):
        spec.validate()


def test_validate_rejects_bad_feature_set_name(genome: Path, bed: Path) -> None:
    spec = BuildSpec(
        id="x",
        version="1",
        sequence=genome,
        feature_sets=[FeatureSetSpec("a/b", bed=bed)],
    )
    with pytest.raises(BuildError, match="path/whitespace"):
        spec.validate()


def test_validate_accepts_variable_k_in_mode_a(genome: Path, bed: Path) -> None:
    # variable_k is supported for BED (mode A): the base is later built from the
    # generated per-feature FASTAs. validate() must not reject it.
    spec = BuildSpec(
        id="x",
        version="1",
        sequence=genome,
        feature_sets=[FeatureSetSpec("r", bed=bed, variable_k=True)],
    )
    spec.validate()  # does not raise
