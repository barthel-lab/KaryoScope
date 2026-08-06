"""End-to-end tests for :func:`karyoscope.core.build.build_database`.

These invoke the real ``hks`` binary and are skipped when it is not on PATH.
The genome is intentionally tiny and ``s`` small so construction is fast.
"""

from __future__ import annotations

import random
import shutil
from pathlib import Path

import pytest

from karyoscope.core.build import build_database
from karyoscope.core.buildspec import BuildSpec, FeatureSetSpec
from karyoscope.core.io.hks import run_hks_lookup
from karyoscope.exceptions import BuildError
from karyoscope.manifest import validate_database_layout

requires_hks = pytest.mark.skipif(shutil.which("hks") is None, reason="hks binary not on PATH")


@pytest.fixture
def tiny_inputs(tmp_path: Path) -> dict:
    rng = random.Random(7)
    chr1 = "".join(rng.choice("ACGT") for _ in range(600))
    chr2 = "".join(rng.choice("ACGT") for _ in range(400))
    genome = tmp_path / "genome.fa"
    genome.write_text(f">chr1\n{chr1}\n>chr2\n{chr2}\n")

    repeat_bed = tmp_path / "repeat.bed"
    repeat_bed.write_text(
        "chr1\t0\t150\tLINE\n"
        "chr1\t100\t250\tSINE\n"  # overlaps LINE
        "chr1\t400\t500\tLINE\n"
        "chr2\t50\t200\tSINE\n"
    )
    priority = tmp_path / "repeat.priority.txt"
    priority.write_text("LINE 1 categorized\nSINE 2 categorized\nnonrepeat 3 categorized\n")

    gene_bed = tmp_path / "gene.bed"
    gene_bed.write_text("chr1\t0\t300\texon\nchr1\t300\t600\tintron\nchr2\t0\t400\texon\n")
    return {
        "genome": genome,
        "repeat_bed": repeat_bed,
        "priority": priority,
        "gene_bed": gene_bed,
    }


@requires_hks
def test_build_database_end_to_end(tmp_path: Path, tiny_inputs: dict) -> None:
    spec = BuildSpec(
        id="HKS_tiny",
        version="1.0.0",
        sequence=tiny_inputs["genome"],
        s=11,
        threads=2,
        mem_gigas=2,
        feature_sets=[
            FeatureSetSpec(
                name="repeat",
                bed=tiny_inputs["repeat_bed"],
                background="nonrepeat",
                priority=tiny_inputs["priority"],
            ),
            FeatureSetSpec(name="gene", bed=tiny_inputs["gene_bed"]),
        ],
    )
    db_root = tmp_path / "db"

    result = build_database(spec, db_root, register=True, force=True)

    db_dir = db_root / "HKS_tiny"
    assert result.db_dir == db_dir.resolve()
    assert result.registered

    # Layout is well-formed and features.tsv is correctly absent.
    manifest = validate_database_layout(db_dir)
    assert manifest.features is None
    assert manifest.feature_sets == ["repeat", "gene"]

    # Expected artifacts exist; the working dir was cleaned up.
    for rel in (
        "manifest.yaml",
        "hierarchy.tsv",
        "colors.tsv",
        "index/features.hksb",
        "index/features.repeat.hksf",
        "index/features.repeat.hierarchy.txt",
        "index/features.gene.hksf",
        "index/features.gene.hierarchy.txt",
    ):
        assert (db_dir / rel).is_file(), rel
    assert not (db_dir / "_build_work").exists()

    # Modes: repeat has a priority file -> priority; gene defaults to fixed.
    modes = {fr.name: fr.mode for fr in result.feature_sets}
    assert modes == {"repeat": "priority", "gene": "fixed"}
    # Gap-fill created the background leaf where the repeat BED left gaps.
    assert "nonrepeat" in dict((fr.name, fr) for fr in result.feature_sets)["repeat"].leaves

    # The index is queryable: annotate the genome and confirm priority resolution.
    out_tsv = tmp_path / "lookup.tsv"
    run_hks_lookup(
        base_path=db_dir / "index/features.hksb",
        feature_set_file=db_dir / "index/features.repeat.hksf",
        k=11,
        input_path=tiny_inputs["genome"],
        output_path=out_tsv,
        report_query_names=True,
        capture=True,
    )
    text = out_tsv.read_text()
    assert "LINE" in text and "SINE" in text and "nonrepeat" in text
    # No k-mer should be a miss: the base index was built from this genome.
    assert "\tnone\n" not in text


@requires_hks
def test_build_variable_k_from_bed_supports_multi_k_query(
    tmp_path: Path, tiny_inputs: dict
) -> None:
    """A variable-k build from mode-A BED input is queryable at k < s.

    This is the whole point of variable-k: one index answers a k-sweep. It works
    for BED input because the base is built from the generated per-feature FASTAs
    (so every feature run starts at a sequence start with a dummy node).
    """
    spec = BuildSpec(
        id="HKS_vk",
        version="1.0.0",
        sequence=tiny_inputs["genome"],
        s=11,
        threads=2,
        mem_gigas=2,
        feature_sets=[FeatureSetSpec(name="gene", bed=tiny_inputs["gene_bed"], variable_k=True)],
    )
    db_root = tmp_path / "db"
    result = build_database(spec, db_root, register=True, force=True)
    assert result.feature_sets[0].mode == "variable-k"

    db_dir = db_root / "HKS_vk"
    manifest = validate_database_layout(db_dir)
    assert manifest.kmer.type == "variable"

    # Query the same index at k = s and at k < s; both must succeed and label.
    for k in (11, 7):
        out_tsv = tmp_path / f"lookup_k{k}.tsv"
        run_hks_lookup(
            base_path=db_dir / "index/features.hksb",
            feature_set_file=db_dir / "index/features.gene.hksf",
            k=k,
            input_path=tiny_inputs["genome"],
            output_path=out_tsv,
            report_query_names=True,
            capture=True,
        )
        text = out_tsv.read_text()
        assert "exon" in text and "intron" in text, f"k={k}: {text[:200]}"


@requires_hks
def test_annotate_k_sweep_on_variable_k_index(tmp_path: Path, tiny_inputs: dict) -> None:
    """`annotate(k=...)` runs a k-sweep on a variable-k index and tags outputs."""
    from karyoscope.core.annotate import annotate
    from karyoscope.exceptions import KaryoscopeError

    spec = BuildSpec(
        id="HKS_vk",
        version="1.0.0",
        sequence=tiny_inputs["genome"],
        s=11,
        threads=2,
        mem_gigas=2,
        feature_sets=[FeatureSetSpec(name="gene", bed=tiny_inputs["gene_bed"], variable_k=True)],
    )
    db_root = tmp_path / "db"
    build_database(spec, db_root, register=True, force=True)

    out = tmp_path / "sweep"
    for query_k in (11, 7):
        annotate(
            input_path=tiny_inputs["genome"],
            output_dir=out,
            db_root=db_root,
            db_id="HKS_vk",
            feature_sets=["gene"],
            smooth=False,
            bgzip=False,
            k=query_k,
        )
        tagged = out / f"genome.HKS_vk.k{query_k}.gene.presmoothed.bed"
        assert tagged.is_file(), f"missing {tagged}"
        assert "exon" in tagged.read_text()

    # A k above the index maximum is rejected.
    with pytest.raises(KaryoscopeError, match="exceeds"):
        annotate(
            input_path=tiny_inputs["genome"],
            output_dir=out,
            db_root=db_root,
            db_id="HKS_vk",
            feature_sets=["gene"],
            smooth=False,
            bgzip=False,
            k=99,
        )


@requires_hks
def test_build_refuses_existing_dir_without_force(tmp_path: Path, tiny_inputs: dict) -> None:
    from karyoscope.exceptions import BuildError

    spec = BuildSpec(
        id="HKS_tiny",
        version="1.0.0",
        sequence=tiny_inputs["genome"],
        s=11,
        threads=2,
        mem_gigas=2,
        feature_sets=[FeatureSetSpec(name="gene", bed=tiny_inputs["gene_bed"])],
    )
    db_root = tmp_path / "db"
    build_database(spec, db_root, register=False, force=False)
    with pytest.raises(BuildError, match="already exists"):
        build_database(spec, db_root, register=False, force=False)


@requires_hks
def test_failed_build_removes_the_partial_database(tmp_path: Path, tiny_inputs: dict) -> None:
    """A failed build must not block the re-run.

    It used to leave an unusable half-database behind, so the obvious next
    move -- fix the input, run the same command again -- hit "database
    directory already exists ... Pass --force", about a directory that never
    held a working database.
    """
    db_root = tmp_path / "db"
    # Siblings whose priorities are neither all-equal nor all-distinct: an
    # HKS constraint, checked during preparation.
    bad_priority = tmp_path / "bad.priority.txt"
    bad_priority.write_text("LINE 1 categorized\nSINE 1 categorized\nnonrepeat 2 categorized\n")
    spec = BuildSpec.from_flags(
        db_id="HKS_fail",
        version="1.0.0",
        sequence=tiny_inputs["genome"],
        feature_beds={"repeat": tiny_inputs["repeat_bed"]},
        priorities={"repeat": bad_priority},
        s=11,
    )
    with pytest.raises(BuildError):
        build_database(spec, db_root, register=False)
    assert not (db_root / "HKS_fail").exists()

    # The re-run with a valid spec now succeeds without --force.
    good = BuildSpec.from_flags(
        db_id="HKS_fail",
        version="1.0.0",
        sequence=tiny_inputs["genome"],
        feature_beds={"repeat": tiny_inputs["repeat_bed"]},
        s=11,
    )
    result = build_database(good, db_root, register=False)
    assert result.db_dir.is_dir()
    validate_database_layout(result.db_dir)


@requires_hks
def test_keep_intermediates_preserves_a_failed_build(tmp_path: Path, tiny_inputs: dict) -> None:
    # The escape hatch for inspecting a build that failed late.
    db_root = tmp_path / "db"
    bad_priority = tmp_path / "bad.priority.txt"
    bad_priority.write_text("LINE 1 categorized\nSINE 1 categorized\nnonrepeat 2 categorized\n")
    spec = BuildSpec.from_flags(
        db_id="HKS_keep",
        version="1.0.0",
        sequence=tiny_inputs["genome"],
        feature_beds={"repeat": tiny_inputs["repeat_bed"]},
        priorities={"repeat": bad_priority},
        s=11,
    )
    with pytest.raises(BuildError):
        build_database(spec, db_root, register=False, keep_intermediates=True)
    assert (db_root / "HKS_keep").is_dir()


def test_priority_file_must_name_every_node(tmp_path: Path, tiny_inputs: dict) -> None:
    """An unlisted node would take priority 0 -- the BEST priority -- and win.

    The node most easily forgotten is the gap-fill leaf, because `build` adds it
    to the hierarchy itself and so it never appears in a hand-written hierarchy
    file. Left unlisted it outranks every real feature.
    """
    prio = tmp_path / "missing_background.priority.txt"
    prio.write_text("LINE 1 categorized\nSINE 2 categorized\n")  # no `nonrepeat`
    spec = BuildSpec.from_flags(
        db_id="HKS_missing_prio",
        version="1.0.0",
        sequence=tiny_inputs["genome"],
        feature_beds={"repeat": tiny_inputs["repeat_bed"]},
        priorities={"repeat": prio},
        s=11,
    )
    with pytest.raises(BuildError) as excinfo:
        build_database(spec, tmp_path / "db", register=False)
    message = str(excinfo.value)
    assert "background" in message
    assert "gap-fill" in message


def test_priority_file_need_not_name_the_root(tmp_path: Path, tiny_inputs: dict) -> None:
    """`categorized` is never a child, so its priority is never compared."""
    prio = tmp_path / "complete.priority.txt"
    prio.write_text("LINE 1 categorized\nSINE 2 categorized\nbackground 3 categorized\n")
    spec = BuildSpec.from_flags(
        db_id="HKS_root_exempt",
        version="1.0.0",
        sequence=tiny_inputs["genome"],
        feature_beds={"repeat": tiny_inputs["repeat_bed"]},
        priorities={"repeat": prio},
        s=11,
    )
    # Preparation must get past the priority check; it may still fail later if
    # hks is absent, which is what @requires_hks guards elsewhere.
    try:
        build_database(spec, tmp_path / "db", register=False)
    except BuildError as e:
        assert "have no entry" not in str(e)
