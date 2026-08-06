"""``prep-bed asat`` — per-array alpha-satellite feature set from a CenSat annotation."""

from __future__ import annotations

import pytest

from karyoscope.core.prep import asat
from karyoscope.core.prep.common import PrepError

CENSAT = "\n".join(
    [
        # a plain single-array HOR interval
        "chr1\t100\t200\thor_1_5(S1C1/5/19H1L)\t100\t.\t100\t200\t255,146,0",
        # bare continuation suffixes: three arrays, not leaves called B and C
        "chr1\t200\t260\thor_1_1(S3C1H2-A,B,C)",
        # two genuinely different arrays on one interval
        "chr16\t300\t340\thor_16_3(S1C16H1L,S1C1/5/19H1L)",
        # a live array followed by a bare suffix -> stem drops the trailing L
        "chr21\t400\t410\thor_21_3(S2C13/21H1L,B)",
        # divergent HOR, both a named array and an SF-only label
        "chr11\t500\t560\tdhor_11_2(S3C11H3d)",
        "chr13\t600\t640\tdhor_13_1(SF2)",
        # monomeric alpha-satellite carries no parenthetical
        "chr1\t700\t760\tmon_1_2",
        # not alpha-satellite -- must not appear
        "chr1\t800\t900\tbsat_1_1",
        "chr1\t900\t950\tct_1_2",
        "chr1\t950\t990\thsat3_1_1",
    ]
)


@pytest.fixture
def censat_bed(tmp_path):
    p = tmp_path / "censat.bed"
    p.write_text(CENSAT + "\n")
    return p


def run(censat_bed, tmp_path, **kw):
    out, hier = tmp_path / "asat.bed", tmp_path / "asat.tsv"
    result = asat.from_censat(input_path=censat_bed, output=out, hierarchy=hier, **kw)
    records = [line.split("\t") for line in out.read_text().splitlines()]
    edges = dict(line.split("\t") for line in hier.read_text().splitlines())
    return result, records, edges


# -- label parsing ----------------------------------------------------


@pytest.mark.parametrize(
    ("names", "expected"),
    [
        (["S3C1H2-A", "B", "C"], ["S3C1H2-A", "S3C1H2-B", "S3C1H2-C"]),
        (["S2C16H2-B", "A"], ["S2C16H2-B", "S2C16H2-A"]),
        # the live-array 'L' is a variant token too, so the stem drops it
        (["S2C13/21H1L", "B"], ["S2C13/21H1L", "S2C13/21H1-B"]),
        # two full names are left alone
        (["S1C16H1L", "S1C1/5/19H1L"], ["S1C16H1L", "S1C1/5/19H1L"]),
        # a suffix with no preceding stem cannot be expanded, so it stands
        (["B"], ["B"]),
    ],
)
def test_expand_labels(names, expected):
    assert asat.expand_labels(names) == expected


def test_class_detection_prefers_dhor_over_hor():
    assert asat.censat_class("dhor_11_2(S3C11H3d)") == "dhor"
    assert asat.censat_class("hor_1_5(S1C1/5/19H1L)") == "hor"
    assert asat.censat_class("mon_X_3") == "mon"
    assert asat.censat_class("hor_Y_1(S4CYH1L)") == "hor"
    assert asat.censat_class("bsat_13_1") is None
    assert asat.censat_class("hsat3_1_1") is None


def test_slash_becomes_underscore():
    """`build` names one FASTA per leaf, so a '/' would read as a path separator."""
    assert asat.safe_label("S1C1/5/19H1L") == "S1C1_5_19H1L"
    assert asat.parse_labels("hor_1_5(S1C1/5/19H1L)", "hor") == ["S1C1_5_19H1L"]


def test_interval_without_parenthetical_is_labelled_by_class():
    assert asat.parse_labels("mon_1_2", "mon") == ["mon"]


# -- record emission --------------------------------------------------


def test_multi_label_interval_emits_one_record_per_label(censat_bed, tmp_path):
    """No winner is picked: HKS resolves the shared k-mers to their ancestor."""
    _result, records, _edges = run(censat_bed, tmp_path)
    shared = [r for r in records if r[0] == "chr16" and r[1] == "300"]
    assert sorted(r[3] for r in shared) == ["S1C16H1L", "S1C1_5_19H1L"]
    # ... and both keep the full interval, rather than splitting it
    assert {(r[1], r[2]) for r in shared} == {("300", "340")}


def test_bare_suffixes_never_become_leaves(censat_bed, tmp_path):
    _result, records, edges = run(censat_bed, tmp_path)
    labels = {r[3] for r in records}
    assert {"S3C1H2-A", "S3C1H2-B", "S3C1H2-C"} <= labels
    assert not [lab for lab in labels if len(lab) == 1]
    assert not [lab for lab in edges if len(lab) == 1]


def test_non_alpha_satellite_is_excluded(censat_bed, tmp_path):
    """bSat/HSat/ct are other satellites; they belong to the region set, not here."""
    _result, records, _edges = run(censat_bed, tmp_path)
    assert all(r[0] != "chr1" or int(r[1]) < 800 for r in records)
    # 7 alpha-sat intervals, of which two name 3 and 2 arrays -> 7 + 2 + 1 + 1 = 11
    assert len(records) == 11


def test_records_are_coordinate_sorted(censat_bed, tmp_path):
    _result, records, _edges = run(censat_bed, tmp_path)
    keys = [(r[0], int(r[1]), int(r[2]), r[3]) for r in records]
    assert keys == sorted(keys)


def test_union_bp_does_not_double_count_multi_label_intervals(censat_bed, tmp_path):
    result, _records, _edges = run(censat_bed, tmp_path)
    # 100 + 60 + 40 + 10 + 60 + 40 + 60, each interval counted once
    assert "370 bp" in result.notes[0]


# -- hierarchy --------------------------------------------------------


def test_scaffold_nests_every_class_under_asat(censat_bed, tmp_path):
    """The whole point: mon and dhor sit inside asat, not in the background.

    The shape mirrors the shipped `region` set, so the two feature sets in one
    database describe the same biology with the same names.
    """
    _result, _records, edges = run(censat_bed, tmp_path)
    assert edges["asat"] == "categorized"
    assert edges["alpha_hor"] == "asat"
    assert edges["active_hor"] == "alpha_hor"
    assert edges["hor"] == "alpha_hor"
    assert edges["dhor"] == "alpha_hor"
    assert edges["mon"] == "asat"


def test_arrays_hang_off_their_class(censat_bed, tmp_path):
    _result, _records, edges = run(censat_bed, tmp_path)
    assert edges["S1C1_5_19H1L"] == "active_hor"  # trailing L = a live array
    assert edges["S3C1H2-B"] == "hor"  # an inactive higher-order repeat
    assert edges["SF2"] == "dhor"


def test_label_named_in_both_classes_is_placed_by_dhor(tmp_path):
    """A divergent array named inside a hor interval is still divergent."""
    src = tmp_path / "c.bed"
    src.write_text(
        "chr11\t0\t10\thor_11_2(S3C11H1L,S3C11H3d)\nchr11\t20\t30\tdhor_11_2(S3C11H3d)\n"
    )
    _result, _records, edges = run(src, tmp_path)
    assert edges["S3C11H3d"] == "dhor"
    assert edges["S3C11H1L"] == "active_hor"


def test_every_emitted_label_appears_in_the_hierarchy(censat_bed, tmp_path):
    """build rejects a produced label that the hierarchy does not name."""
    _result, records, edges = run(censat_bed, tmp_path)
    assert {r[3] for r in records} <= set(edges)


# -- class selection --------------------------------------------------


def test_excluding_a_class_drops_it_and_warns(censat_bed, tmp_path):
    result, records, edges = run(censat_bed, tmp_path, classes=("hor",))
    assert "mon" not in {r[3] for r in records}
    assert "dhor" not in edges
    assert any("resolve to the root and paint nothing" in n for n in result.notes)


def test_selecting_active_hor_works_in_either_dialect(tmp_path):
    """Selection is on the label's resolved class, not the record's.

    In the older dialect a live array arrives inside a `hor` record and only
    becomes `active_hor` once the `L` rule is applied, so filtering on the
    record's class would make `--class active_hor` match nothing there.
    """
    legacy = tmp_path / "legacy.bed"
    legacy.write_text("chr1\t100\t200\thor_1_5(S1C1/5/19H1L)\nchr1\t200\t260\thor_1_1(S3C1H2-A)\n")
    canonical = tmp_path / "canon.bed"
    canonical.write_text(
        "chr1\t100\t200\tactive_hor(S1C1/5/19H1L)\nchr1\t200\t260\thor(S3C1H2-A)\n"
    )
    for src in (legacy, canonical):
        _result, records, edges = run(src, tmp_path, classes=("active_hor",))
        assert {r[3] for r in records} == {"S1C1_5_19H1L"}
        assert edges["S1C1_5_19H1L"] == "active_hor"
        assert "S3C1H2-A" not in edges  # inactive, not selected


def test_unknown_class_is_rejected(censat_bed, tmp_path):
    with pytest.raises(PrepError, match="unknown"):
        run(censat_bed, tmp_path, classes=("hsat3",))


def test_no_alpha_satellite_records_is_an_error(tmp_path):
    src = tmp_path / "empty.bed"
    src.write_text("chr1\t0\t10\tbsat_1_1\n")
    with pytest.raises(PrepError, match="no alpha-satellite records"):
        run(src, tmp_path)


def test_result_declares_a_background(censat_bed, tmp_path):
    """The non-satellite genome still needs a gap-fill; it just has no alpha-sat left in it."""
    result, _records, _edges = run(censat_bed, tmp_path)
    assert result.background == "background"


# -- priorities -------------------------------------------------------


def prio_rows(tmp_path):
    return {
        line.split("\t")[0]: (int(line.split("\t")[1]), line.split("\t")[2])
        for line in (tmp_path / "asat.prio.txt").read_text().splitlines()
    }


def test_priority_ranks_the_two_sibling_groups(censat_bed, tmp_path):
    """Live beats inactive beats divergent; any HOR beats monomeric.

    Priorities are only compared within a sibling group, so these are two
    independent orderings rather than one global ranking.
    """
    run(censat_bed, tmp_path, priority=tmp_path / "asat.prio.txt")
    rows = prio_rows(tmp_path)
    assert rows["alpha_hor"][0] < rows["mon"][0]  # under asat
    assert rows["active_hor"][0] < rows["hor"][0] < rows["dhor"][0]  # under alpha_hor


def test_priority_names_the_background_explicitly(censat_bed, tmp_path):
    """build defaults an unlisted node to 0 -- the BEST priority.

    Omitting the gap-fill would therefore make it outrank every feature, and an
    array sharing all its k-mers with the background would stop painting.
    """
    run(censat_bed, tmp_path, priority=tmp_path / "asat.prio.txt")
    rows = prio_rows(tmp_path)
    assert "background" in rows
    assert rows["background"] == (asat.ROOT_PRIORITY, "categorized")


def test_background_and_asat_tie_so_shared_kmers_fall_back_to_the_root(censat_bed, tmp_path):
    """Neither side should claim a k-mer shared with the rest of the genome."""
    run(censat_bed, tmp_path, priority=tmp_path / "asat.prio.txt")
    rows = prio_rows(tmp_path)
    assert rows["background"][0] == rows["asat"][0]


def test_priority_honours_a_renamed_background(censat_bed, tmp_path):
    run(censat_bed, tmp_path, priority=tmp_path / "asat.prio.txt", background="nonsatellite")
    rows = prio_rows(tmp_path)
    assert "nonsatellite" in rows and "background" not in rows


def test_priority_siblings_are_all_equal_or_all_distinct(censat_bed, tmp_path):
    """HKS rejects a mixed sibling group; validate the same rule build does."""
    from karyoscope.core.io.hks import validate_sibling_priorities

    run(censat_bed, tmp_path, priority=tmp_path / "asat.prio.txt")
    rows = prio_rows(tmp_path)
    parent_of = {child: parent for child, (_p, parent) in rows.items()}
    priority = {child: p for child, (p, _parent) in rows.items()}
    assert validate_sibling_priorities(parent_of, priority) == []


def test_leaves_share_one_priority(censat_bed, tmp_path):
    """No basis exists for ranking arrays against each other, so they tie."""
    run(censat_bed, tmp_path, priority=tmp_path / "asat.prio.txt")
    rows = prio_rows(tmp_path)
    classes = set(asat.CLASSES) | {asat.ALPHA_HOR, asat.ASAT}
    leaves = [p for child, (p, _parent) in rows.items() if child not in classes]
    assert leaves and len(set(leaves)) == 1


def test_no_priority_file_written_unless_asked(censat_bed, tmp_path):
    result, _records, _edges = run(censat_bed, tmp_path)
    assert result.priority is None
    assert not (tmp_path / "asat.prio.txt").exists()
