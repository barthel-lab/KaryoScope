"""``prep-bed pave`` — a papillomavirus `gene` set from PaVE genome records."""

from __future__ import annotations

import json

import pytest

from karyoscope.core.prep import pave
from karyoscope.core.prep.common import PrepError, SeqidRewriter


def feature(name, locations, ftype="CDS"):
    return {
        "featureName": name,
        "featureType": ftype,
        "featureIsSpliced": len(locations) > 1,
        "featureLocations": [{"start": s, "end": e} for s, e in locations],
    }


def genome(pave_id, species, genus, sequence, features):
    return {
        "genome": {"pave_id": pave_id, "description": f"{pave_id}, complete genome"},
        "itemSequence": {"sequence": sequence, "length": len(sequence)},
        "virus_lineage": [
            {"name": pave_id},
            {"name": species},
            {"name": genus},
            {"name": "Papillomaviridae"},
        ],
        "features": features,
    }


#: Two genomes exercising the cases that matter: a spliced transcript nested in
#: its parent ORF, a binding motif, a URR spanning the origin, an E1/E2 overlap,
#: a bare ``E5``, and a gap no ORF covers.
RECORDS = [
    genome(
        "HPV16REF",
        "Alphapapillomavirus 9",
        "Alphapapillomavirus",
        "ACGT" * 50,  # 200 bp
        [
            feature("E2BS", [(1, 12)], ftype="BindingMotif"),
            feature("E6", [(11, 40)]),
            feature("E6*", [(11, 20), (30, 40)]),  # nested in E6
            feature("E7", [(38, 60)]),  # overlaps E6 by 3
            feature("E1", [(60, 120)]),
            feature("E1^E4", [(60, 70), (100, 120)]),  # nested in E1
            feature("E2", [(115, 150)]),  # overlaps E1 by 6
            feature("L2", [(155, 175)]),  # leaves 150..155 uncovered
            feature("L1", [(175, 190)]),
            feature("URR", [(190, 200), (1, 10)]),  # spans the origin
        ],
    ),
    genome(
        "HPV6REF",
        "Alphapapillomavirus 10",
        "Alphapapillomavirus",
        "TTGC" * 25,  # 100 bp
        [
            feature("E5", [(1, 30)]),  # the bare form
            feature("E5_ALPHA", [(40, 60)]),
            feature("URR", [(60, 100)]),
        ],
    ),
]


@pytest.fixture
def records_json(tmp_path):
    p = tmp_path / "pave.json"
    p.write_text(json.dumps(RECORDS))
    return p


def run(records_json, tmp_path, **kw):
    out, hier = tmp_path / "gene.bed", tmp_path / "gene.tsv"
    result = pave.from_pave(input_path=records_json, output=out, hierarchy=hier, **kw)
    records = [line.split("\t") for line in out.read_text().splitlines()]
    edges = dict(line.split("\t") for line in hier.read_text().splitlines())
    return result, records, edges


def test_spliced_transcripts_and_binding_motifs_are_dropped(records_json, tmp_path):
    """They are wholly nested in their parent ORF / shorter than any usable k."""
    _result, records, _edges = run(records_json, tmp_path)
    labels = {r[3] for r in records}
    assert not labels & {"E1^E4", "E8^E2", "E6*", "E2BS", "E1BS"}
    assert {"E6", "E7", "E1", "E2", "L2", "L1", "URR"} <= labels


def test_origin_spanning_feature_becomes_two_records(records_json, tmp_path):
    """A circular genome's URR arrives as two locations and stays two intervals."""
    _result, records, _edges = run(records_json, tmp_path)
    urr = [r for r in records if r[0] == "HPV16REF" and r[3] == "URR"]
    assert sorted((int(r[1]), int(r[2])) for r in urr) == [(0, 10), (189, 200)]


def test_coordinates_convert_to_bed_half_open(records_json, tmp_path):
    """PaVE is 1-based inclusive; BED is 0-based half-open."""
    _result, records, _edges = run(records_json, tmp_path)
    e6 = next(r for r in records if r[0] == "HPV16REF" and r[3] == "E6")
    assert (int(e6[1]), int(e6[2])) == (10, 40)


def test_overlaps_are_kept_for_priority_to_resolve(records_json, tmp_path):
    """The set is built in priority mode, so overlapping ORFs are not flattened."""
    _result, records, _edges = run(records_json, tmp_path)
    e1 = next(r for r in records if r[0] == "HPV16REF" and r[3] == "E1")
    e2 = next(r for r in records if r[0] == "HPV16REF" and r[3] == "E2")
    assert int(e2[1]) < int(e1[2])  # they still overlap


def test_bare_e5_is_renamed_so_leaf_and_node_names_stay_distinct(records_json, tmp_path):
    _result, records, edges = run(records_json, tmp_path)
    assert any(r[3] == pave.E5_OTHER for r in records)
    assert not any(r[3] == "E5" for r in records)
    assert edges[pave.E5_OTHER] == "E5"  # E5 is the interior node
    assert edges["E5"] == "early"


def test_hierarchy_groups_early_and_late(records_json, tmp_path):
    _result, _records, edges = run(records_json, tmp_path)
    assert edges["E6"] == edges["E1"] == "early"
    assert edges["L1"] == edges["L2"] == "late"
    assert edges["early"] == edges["late"] == edges["URR"] == "categorized"


def test_priority_file_names_every_node_including_the_background(records_json, tmp_path):
    """`build` rejects a priority file with a gap, and the gap-fill leaf is the
    one the hierarchy file does not carry."""
    prio = tmp_path / "gene.priority.txt"
    result, _records, edges = run(records_json, tmp_path, priority=prio)
    rows = [line.split("\t") for line in prio.read_text().splitlines()]
    named = {r[0] for r in rows}
    assert named == set(edges) | {"intergenic"}
    assert result.n_priority == len(rows)
    values = {r[0]: int(r[1]) for r in rows}
    # Genome order: the earlier ORF wins every overlap it takes part in.
    assert values["E6"] < values["E7"] < values["E1"] < values["E2"]
    assert values["L2"] < values["L1"]
    assert values["early"] < values["late"] < values["URR"] < values["intergenic"]


def test_siblings_have_distinct_priorities(records_json, tmp_path):
    """HKS requires siblings all-equal or all-distinct; ties would resolve to
    the parent instead of naming a gene."""
    prio = tmp_path / "gene.priority.txt"
    run(records_json, tmp_path, priority=prio)
    rows = [line.split("\t") for line in prio.read_text().splitlines()]
    by_parent: dict[str, list[int]] = {}
    for _child, value, parent in rows:
        by_parent.setdefault(parent, []).append(int(value))
    for parent, values in by_parent.items():
        assert len(values) == len(set(values)), f"tied siblings under {parent}"


def test_taxonomy_hierarchy_drops_the_constant_family_level(records_json, tmp_path):
    tax = tmp_path / "type.hierarchy.txt"
    run(records_json, tmp_path, taxonomy=tax)
    edges = dict(line.split("\t") for line in tax.read_text().splitlines())
    assert "Papillomaviridae" not in edges  # a single child of the root says nothing
    assert edges["HPV16REF"] == "Alphapapillomavirus_9"
    assert edges["Alphapapillomavirus_9"] == "Alphapapillomavirus"
    assert edges["Alphapapillomavirus"] == "categorized"


def test_taxon_names_lose_their_spaces(records_json, tmp_path):
    """Hierarchy files are whitespace-separated, so 'Alphapapillomavirus 9'
    would parse as two columns and be silently truncated."""
    tax = tmp_path / "type.hierarchy.txt"
    run(records_json, tmp_path, taxonomy=tax)
    for line in tax.read_text().splitlines():
        assert len(line.split()) == 2


def test_fasta_carries_the_sequences_and_matches_the_bed_seqids(records_json, tmp_path):
    fasta = tmp_path / "hpv.fasta"
    _result, records, _edges = run(records_json, tmp_path, fasta=fasta)
    names = [ln[1:].split()[0] for ln in fasta.read_text().splitlines() if ln.startswith(">")]
    assert names == ["HPV16REF", "HPV6REF"]
    assert {r[0] for r in records} <= set(names)
    body = "".join(ln for ln in fasta.read_text().splitlines() if not ln.startswith(">"))
    assert len(body) == 300  # 200 + 100


def test_uncovered_bases_are_reported_not_filled(records_json, tmp_path):
    """`build` owns the gap-fill; the converter only names the label."""
    result, records, _edges = run(records_json, tmp_path)
    assert result.background == "intergenic"
    assert not any(r[3] == "intergenic" for r in records)
    assert any("uncovered" in note for note in result.notes)


def test_seqid_rewriting_reaches_bed_fasta_and_taxonomy(records_json, tmp_path):
    fasta, tax = tmp_path / "hpv.fasta", tmp_path / "type.hierarchy.txt"
    _result, records, _edges = run(
        records_json,
        tmp_path,
        fasta=fasta,
        taxonomy=tax,
        rename=SeqidRewriter(table={"HPV16REF": "HPV16"}),
    )
    assert "HPV16" in {r[0] for r in records}
    assert "HPV16REF" not in {r[0] for r in records}
    assert ">HPV16 " in fasta.read_text()
    assert "HPV16\tAlphapapillomavirus_9" in tax.read_text()


def test_unknown_feature_name_is_an_error_not_a_silent_drop(tmp_path):
    """A new PaVE feature needs a leaf, a priority and a colour before it can be
    built; dropping it would quietly hand its bases to the gap-fill."""
    rec = genome("HPVxREF", "Sp 1", "Alphapapillomavirus", "ACGT" * 10, [feature("E99", [(1, 20)])])
    p = tmp_path / "pave.json"
    p.write_text(json.dumps([rec]))
    with pytest.raises(PrepError, match="E99"):
        run(p, tmp_path)


def test_summary_endpoint_records_are_rejected_with_a_usable_message(tmp_path):
    """/api/genome returns summaries with no sequence — a plausible mistake."""
    p = tmp_path / "pave.json"
    p.write_text(json.dumps({"total": 1, "data": [{"locus_id": "HPV16REF", "regions": ["E6"]}]}))
    with pytest.raises(PrepError, match="itemSequence"):
        run(p, tmp_path)


def test_a_single_record_is_accepted(tmp_path):
    """`/api/genome/{id}` returns one object, not a list."""
    p = tmp_path / "pave.json"
    p.write_text(json.dumps(RECORDS[0]))
    _result, records, _edges = run(p, tmp_path)
    assert {r[0] for r in records} == {"HPV16REF"}


def test_colors_group_the_e5_variants_into_one_legend_row(records_json, tmp_path):
    colors = tmp_path / "gene.colors.tsv"
    run(records_json, tmp_path, colors=colors)
    rows = [ln.split("\t") for ln in colors.read_text().splitlines()[1:]]
    groups = {r[1]: r[3] for r in rows}
    assert groups["E5_ALPHA"] == groups["E5_BETA"] == "E5"
    assert groups["E6"] == ""
    hexes = {r[1]: r[2] for r in rows}
    assert hexes["E5_ALPHA"] == hexes["E5_BETA"]  # one colour, so one legend row
