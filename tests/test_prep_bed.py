"""Tests for ``karyoscope prep-bed`` and the converters behind it.

The converters were checked against the real reference feature sets during
development — the EDTA, ``.fai``, satellite and UCSC-cytoband outputs are
byte-identical to the files the bespoke scripts produced from the same inputs.
What is pinned here is the behaviour that regression would silently break:
coordinate conventions, the precedence rule, which rows survive, and the fact
that the pasteable stanza goes to stdout while commentary goes to stderr.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest
from click.testing import CliRunner

from karyoscope.cli import main
from karyoscope.core.prep import cytoband as cytoband_prep
from karyoscope.core.prep import genes as genes_prep
from karyoscope.core.prep import repeats as repeats_prep
from karyoscope.core.prep import structural as structural_prep
from karyoscope.core.prep.common import (
    PrepError,
    PrepResult,
    SeqidRewriter,
    coalesce,
    open_text,
    read_fai,
    render_stanza,
)


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def fai(tmp_path: Path) -> Path:
    """A .fai deliberately NOT in natural order, to catch re-sorting."""
    path = tmp_path / "ref.fa.fai"
    path.write_text("chr2\t2000\t0\t60\t61\nchr1\t1000\t0\t60\t61\nchrM\t100\t0\t60\t61\n")
    return path


def bed_rows(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text().splitlines()]


# -- shared plumbing --------------------------------------------------


def test_open_text_detects_gzip_by_content_not_suffix(tmp_path: Path) -> None:
    """A gzipped file with a misleading name still reads: we sniff magic bytes."""
    path = tmp_path / "annotation.gff3"  # no .gz, but gzipped
    path.write_bytes(gzip.compress(b"chr1\tx\texon\t1\t10\t.\t+\t.\tID=a\n"))
    with open_text(path) as fh:
        assert fh.read().startswith("chr1")


def test_read_fai_preserves_file_order(fai: Path) -> None:
    assert list(read_fai(fai)) == ["chr2", "chr1", "chrM"]


def test_read_fai_rejects_a_non_integer_length(tmp_path: Path) -> None:
    path = tmp_path / "bad.fai"
    path.write_text("chr1\tlots\n")
    with pytest.raises(PrepError, match="not an integer"):
        read_fai(path)


def test_seqid_rewriter_table_beats_prefix(tmp_path: Path) -> None:
    table = tmp_path / "map.tsv"
    table.write_text("chrA_1\trenamed\n")
    rewriter = SeqidRewriter.build(rename_prefix="chrA_:Chr", seqid_map=table)
    assert rewriter("chrA_1") == "renamed"  # exact table entry wins
    assert rewriter("chrA_2") == "Chr2"  # falls through to the prefix rule
    assert rewriter("scaffold_9") == "scaffold_9"  # matched by neither


def test_seqid_rewriter_rejects_a_prefix_without_a_colon() -> None:
    with pytest.raises(PrepError, match="OLD:NEW"):
        SeqidRewriter.build(rename_prefix="nocolon", seqid_map=None)


def test_coalesce_joins_only_abutting_same_label_runs() -> None:
    assert coalesce([("c", 0, 5, "a"), ("c", 5, 9, "a"), ("c", 9, 12, "b"), ("d", 0, 4, "a")]) == [
        ("c", 0, 9, "a"),
        ("c", 9, 12, "b"),
        ("d", 0, 4, "a"),
    ]


def test_render_stanza_marks_a_fully_tiled_set_and_lists_exclusions(tmp_path: Path) -> None:
    stanza = render_stanza(
        PrepResult(
            name="cytoband",
            bed=Path("cytoband.bed"),
            n_records=3,
            background=None,
            hierarchy=Path("cytoband.tsv"),
            exclude=["chrM", "chrUn_x"],
        )
    )
    assert "background: null" in stanza
    assert "exclude:\n  - chrM\n  - chrUn_x" in stanza


def test_render_stanza_names_the_background_when_there_is_one() -> None:
    stanza = render_stanza(
        PrepResult(name="repeat", bed=Path("r.bed"), n_records=1, background="nonrepeat")
    )
    assert "background: nonrepeat" in stanza
    assert "exclude:" not in stanza


# -- RepeatMasker -----------------------------------------------------

RM_OUT = """\
   SW   perc perc perc  query     position in query      matching  repeat  position in repeat
score   div. del. ins.  sequence  begin end (left) strand repeat    class/family  begin end (left)  ID

  463   1.3  0.6  1.7  chr1        11     60  (940) +   L1MA        LINE/L1         1   50  (0)  1
  200  10.0  0.0  0.0  chr1       101    150  (850) C   AluSx       SINE/Alu        1   50  (0)  2
  150  12.0  0.0  0.0  chr1       201    250  (750) +   (CA)n       Simple_repeat   1   50  (0)  3
  120  15.0  0.0  0.0  chr1       301    350  (650) +   Weird       Novel_class     1   50  (0)  4
"""

RM_BED = "\n".join(
    [
        "chr1\t10\t60\tL1MA#LINE/L1\t.\t+\t10\t60\t0\t1\t50\t0",
        "chr1\t100\t150\tAluSx#SINE/Alu\t.\t-\t100\t150\t0\t1\t50\t0",
        "chr1\t200\t250\t(CA)n#Simple_repeat\t.\t+\t200\t250\t0\t1\t50\t0",
        "chr1\t300\t350\tWeird#Novel_class\t.\t+\t300\t350\t0\t1\t50\t0",
    ]
)


def _run_rm(tmp_path: Path, text: str, name: str, **kwargs: object) -> tuple[PrepResult, Path]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    source = tmp_path / name
    source.write_text(text + "\n")
    output = tmp_path / "repeat.bed"
    result = repeats_prep.from_repeatmasker(
        input_path=source,
        output=output,
        hierarchy=tmp_path / "repeat.tsv",
        **kwargs,  # type: ignore[arg-type]
    )
    return result, output


def test_repeatmasker_out_dialect_converts_1_based_inclusive_to_bed(tmp_path: Path) -> None:
    result, output = _run_rm(tmp_path, RM_OUT, "rm.out")
    assert "'out' dialect" in result.notes[0]
    assert bed_rows(output)[0] == ["chr1", "10", "60", "LINE"]


def test_repeatmasker_bed_dialect_gives_the_same_labels(tmp_path: Path) -> None:
    """Both dialects of the same tool's output must agree, coordinates included."""
    _, from_out = _run_rm(tmp_path / "a", RM_OUT, "rm.out")
    _, from_bed = _run_rm(tmp_path / "b", RM_BED, "rm.bed")
    assert from_out.read_text() == from_bed.read_text()


def test_repeatmasker_labels_unknown_classes_rather_than_dropping_them(tmp_path: Path) -> None:
    """Dropping them would leave those bases to the nonrepeat gap-fill, asserting
    the opposite of what the annotation says."""
    result, output = _run_rm(tmp_path, RM_OUT, "rm.out")
    labels = [row[3] for row in bed_rows(output)]
    assert labels == ["LINE", "SINE", "Simple_repeat", "Unknown"]
    assert any("Novel_class" in note for note in result.notes)


def test_repeatmasker_strips_the_uncertainty_marker(tmp_path: Path) -> None:
    """RepeatMasker's ``DNA?`` is still DNA, not an unrecognised class."""
    text = RM_OUT + "  100  1.0  0.0  0.0  chr1  401 450 (550) + X  DNA?/hAT  1 50 (0) 5\n"
    _, output = _run_rm(tmp_path, text, "rm.out")
    assert bed_rows(output)[-1][3] == "DNA"


def test_repeatmasker_colors_group_the_rna_leaves_into_one_legend_row(tmp_path: Path) -> None:
    colors = tmp_path / "colors.tsv"
    _run_rm(tmp_path, RM_OUT, "rm.out", colors=colors)
    rows = [line.split("\t") for line in colors.read_text().splitlines()[1:]]
    groups = {feature: group for _set, feature, _color, group in rows}
    assert groups["rRNA"] == groups["tRNA"] == "RNA"
    assert groups["LINE"] == ""  # a colour of its own, so no grouping needed
    # build enforces "same group => same colour"; make sure we honour it.
    by_group: dict[str, set[str]] = {}
    for _set, _feature, color, group in rows:
        if group:
            by_group.setdefault(group, set()).add(color)
    assert all(len(colors_used) == 1 for colors_used in by_group.values())


def test_repeatmasker_hierarchy_covers_every_leaf_it_can_emit(tmp_path: Path) -> None:
    """A leaf missing from the hierarchy fails the build, so check it here."""
    nodes = {child for child, _parent in repeats_prep.RM_HIERARCHY}
    assert set(repeats_prep.RM_CLASS_TO_LEAF.values()) <= nodes
    assert repeats_prep.RM_FALLBACK_LEAF in nodes
    assert set(repeats_prep.RM_COLORS) == nodes


# -- EDTA -------------------------------------------------------------

EDTA_GFF = """\
##gff-version 3
Chr1\tEDTA\trepeat_region\t101\t200\t.\t+\t.\tID=TE_1;Classification=DNA/HAT
Chr1\tEDTA\trepeat_region\t301\t400\t.\t+\t.\tID=TE_2;Classification=DNA/DTA
Chr1\tEDTA\trepeat_region\t501\t600\t.\t+\t.\tID=TE_3;Classification=Nonsense/Made-up
"""


def test_edta_normalises_wicker_aliases_onto_one_leaf(tmp_path: Path) -> None:
    """DNA/HAT and DNA/DTA are both hAT; they must not become two leaves."""
    source = tmp_path / "te.gff3"
    source.write_text(EDTA_GFF)
    output = tmp_path / "repeat.bed"
    result = repeats_prep.from_edta(
        input_path=source, output=output, hierarchy=tmp_path / "repeat.tsv"
    )
    rows = bed_rows(output)
    assert [row[3] for row in rows[:2]] == ["hAT", "hAT"]
    assert rows[0][1:3] == ["100", "200"]  # GFF3 1-based inclusive -> BED half-open
    assert rows[2][3] == repeats_prep.EDTA_FALLBACK_LEAF
    assert any("Nonsense/Made-up" in note for note in result.notes)
    assert result.background == "nonrepeat"


def test_edta_hierarchy_covers_every_alias_target() -> None:
    nodes = {child for child, _parent in repeats_prep.EDTA_HIERARCHY}
    assert set(repeats_prep.EDTA_ALIAS.values()) <= nodes


# -- gene -------------------------------------------------------------

GTF = """\
chr1\tsrc\texon\t101\t200\t.\t+\t.\ttranscript_id "t1"; gene_id "g1";
chr1\tsrc\texon\t301\t400\t.\t+\t.\ttranscript_id "t1"; gene_id "g1";
"""

GFF3 = """\
##gff-version 3
chr1\tsrc\texon\t101\t200\t.\t+\t.\tID=e1;Parent=transcript:t1
chr1\tsrc\texon\t301\t400\t.\t+\t.\tID=e2;Parent=transcript:t1
"""


def _run_gene(tmp_path: Path, text: str, fai: Path, **kwargs: object) -> Path:
    source = tmp_path / "genes.gff3"
    source.write_text(text)
    output = tmp_path / "gene.bed"
    genes_prep.from_gff(
        input_path=source,
        lengths=fai,
        output=output,
        hierarchy=tmp_path / "gene.tsv",
        **kwargs,  # type: ignore[arg-type]
    )
    return output


def test_transcript_ids_reads_both_attribute_dialects() -> None:
    assert genes_prep.transcript_ids("ID=e1;Parent=transcript:t1") == ["transcript:t1"]
    assert genes_prep.transcript_ids('transcript_id "t1"; gene_id "g1";') == ["t1"]
    assert genes_prep.transcript_ids("Parent=a,b") == ["a", "b"]
    assert genes_prep.transcript_ids("nothing_useful") == []


@pytest.mark.parametrize("text", [GTF, GFF3], ids=["gtf", "gff3"])
def test_gene_derives_introns_between_a_transcripts_exons(
    tmp_path: Path, fai: Path, text: str
) -> None:
    output = _run_gene(tmp_path, text, fai)
    chr1 = [row for row in bed_rows(output) if row[0] == "chr1"]
    assert chr1 == [
        ["chr1", "0", "100", "intergenic"],
        ["chr1", "100", "200", "exon"],
        ["chr1", "200", "300", "intron"],
        ["chr1", "300", "400", "exon"],
        ["chr1", "400", "1000", "intergenic"],
    ]


def test_gene_does_not_call_the_gap_between_two_genes_an_intron(tmp_path: Path, fai: Path) -> None:
    """The gap between separate single-exon transcripts is intergenic. Deriving
    introns from a chromosome-wide exon list instead gets this wrong, and it is
    the difference between ~20% and ~42% intron on a compact genome."""
    text = (
        'chr1\tsrc\texon\t101\t200\t.\t+\t.\ttranscript_id "t1";\n'
        'chr1\tsrc\texon\t301\t400\t.\t+\t.\ttranscript_id "t2";\n'
    )
    output = _run_gene(tmp_path, text, fai)
    labels = {(row[1], row[2]): row[3] for row in bed_rows(output) if row[0] == "chr1"}
    assert labels[("200", "300")] == "intergenic"


def test_gene_prefers_exon_over_intron_where_transcripts_disagree(
    tmp_path: Path, fai: Path
) -> None:
    """One transcript's intron overlapping another's exon must read as exon."""
    text = (
        'chr1\tsrc\texon\t101\t200\t.\t+\t.\ttranscript_id "t1";\n'
        'chr1\tsrc\texon\t401\t500\t.\t+\t.\ttranscript_id "t1";\n'
        'chr1\tsrc\texon\t251\t300\t.\t+\t.\ttranscript_id "t2";\n'
    )
    output = _run_gene(tmp_path, text, fai)
    labels = {(row[1], row[2]): row[3] for row in bed_rows(output) if row[0] == "chr1"}
    assert labels[("250", "300")] == "exon"
    assert labels[("200", "250")] == "intron"


def test_gene_tiles_every_sequence_including_unannotated_ones(tmp_path: Path, fai: Path) -> None:
    output = _run_gene(tmp_path, GTF, fai)
    covered: dict[str, int] = {}
    for seqid, start, end, _label in bed_rows(output):
        assert covered.setdefault(seqid, 0) == int(start), f"{seqid} is not contiguous"
        covered[seqid] = int(end)
    assert covered == read_fai(fai)


def test_gene_merges_abutting_spans_of_the_same_label(tmp_path: Path, fai: Path) -> None:
    """Overlapping exons must yield one record, not one per boundary."""
    text = (
        'chr1\tsrc\texon\t101\t300\t.\t+\t.\ttranscript_id "t1";\n'
        'chr1\tsrc\texon\t201\t400\t.\t+\t.\ttranscript_id "t2";\n'
    )
    output = _run_gene(tmp_path, text, fai)
    exons = [row for row in bed_rows(output) if row[3] == "exon"]
    assert exons == [["chr1", "100", "400", "exon"]]


def test_gene_reports_a_useful_error_when_nothing_matches(tmp_path: Path, fai: Path) -> None:
    text = 'other\tsrc\texon\t101\t200\t.\t+\t.\ttranscript_id "t1";\n'
    with pytest.raises(PrepError, match="--seqid-map"):
        _run_gene(tmp_path, text, fai)


# -- cytoband ---------------------------------------------------------

CYTOBAND = "\n".join(
    [
        "chr1\t0\t1000\tp36.33\tgneg",
        "chr1\t1000\t1500\tp36.32\tgpos25",
        "chr1\t1500\t2000\tp33\tacen",
        "chr2\t0\t2000\tq11.1\tgvar",
        "chrUn_x\t0\t50\tp11\tgneg",
    ]
)


def _run_cytoband(tmp_path: Path, fai: Path, **kwargs: object) -> tuple[PrepResult, Path]:
    source = tmp_path / "cytoBand.txt"
    source.write_text(CYTOBAND + "\n")
    output = tmp_path / "cytoband.bed"
    result = cytoband_prep.from_ucsc(
        input_path=source,
        lengths=fai,
        output=output,
        hierarchy=tmp_path / "cytoband.tsv",
        **kwargs,  # type: ignore[arg-type]
    )
    return result, output


def test_cytoband_labels_keep_their_chromosome(tmp_path: Path, fai: Path) -> None:
    """Bare band names collide across chromosomes; ``1p36.33`` cannot."""
    _, output = _run_cytoband(tmp_path, fai)
    assert [row[3] for row in bed_rows(output)] == ["2q11.1", "1p36.33", "1p36.32", "1p33"]


def test_cytoband_nests_sub_banded_and_singleton_bands_differently(
    tmp_path: Path, fai: Path
) -> None:
    result, _ = _run_cytoband(tmp_path, fai)
    edges = dict(line.split("\t") for line in (result.hierarchy or Path()).read_text().splitlines())
    assert edges["1p36.33"] == "1p36"  # sub-banded: under its group
    assert edges["1p36"] == "1"
    assert edges["1p33"] == "1"  # singleton: straight under the chromosome
    assert edges["1"] == "categorized"


def test_cytoband_reports_unbanded_sequences_for_exclude_rather_than_labelling_them(
    tmp_path: Path, fai: Path
) -> None:
    """The archive wrote chrM as a literal 'exclude' label; build has a real
    exclude: list now, and a placeholder label would claim the sequence."""
    result, output = _run_cytoband(tmp_path, fai)
    assert "chrM" in result.exclude
    assert all(row[0] != "chrM" for row in bed_rows(output))


def test_cytoband_colors_come_from_the_stain_and_group_the_legend(
    tmp_path: Path, fai: Path
) -> None:
    colors = tmp_path / "colors.tsv"
    _run_cytoband(tmp_path, fai, colors=colors)
    lines = colors.read_text().splitlines()
    assert lines[0].endswith("legend_group")
    rows = {parts[1]: (parts[2], parts[3]) for parts in (line.split("\t") for line in lines[1:])}
    assert rows["1p36.33"] == (cytoband_prep.STAIN_COLORS["gneg"], "gneg")
    assert rows["1p33"] == (cytoband_prep.STAIN_COLORS["acen"], "acen")
    # Interior nodes are not bands and must not join a stain group.
    assert rows["1"] == (cytoband_prep.INTERIOR_COLOR, cytoband_prep.INTERIOR_GROUP)


def test_cytoband_flags_a_stain_it_has_no_colour_for(tmp_path: Path, fai: Path) -> None:
    source = tmp_path / "cytoBand.txt"
    source.write_text("chr1\t0\t1000\tp36.33\tglowing\n")
    result = cytoband_prep.from_ucsc(
        input_path=source,
        lengths=fai,
        output=tmp_path / "c.bed",
        hierarchy=tmp_path / "c.tsv",
        colors=tmp_path / "c_colors.tsv",
    )
    assert any("glowing" in note for note in result.notes)


def test_cytoband_refuses_a_pattern_that_matches_nothing(tmp_path: Path, fai: Path) -> None:
    with pytest.raises(PrepError, match="widen --primary-pattern"):
        _run_cytoband(tmp_path, fai, primary_pattern=r"^nope$")


# -- chromosome / satellite -------------------------------------------


def test_fai_emits_one_whole_length_record_per_sequence_in_file_order(
    tmp_path: Path, fai: Path
) -> None:
    output = tmp_path / "chromosome.bed"
    result = structural_prep.from_fai(lengths=fai, output=output)
    assert bed_rows(output) == [
        ["chr2", "0", "2000", "chr2"],
        ["chr1", "0", "1000", "chr1"],
        ["chrM", "0", "100", "chrM"],
    ]
    assert result.background is None
    assert result.hierarchy is None  # grouping is curation, not derivable


def test_core_cluster_picks_the_densest_cluster_not_the_widest_extent() -> None:
    """A stray distal monomer must not drag the centromere across the arm."""
    bands = [(0, 10), (1_000_000, 1_010_000), (1_010_500, 1_020_000)]
    assert structural_prep.core_cluster(bands, cluster_gap=500_000) == (1_000_000, 1_020_000)


def test_merge_spans_bridges_only_gaps_within_the_limit() -> None:
    assert structural_prep.merge_spans([(0, 10), (11, 20), (100, 110)], gap=10) == [
        (0, 20),
        (100, 110),
    ]
    assert structural_prep.merge_spans([(0, 10), (11, 20)], gap=0) == [(0, 10), (11, 20)]


def test_satellite_splits_the_arms_around_the_core_and_tiles_fully(
    tmp_path: Path, fai: Path
) -> None:
    source = tmp_path / "sat.gff3"
    source.write_text("chr1\ttrf\tmonomer\t501\t600\t.\t+\t.\tID=m1\n")
    output = tmp_path / "region.bed"
    result = structural_prep.from_satellite(
        input_path=source,
        lengths=fai,
        output=output,
        hierarchy=tmp_path / "region.tsv",
        satellite="CEN180",
    )
    chr1 = [row for row in bed_rows(output) if row[0] == "chr1"]
    assert chr1 == [
        ["chr1", "0", "500", "p_arm"],
        ["chr1", "500", "600", "CEN180"],
        ["chr1", "600", "1000", "q_arm"],
    ]
    # Sequences with no satellite are still covered, so the set tiles the genome.
    assert [row for row in bed_rows(output) if row[0] == "chrM"] == [["chrM", "0", "100", "p_arm"]]
    assert result.background is None


def test_satellite_reads_bed_input_as_zero_based(tmp_path: Path, fai: Path) -> None:
    """The coordinate convention follows the suffix; a BED must not be shifted."""
    source = tmp_path / "sat.bed"
    source.write_text("chr1\t500\t600\tmonomer\n")
    output = tmp_path / "region.bed"
    structural_prep.from_satellite(
        input_path=source, lengths=fai, output=output, hierarchy=tmp_path / "region.tsv"
    )
    assert ["chr1", "500", "600", "satellite"] in bed_rows(output)


# -- CLI wiring -------------------------------------------------------


def test_stanza_goes_to_stdout_and_commentary_to_stderr(
    tmp_path: Path, fai: Path, runner: CliRunner
) -> None:
    """So `prep-bed ... >> spec.yaml` captures the stanza and nothing else."""
    result = runner.invoke(
        main,
        ["prep-bed", "fai", "--lengths", str(fai), "--output", str(tmp_path / "c.bed")],
        catch_exceptions=False,
    )
    assert result.exit_code == 0
    assert result.stdout.startswith("# add to your build spec:")
    assert "feature_sets:\n  - name: chromosome" in result.stdout
    assert "wrote" not in result.stdout


def test_existing_outputs_are_not_clobbered_without_force(
    tmp_path: Path, fai: Path, runner: CliRunner
) -> None:
    existing = tmp_path / "c.bed"
    existing.write_text("precious\n")
    args = ["prep-bed", "fai", "--lengths", str(fai), "--output", str(existing)]
    result = runner.invoke(main, args)
    assert result.exit_code != 0
    assert "refusing to overwrite" in result.output
    assert existing.read_text() == "precious\n"

    assert runner.invoke(main, [*args, "--force"]).exit_code == 0
    assert existing.read_text() != "precious\n"


def test_every_format_subcommand_is_reachable(runner: CliRunner) -> None:
    listing = runner.invoke(main, ["prep-bed", "--help"]).output
    for leaf in ("repeatmasker", "edta", "gff-gene", "cytoband", "fai", "satellite"):
        assert leaf in listing


def test_a_bad_rename_prefix_is_a_clean_error_not_a_traceback(
    tmp_path: Path, fai: Path, runner: CliRunner
) -> None:
    result = runner.invoke(
        main,
        [
            "prep-bed",
            "fai",
            "--lengths",
            str(fai),
            "--output",
            str(tmp_path / "c.bed"),
            "--rename-prefix",
            "missing-colon",
        ],
    )
    assert result.exit_code != 0
    assert "OLD:NEW" in result.output
    assert "Traceback" not in result.output
