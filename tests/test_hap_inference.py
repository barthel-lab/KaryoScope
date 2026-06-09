"""Unit tests for :mod:`karyoscope.core.hap_inference`."""

from __future__ import annotations

from pathlib import Path

from karyoscope.core.hap_inference import (
    assign_per_input_labels,
    classify_contigs,
    infer_hap_from_contig,
    infer_hap_from_filename,
    read_fasta_contig_names,
)

# --- pattern matching -----------------------------------------------


class TestInferFromContig:
    def test_hifiasm_h1(self) -> None:
        assert infer_hap_from_contig("h1tg000001l") == "hap1"

    def test_hifiasm_h2(self) -> None:
        assert infer_hap_from_contig("h2tg000003l") == "hap2"

    def test_explicit_hap1(self) -> None:
        assert infer_hap_from_contig("chr1_hap1_contig") == "hap1"

    def test_explicit_hap2_case_insensitive(self) -> None:
        assert infer_hap_from_contig("chr1.HAP2.001") == "hap2"

    def test_full_word_haplotype1(self) -> None:
        # Some long-read assemblers emit contig names like
        # "haplotype1-0000001"; the bare "hap([12])" rule can't catch
        # these ("hap" is followed by "l"), so they need the full-word
        # pattern.
        assert infer_hap_from_contig("haplotype1-0000001") == "hap1"

    def test_full_word_haplotype2(self) -> None:
        assert infer_hap_from_contig("haplotype2-0000029") == "hap2"

    def test_full_word_haplotype_case_insensitive(self) -> None:
        assert infer_hap_from_contig("HAPLOTYPE2-5") == "hap2"

    def test_full_word_haplotype_does_not_overmatch(self) -> None:
        # "haplotype10" must not be read as hap1 (trailing boundary).
        assert infer_hap_from_contig("haplotype10-0001") is None

    def test_maternal_short(self) -> None:
        assert infer_hap_from_contig("chr1_MAT_contig01") == "maternal"

    def test_maternal_long(self) -> None:
        assert infer_hap_from_contig("HG002.MATERNAL.0001") == "maternal"

    def test_paternal_short(self) -> None:
        assert infer_hap_from_contig("chr1_pat_001") == "paternal"

    def test_paternal_long_lower(self) -> None:
        assert infer_hap_from_contig("hg002.paternal.001") == "paternal"

    def test_no_match(self) -> None:
        assert infer_hap_from_contig("JAHEPM010000001.1") is None

    def test_no_match_for_random_contig(self) -> None:
        assert infer_hap_from_contig("scaffold_001") is None


class TestInferFromFilename:
    def test_simple_stem(self) -> None:
        assert infer_hap_from_filename("hap1") == "hap1"

    def test_stem_with_sample_prefix(self) -> None:
        assert infer_hap_from_filename("HG002.maternal") == "maternal"

    def test_no_match(self) -> None:
        assert infer_hap_from_filename("HG002") is None


# --- assign_per_input_labels ----------------------------------------


class TestAssignPerInputLabels:
    def test_all_explicit_names(self) -> None:
        result = assign_per_input_labels([("hap1", Path("a.fa")), ("hap2", Path("b.fa"))])
        assert result == ["hap1", "hap2"]

    def test_stem_inference(self) -> None:
        result = assign_per_input_labels([(None, Path("hap1.fa.gz")), (None, Path("hap2.fa.gz"))])
        assert result == ["hap1", "hap2"]

    def test_positional_fallback(self) -> None:
        result = assign_per_input_labels(
            [(None, Path("HG002.fa.gz")), (None, Path("unknown.fa.gz"))]
        )
        assert result == ["input1", "input2"]

    def test_mixed_explicit_and_inferred(self) -> None:
        result = assign_per_input_labels(
            [("unassigned", Path("other.fa")), (None, Path("HG002.maternal.fa.gz"))]
        )
        assert result == ["unassigned", "maternal"]

    def test_collision_falls_back_to_positional(self) -> None:
        # Two files both inferring "hap1" → the second gets a positional fallback.
        result = assign_per_input_labels([(None, Path("hap1.fa.gz")), (None, Path("X.hap1.fa.gz"))])
        assert result[0] == "hap1"
        assert result[1] == "input1"


# --- classify_contigs (per input) -----------------------------------


class TestClassifyContigsExplicit:
    def test_explicit_name_wins(self) -> None:
        result = classify_contigs(
            ["h1tg000001l", "h2tg000003l"],
            file_level_label="unassigned",
            explicit_name_given=True,
        )
        assert result == {"h1tg000001l": "unassigned", "h2tg000003l": "unassigned"}


class TestClassifyContigsMultiInput:
    def test_no_split_haps_no_explicit_uses_file_label(self) -> None:
        # In the multi-input case we use the file-level label for
        # everything that doesn't match a built-in pattern.
        result = classify_contigs(
            ["JAHEPM010000001.1", "JAHEPM010000002.1"],
            file_level_label="hap2",
        )
        assert result == {
            "JAHEPM010000001.1": "hap2",
            "JAHEPM010000002.1": "hap2",
        }

    def test_per_contig_pattern_can_override_file_label(self) -> None:
        # When some contigs have a recognisable hap tag, they get
        # that label even though the file-level label is different.
        result = classify_contigs(
            ["random_ctg_001", "h1tg000001l"],
            file_level_label="hap2",
        )
        # random gets the file label; h1tg matches built-in pattern → hap1.
        assert result == {"random_ctg_001": "hap2", "h1tg000001l": "hap1"}


class TestClassifyContigsSingleInput:
    def test_no_matches_all_become_hap1(self) -> None:
        result = classify_contigs(
            ["JAHEPM010000001.1", "JAHEPM010000002.1"],
            file_level_label="hap1",
            is_only_input=True,
        )
        assert result == {
            "JAHEPM010000001.1": "hap1",
            "JAHEPM010000002.1": "hap1",
        }

    def test_one_label_all_get_it(self) -> None:
        result = classify_contigs(
            ["h1tg000001l", "h1tg000002l", "random_ctg"],
            file_level_label="hap1",
            is_only_input=True,
        )
        assert result == {
            "h1tg000001l": "hap1",
            "h1tg000002l": "hap1",
            "random_ctg": "hap1",
        }

    def test_packed_diploid_split(self) -> None:
        # Combined HG002-style file: h1tg / h2tg.
        result = classify_contigs(
            ["h1tg000001l", "h1tg000002l", "h2tg000003l", "h2tg000004l"],
            file_level_label="hap1",  # gets ignored when patterns split
            is_only_input=True,
        )
        assert result == {
            "h1tg000001l": "hap1",
            "h1tg000002l": "hap1",
            "h2tg000003l": "hap2",
            "h2tg000004l": "hap2",
        }

    def test_packed_diploid_split_full_word_haplotype(self) -> None:
        # Combined file whose contigs are named "haplotype1-..." /
        # "haplotype2-..." (the case that previously collapsed to a
        # single hap1, so the karyotype drew no h1/h2 split).
        result = classify_contigs(
            ["haplotype1-0000001", "haplotype1-0000002", "haplotype2-0000029"],
            file_level_label="hap1",
            is_only_input=True,
        )
        assert result == {
            "haplotype1-0000001": "hap1",
            "haplotype1-0000002": "hap1",
            "haplotype2-0000029": "hap2",
        }

    def test_packed_diploid_with_unmatched_gets_default(self) -> None:
        result = classify_contigs(
            ["h1tg000001l", "h2tg000002l", "random_ctg"],
            file_level_label="hap1",
            is_only_input=True,
        )
        # Lexically first matched label is "hap1"; the unmatched
        # contig takes it.
        assert result == {
            "h1tg000001l": "hap1",
            "h2tg000002l": "hap2",
            "random_ctg": "hap1",
        }


class TestClassifyContigsSplitHaps:
    def test_split_haps_regex_captures_label(self) -> None:
        # A custom regex; capture group 1 is the label.
        result = classify_contigs(
            ["chr1_A_01", "chr1_B_02", "chr2_A_01"],
            file_level_label="hap1",
            split_haps_regex=r"_(A|B)_",
        )
        assert result == {
            "chr1_A_01": "A",
            "chr1_B_02": "B",
            "chr2_A_01": "A",
        }

    def test_split_haps_unmatched_get_file_label(self) -> None:
        result = classify_contigs(
            ["chr1_A_01", "weird_thing"],
            file_level_label="hap1",
            split_haps_regex=r"_(A|B)_",
        )
        assert result == {"chr1_A_01": "A", "weird_thing": "hap1"}


# --- read_fasta_contig_names ----------------------------------------


class TestReadFastaContigNames:
    def test_plain_fasta(self, tmp_path: Path) -> None:
        p = tmp_path / "x.fa"
        p.write_text(">seq_a description line\nACGT\n>seq_b\nGCTA\n")
        assert read_fasta_contig_names(p) == ["seq_a", "seq_b"]

    def test_gzip_fasta(self, tmp_path: Path) -> None:
        import gzip

        p = tmp_path / "x.fa.gz"
        with gzip.open(p, "wt") as h:
            h.write(">seq_a\nACGT\n>seq_b\nGCTA\n")
        assert read_fasta_contig_names(p) == ["seq_a", "seq_b"]
