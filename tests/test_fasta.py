"""Unit tests for :mod:`karyoscope.core.io.fasta`."""

from __future__ import annotations

import gzip
from collections import OrderedDict
from pathlib import Path

from karyoscope.core.io.fasta import (
    read_fasta_contig_names,
    read_fasta_records,
    reverse_complement,
    write_fasta_records,
)

# --- reverse_complement ---------------------------------------------


class TestReverseComplement:
    def test_basic_acgt(self) -> None:
        assert reverse_complement("AATTCG") == "CGAATT"

    def test_case_preserved(self) -> None:
        # Lowercase stays lowercase (soft-masked sequence round-trips).
        assert reverse_complement("aattCG") == "CGaatt"

    def test_iupac_ambiguity_codes(self) -> None:
        # Per the IUPAC table: R<->Y, S<->S, W<->W, K<->M, B<->V, D<->H, N<->N.
        assert reverse_complement("RYSWKMBDHVN") == "NBDHVKMWSRY"

    def test_palindrome_returns_self(self) -> None:
        # GGCC is its own reverse complement.
        assert reverse_complement("GGCC") == "GGCC"

    def test_empty_string(self) -> None:
        assert reverse_complement("") == ""

    def test_rna_u_becomes_a(self) -> None:
        # Per the table, U complements to A (one-way for RNA->DNA convenience).
        assert reverse_complement("UU") == "AA"

    def test_unknown_chars_pass_through(self) -> None:
        # Gap characters and other markers don't error.
        assert reverse_complement("A-T") == "A-T"  # A<->T flip + gap stays


# --- read_fasta_records --------------------------------------------


class TestReadFastaRecords:
    def test_single_record(self, tmp_path: Path) -> None:
        p = tmp_path / "x.fa"
        p.write_text(">seqA\nACGT\n")
        assert read_fasta_records(p) == OrderedDict([("seqA", "ACGT")])

    def test_multi_record_preserves_order(self, tmp_path: Path) -> None:
        p = tmp_path / "x.fa"
        p.write_text(">a\nAC\n>b\nGT\n>c\nTT\n")
        recs = read_fasta_records(p)
        assert list(recs.keys()) == ["a", "b", "c"]

    def test_multiline_sequence_joined(self, tmp_path: Path) -> None:
        p = tmp_path / "x.fa"
        p.write_text(">a\nACGT\nNNNN\nGGGG\n")
        assert read_fasta_records(p)["a"] == "ACGTNNNNGGGG"

    def test_header_takes_first_token_only(self, tmp_path: Path) -> None:
        p = tmp_path / "x.fa"
        p.write_text(">seqA some description line\nACGT\n")
        assert "seqA" in read_fasta_records(p)

    def test_gzip(self, tmp_path: Path) -> None:
        p = tmp_path / "x.fa.gz"
        with gzip.open(p, "wt") as h:
            h.write(">a\nACGT\n>b\nGCTA\n")
        recs = read_fasta_records(p)
        assert recs == OrderedDict([("a", "ACGT"), ("b", "GCTA")])

    def test_strips_carriage_returns(self, tmp_path: Path) -> None:
        # Windows-style line endings.
        p = tmp_path / "x.fa"
        p.write_text(">a\r\nACGT\r\n")
        assert read_fasta_records(p)["a"] == "ACGT"


# --- read_fasta_contig_names ---------------------------------------


class TestReadFastaContigNames:
    def test_returns_names_only(self, tmp_path: Path) -> None:
        p = tmp_path / "x.fa"
        p.write_text(">a desc\nACGT\n>b\nGCTA\n")
        assert read_fasta_contig_names(p) == ["a", "b"]


# --- write_fasta_records -------------------------------------------


class TestWriteFastaRecords:
    def test_plain(self, tmp_path: Path) -> None:
        p = tmp_path / "x.fa"
        write_fasta_records({"a": "ACGT", "b": "GCTA"}, p)
        assert p.read_text() == ">a\nACGT\n>b\nGCTA\n"

    def test_gzip_when_path_ends_gz(self, tmp_path: Path) -> None:
        p = tmp_path / "x.fa.gz"
        write_fasta_records({"a": "ACGT"}, p)
        with gzip.open(p, "rt") as h:
            assert h.read() == ">a\nACGT\n"

    def test_explicit_gzip_false_overrides_extension(self, tmp_path: Path) -> None:
        # Path ends in .gz but we say "don't gzip" -- the file is plain text.
        p = tmp_path / "x.fa.gz"
        write_fasta_records({"a": "ACGT"}, p, gzip_out=False)
        assert p.read_bytes().startswith(b">a\n")

    def test_line_wrap(self, tmp_path: Path) -> None:
        p = tmp_path / "x.fa"
        write_fasta_records({"a": "ACGTACGTAC"}, p, line_width=4)
        # 10-char sequence wrapped at 4: "ACGT", "ACGT", "AC".
        assert p.read_text() == ">a\nACGT\nACGT\nAC\n"

    def test_iteration_order(self, tmp_path: Path) -> None:
        p = tmp_path / "x.fa"
        write_fasta_records(OrderedDict([("z", "A"), ("a", "C"), ("m", "G")]), p)
        # Should iterate in insertion order, not alphabetical.
        assert p.read_text() == ">z\nA\n>a\nC\n>m\nG\n"


# --- round-trip -----------------------------------------------------


class TestRoundTrip:
    def test_roundtrip_plain(self, tmp_path: Path) -> None:
        src = {"a": "ACGT", "b": "NNN"}
        p = tmp_path / "x.fa"
        write_fasta_records(src, p)
        assert dict(read_fasta_records(p)) == src

    def test_roundtrip_gzip(self, tmp_path: Path) -> None:
        src = {"a": "ACGTACGT", "b": "GCATGCAT"}
        p = tmp_path / "x.fa.gz"
        write_fasta_records(src, p)
        assert dict(read_fasta_records(p)) == src
