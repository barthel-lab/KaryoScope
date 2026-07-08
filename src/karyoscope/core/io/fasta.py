"""Minimal FASTA reader / writer and DNA helpers.

KaryoScope only needs three operations on FASTA files:

* Read every sequence into memory keyed by its header name. Done by
  :func:`read_fasta_records`. Returns an :class:`collections.OrderedDict`
  so the writer can iterate in source order when needed.
* Reverse-complement a DNA string for the orientation step. Done by
  :func:`reverse_complement`. Handles the IUPAC ambiguity codes in
  addition to the canonical ``ACGT``, and preserves case (lower
  stays lower) so users who keep soft-masked sequences get the same
  masking back after a reverse complement.
* Write a dict of sequences back out. Done by
  :func:`write_fasta_records`. Plain or gzip output, with optional
  line wrapping (default: unwrapped, one sequence per line, matching
  what the archive produced).

This module deliberately does not implement indexed (samtools-faidx)
random access. The scaffold use case loads each input FASTA exactly
once and emits a renamed/oriented copy; in-memory is fast enough
for any current real-world genome assembly (~3 GB worst case for
human) and avoids a samtools dependency at scaffold time.
"""

from __future__ import annotations

import gzip
from collections import OrderedDict
from pathlib import Path
from typing import IO

#: IUPAC complement table. Includes lowercase so soft-masked
#: sequence remains soft-masked after reverse-complementing.
_COMPLEMENT = str.maketrans(
    "ACGTUMRWSYKVHDBNacgtumrwsykvhdbn",
    "TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn",
)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of ``seq``.

    Recognises ACGT, IUPAC ambiguity codes (RYSWKMBDHVN), and U
    (RNA, mapped to A). Unknown characters pass through unchanged
    rather than raising, so the function works on annotated
    sequences containing gap characters (``-``) or other markers.
    """
    return seq.translate(_COMPLEMENT)[::-1]


def _open_in(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def _open_out(path: Path, *, gzip_out: bool) -> IO[str]:
    if gzip_out:
        return gzip.open(path, "wt")
    return path.open("w")


def read_fasta_records(path: Path) -> OrderedDict[str, str]:
    """Read a FASTA file into ``{name: sequence}``.

    The name is the first whitespace-delimited token of the ``>``
    header line. Sequences are concatenated across multi-line FASTA
    bodies. Insertion order is preserved.

    Plain and ``.gz`` inputs are both supported. The whole file is
    held in memory; for the FASTA mode of ``karyoscope scaffold``
    this is fine (the assembly is read once per input).
    """
    records: OrderedDict[str, str] = OrderedDict()
    current_name: str | None = None
    current_chunks: list[str] = []

    with _open_in(path) as h:
        for raw in h:
            line = raw.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records[current_name] = "".join(current_chunks)
                head = line[1:].lstrip()
                current_name = head.split()[0] if head else ""
                current_chunks = []
            else:
                if current_name is not None:
                    current_chunks.append(line)

    if current_name is not None:
        records[current_name] = "".join(current_chunks)
    return records


def read_fasta_contig_names(path: Path) -> list[str]:
    """Return just the contig names from ``path``.

    Lighter-weight alternative to :func:`read_fasta_records` when the
    caller doesn't need sequence bodies. The reader stops at the
    end of each header line and skips sequence lines entirely.
    """
    names: list[str] = []
    with _open_in(path) as h:
        for line in h:
            if line.startswith(">"):
                head = line[1:].lstrip()
                if head:
                    names.append(head.split()[0])
    return names


def read_fasta_lengths(path: Path) -> OrderedDict[str, int]:
    """Read a FASTA file into ``{name: sequence_length}`` in source order.

    Streaming counterpart to :func:`read_fasta_records` for callers that
    only need lengths (e.g. the combined-scaffold layout): sequence
    bodies are counted, never held, so memory stays O(number of
    contigs). Name parsing and the length computation match
    :func:`read_fasta_records` exactly (first whitespace token; blank
    lines skipped; concatenated body length).
    """
    lengths: OrderedDict[str, int] = OrderedDict()
    current_name: str | None = None
    current_len = 0

    with _open_in(path) as h:
        for raw in h:
            line = raw.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    lengths[current_name] = current_len
                head = line[1:].lstrip()
                current_name = head.split()[0] if head else ""
                current_len = 0
            elif current_name is not None:
                current_len += len(line)

    if current_name is not None:
        lengths[current_name] = current_len
    return lengths


def write_fasta_records(
    records: dict[str, str] | OrderedDict[str, str],
    path: Path,
    *,
    gzip_out: bool | None = None,
    line_width: int | None = None,
) -> None:
    """Write ``records`` to ``path``.

    ``gzip_out=None`` (the default) gzips iff ``path`` ends in
    ``.gz``. ``line_width=None`` writes each sequence on a single
    line (the archive's convention); pass an integer to wrap at
    that width.

    Iteration order follows ``records`` -- callers that need a
    specific output ordering should pass an :class:`OrderedDict`
    or rely on the ``dict`` insertion-order guarantee.
    """
    if gzip_out is None:
        gzip_out = str(path).endswith(".gz")

    with _open_out(path, gzip_out=gzip_out) as out:
        for name, seq in records.items():
            out.write(f">{name}\n")
            if line_width is None or line_width <= 0:
                out.write(seq)
                out.write("\n")
            else:
                for i in range(0, len(seq), line_width):
                    out.write(seq[i : i + line_width])
                    out.write("\n")
