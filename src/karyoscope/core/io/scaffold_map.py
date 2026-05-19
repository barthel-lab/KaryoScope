"""The authoritative ``scaffold_map.tsv`` format.

After ``karyoscope scaffold`` classifies and orients each contig, it
writes a per-input map TSV that pins every encoded contig name back to
its source. The encoded name (``<chrom>_<hap>_<contig>[_rc]``) is a
human-readable convenience; the map file is the contract.

Every consumer downstream of ``scaffold`` (the BED rewriter in this
stage, ``centromeres`` and ``karyotype`` in later stages, future FASTA
mode in Stage 5d-1b) reads ``scaffold_map.tsv`` to learn the new -> old
mapping rather than parsing the encoded name. This makes the encoded
format itself a presentation choice that we can change between
releases without breaking the pipeline.

Schema (tab-separated, header row required, columns in this order):

    new_name        chr1_hap1_h1tg000001l_rc
    original_name   h1tg000001l
    input_file      hap1.fa.gz                 (basename, not path)
    hap             hap1
    chromosome      chr1
    flipped         yes | no
    length          248956422
    stats           TPCQT                      (legacy P/C/Q/T summary)

A separate, smaller ``scaffold_stats.tsv`` is written alongside in the
archive's 2-column format (``new_name\\tstats``) for back-compat with
any script that consumed the old output. See
:func:`write_legacy_stats`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from karyoscope.exceptions import ScaffoldError

_HEADER: tuple[str, ...] = (
    "new_name",
    "original_name",
    "input_file",
    "hap",
    "chromosome",
    "flipped",
    "length",
    "stats",
)


@dataclass(frozen=True)
class MapRow:
    """One scaffolded contig's source-of-truth record.

    ``flipped`` is a real ``bool``; serialisation maps it to ``yes`` /
    ``no`` in the TSV so a human can scan the file.

    ``input_file`` is the *basename* of the source FASTA (not the full
    path) so the map file is portable across machines.
    """

    new_name: str
    original_name: str
    input_file: str
    hap: str
    chromosome: str
    flipped: bool
    length: int
    stats: str


def write_map(rows: list[MapRow], path: Path) -> None:
    """Write a list of :class:`MapRow` records to ``path`` as TSV."""
    with path.open("w") as h:
        h.write("\t".join(_HEADER) + "\n")
        for r in rows:
            h.write(
                "\t".join(
                    (
                        r.new_name,
                        r.original_name,
                        r.input_file,
                        r.hap,
                        r.chromosome,
                        "yes" if r.flipped else "no",
                        str(r.length),
                        r.stats,
                    )
                )
                + "\n"
            )


def read_map(path: Path) -> list[MapRow]:
    """Read a ``scaffold_map.tsv`` previously written by :func:`write_map`.

    Validates header, column count per row, and the ``flipped`` /
    ``length`` types. Raises :class:`ScaffoldError` on any malformation.
    """
    if not path.is_file():
        raise ScaffoldError(f"scaffold map not found: {path}")

    rows: list[MapRow] = []
    with path.open("r") as h:
        header_line = h.readline().rstrip("\n")
        header = tuple(header_line.split("\t"))
        if header != _HEADER:
            raise ScaffoldError(f"{path}: unexpected header. Got {header!r}, expected {_HEADER!r}")
        for i, raw in enumerate(h, start=2):
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != len(_HEADER):
                raise ScaffoldError(
                    f"{path}:{i}: expected {len(_HEADER)} columns, got {len(parts)}"
                )
            flipped_str = parts[5].strip().lower()
            if flipped_str not in ("yes", "no"):
                raise ScaffoldError(
                    f"{path}:{i}: 'flipped' must be 'yes' or 'no', got {parts[5]!r}"
                )
            try:
                length = int(parts[6])
            except ValueError as e:
                raise ScaffoldError(
                    f"{path}:{i}: 'length' must be an integer, got {parts[6]!r}"
                ) from e
            rows.append(
                MapRow(
                    new_name=parts[0],
                    original_name=parts[1],
                    input_file=parts[2],
                    hap=parts[3],
                    chromosome=parts[4],
                    flipped=(flipped_str == "yes"),
                    length=length,
                    stats=parts[7],
                )
            )
    return rows


def write_legacy_stats(rows: list[MapRow], path: Path) -> None:
    """Write the archive's 2-column ``<name>\\t<stats>\\n`` format.

    Kept for back-compat with any pre-CLI script that consumed
    ``scaffold_stats.tsv`` directly. The richer information lives in
    the map file.
    """
    with path.open("w") as h:
        for r in rows:
            h.write(f"{r.new_name}\t{r.stats}\n")
