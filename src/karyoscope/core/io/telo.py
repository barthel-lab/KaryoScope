"""``seqtk telo`` output: parsing existing files, and producing new ones.

``seqtk telo`` reports candidate telomeric repeat regions per sequence
in a 3-column BED-like format::

    seq_a   0       500
    seq_b   100     600
    seq_c   95000   100000

Each row marks one telomere-like region on one sequence. The same
sequence can appear on multiple rows when both ends have telomeric
repeats. A telomere at the very beginning of the sequence (``start ==
0``) is interpreted as a *start telomere*; any other row on the same
sequence is interpreted as a *stop telomere*. ``karyoscope scaffold``
uses only these two binary flags per sequence — the exact coordinates
are irrelevant to the orientation logic.

This module provides:

* :func:`parse_telo_file` — read an existing file and return the
  per-sequence ``{start: bool, stop: bool}`` flags.
* :func:`run_seqtk_telo` — shell out to ``seqtk telo`` when the user
  hasn't provided a precomputed file. Used by the auto-derive cascade
  in :mod:`karyoscope.core.scaffold`.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

from karyoscope.core.external import require_tool, run_tool

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class TeloFlags:
    """Per-sequence start/stop telomere presence.

    ``start`` is True when the sequence has a telomeric region starting
    at coordinate 0; ``stop`` is True when it has any other telomeric
    region (assumed to be at the 3' end).
    """

    start: bool
    stop: bool


def parse_telo_file(path: Path) -> dict[str, TeloFlags]:
    """Parse a ``seqtk telo`` output file.

    Sequences not present in the file are absent from the returned
    dict — callers should treat ``KeyError`` (or :meth:`dict.get`
    returning ``None``) as "no telomeres on either end".

    The file's coordinate fields are validated as integers but their
    exact values are not retained; only the ``start == 0`` boolean is
    consulted.
    """
    if not path.is_file():
        raise FileNotFoundError(f"telo file not found: {path}")

    flags: dict[str, list[bool]] = {}
    # Two-element list per seq: [has_start, has_stop]. Avoids
    # allocating new tuples on every update.

    with path.open("r") as h:
        for i, raw in enumerate(h, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                # seqtk telo always emits at least 3 columns, but be
                # lenient about extra blank or short rows.
                continue
            seq = parts[0]
            try:
                start = int(parts[1])
            except ValueError as e:
                raise ValueError(
                    f"{path}:{i}: 'start' column must be an integer, got {parts[1]!r}"
                ) from e
            cur = flags.setdefault(seq, [False, False])
            if start == 0:
                cur[0] = True
            else:
                cur[1] = True

    return {seq: TeloFlags(start=s, stop=e) for seq, (s, e) in flags.items()}


def run_seqtk_telo(fasta_path: Path, out_path: Path) -> None:
    """Run ``seqtk telo`` on ``fasta_path`` and write its stdout to ``out_path``.

    Raises :class:`karyoscope.core.external.ToolNotFoundError` if seqtk
    isn't on ``$PATH``, and
    :class:`karyoscope.core.external.ExternalToolError` on non-zero
    exit.
    """
    seqtk = require_tool(
        "seqtk",
        install_hint="Install with: conda install -c bioconda seqtk",
    )
    logger.info("running seqtk telo on %s -> %s", fasta_path, out_path)
    result = run_tool([seqtk, "telo", str(fasta_path)], capture=True)
    out_path.write_text(result.stdout or "")
