"""AGP (A Golden Path) writer for combined-chromosome scaffolds.

When ``karyoscope scaffold --combine-chromosomes`` concatenates the
contigs of one ``(chromosome, haplotype)`` into a single sequence, the
simplified record name (``chr4_hap2``) drops the per-contig provenance
that the encoded ``chr4_hap2_<contig>[_rc]`` name carried. The AGP file
restores it: it describes, for every sequence in the output FASTA, the
exact placement of each original contig and the coordinates of the N
gaps inserted between them.

We emit AGP version 2.1. The format is tab-separated, 1-based, with two
row kinds distinguished by column 5 (``component_type``):

* **component rows** (``W`` -- WGS contig)::

      object  obj_beg  obj_end  part_num  W  component_id  comp_beg  comp_end  orientation

  ``comp_beg``/``comp_end`` span the full original contig (``1`` ..
  ``length``); ``orientation`` is ``+`` for a contig used as-is and
  ``-`` for one that scaffold reverse-complemented (the ``_rc``
  contigs).

* **gap rows** (``N`` -- gap of known length)::

      object  obj_beg  obj_end  part_num  N  gap_length  gap_type  linkage  linkage_evidence

  KaryoScope writes ``gap_type=scaffold``, ``linkage=yes``,
  ``linkage_evidence=align_genus``. The evidence value reflects that,
  although KaryoScope is alignment-free, the contig order and
  orientation are asserted from each contig's k-mer feature profile
  against a human reference database (a reference within the same
  genus).

An AGP is meant to reconstruct its FASTA exactly: every sequence and
every base must be accounted for. The writer therefore expects the
caller to pass one :class:`AgpObject` per output FASTA record,
including singleton (uncombined) contigs and any unscaffolded
leftovers, each as a one-component object with no gaps.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

#: Gap metadata KaryoScope writes for every inserted scaffold gap.
#: See the module docstring for the rationale behind ``align_genus``.
GAP_TYPE = "scaffold"
GAP_LINKAGE = "yes"
GAP_LINKAGE_EVIDENCE = "align_genus"


@dataclass(frozen=True)
class AgpComponent:
    """One placed original contig (an AGP ``W`` row).

    ``object_start`` / ``object_end`` are 0-based, half-open, in the
    coordinate system of the combined output sequence; the writer
    converts them to AGP's 1-based inclusive convention. ``length`` is
    the full length of the original contig in bp. ``orientation`` is
    ``"+"`` or ``"-"``.
    """

    component_id: str
    object_start: int
    object_end: int
    length: int
    orientation: str


@dataclass(frozen=True)
class AgpGap:
    """One inserted N gap (an AGP ``N`` row).

    ``object_start`` / ``object_end`` are 0-based half-open in the
    combined sequence; ``length`` equals ``object_end - object_start``
    and is the number of N bases.
    """

    object_start: int
    object_end: int
    length: int


@dataclass(frozen=True)
class AgpObject:
    """One output FASTA record and the parts (components + gaps) it is built from.

    ``parts`` is the ordered list of :class:`AgpComponent` and
    :class:`AgpGap` that tile the object from base 0 to its end with no
    gaps or overlaps. A singleton contig is a single-component object
    with an empty-of-gaps part list of length one.
    """

    name: str
    parts: list[AgpComponent | AgpGap]


def write_agp(objects: list[AgpObject], path: Path) -> None:
    """Write ``objects`` to ``path`` as an AGP 2.1 file.

    Part numbers restart at 1 for each object and increment across both
    component and gap rows (the AGP spec numbers every part, gaps
    included). Coordinates are emitted 1-based inclusive.
    """
    with path.open("w") as h:
        h.write("##agp-version\t2.1\n")
        for obj in objects:
            for part_num, part in enumerate(obj.parts, start=1):
                obj_beg = part.object_start + 1
                obj_end = part.object_end
                if isinstance(part, AgpComponent):
                    h.write(
                        "\t".join(
                            (
                                obj.name,
                                str(obj_beg),
                                str(obj_end),
                                str(part_num),
                                "W",
                                part.component_id,
                                "1",
                                str(part.length),
                                part.orientation,
                            )
                        )
                        + "\n"
                    )
                else:
                    h.write(
                        "\t".join(
                            (
                                obj.name,
                                str(obj_beg),
                                str(obj_end),
                                str(part_num),
                                "N",
                                str(part.length),
                                GAP_TYPE,
                                GAP_LINKAGE,
                                GAP_LINKAGE_EVIDENCE,
                            )
                        )
                        + "\n"
                    )
