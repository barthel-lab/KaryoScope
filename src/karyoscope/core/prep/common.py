"""Shared plumbing for the ``prep-bed`` converters.

Readers (:func:`open_text`, :func:`read_fai`), the seqid rewriting every
converter offers, the file writers, and :class:`PrepResult` — the description of
what a converter produced, which :func:`render_stanza` turns into the
build-spec fragment ``prep-bed`` prints on stdout.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO

from karyoscope.exceptions import KaryoscopeError

#: One BED record: ``(seqid, start, end, label)``, start/end 0-based half-open.
Record = tuple[str, int, int, str]

#: One ``child parent`` hierarchy edge.
Edge = tuple[str, str]

#: One colours row: ``(feature, hex_color, legend_group)``; group ``""`` = none.
ColorRow = tuple[str, str, str]


class PrepError(KaryoscopeError):
    """A source annotation could not be converted."""


# -- readers ----------------------------------------------------------


def open_text(path: Path) -> IO[str]:
    """Open ``path`` as text, transparently handling gzip.

    Detects gzip by magic bytes rather than by suffix: annotation downloads are
    routinely renamed, and a mis-suffixed file would otherwise fail deep inside
    a parser with an unreadable error.
    """
    try:
        with path.open("rb") as probe:
            is_gzip = probe.read(2) == b"\x1f\x8b"
    except OSError as e:
        raise PrepError(f"could not read {path}: {e}") from e
    return gzip.open(path, "rt") if is_gzip else path.open()


def read_fai(path: Path) -> dict[str, int]:
    """Read sequence lengths from a samtools ``.fai`` or a 2-column sizes file.

    Insertion order is the file's order, and every converter iterates it as
    given. The ``.fai`` is the assembly's own declaration of what order its
    sequences come in; re-sorting it would invent an ordering the assembly never
    stated, and any rule for doing so is guesswork outside human chromosome
    naming.
    """
    sizes: dict[str, int] = {}
    with path.open() as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise PrepError(f"{path}:{lineno}: expected 'name length', got {raw.rstrip()!r}")
            try:
                sizes[parts[0]] = int(parts[1])
            except ValueError as e:
                raise PrepError(f"{path}:{lineno}: length is not an integer: {parts[1]!r}") from e
    if not sizes:
        raise PrepError(f"{path}: no sequences found")
    return sizes


def iter_gff(path: Path) -> Iterator[tuple[list[str], int]]:
    """Yield ``(fields, lineno)`` for each data line of a GFF3/GTF file.

    Comments, blank lines, and the FASTA section some GFF3 files append are
    skipped; short lines raise, since a truncated annotation should not silently
    produce a partial feature set.
    """
    with open_text(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            if raw.startswith("#"):
                if raw.startswith("##FASTA"):
                    return
                continue
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                raise PrepError(
                    f"{path}:{lineno}: expected 9 tab-separated columns, got {len(fields)}"
                )
            yield fields, lineno


# -- seqid rewriting --------------------------------------------------


@dataclass(frozen=True)
class SeqidRewriter:
    """Rewrites annotation seqids so they match the assembly FASTA.

    Two mechanisms, applied in order: an exact ``old -> new`` table (for
    accession-style names like ``NC_060925.1`` -> ``chr1``), then a prefix
    substitution (for systematic differences like ``Col-CEN_chr`` -> ``Chr``).
    A name matched by neither is passed through unchanged.
    """

    table: dict[str, str] = field(default_factory=dict)
    old_prefix: str = ""
    new_prefix: str = ""

    def __call__(self, seqid: str) -> str:
        mapped = self.table.get(seqid)
        if mapped is not None:
            return mapped
        if self.old_prefix and seqid.startswith(self.old_prefix):
            return self.new_prefix + seqid[len(self.old_prefix) :]
        return seqid

    @classmethod
    def build(cls, *, rename_prefix: str | None, seqid_map: Path | None) -> SeqidRewriter:
        """Assemble a rewriter from the ``--rename-prefix`` / ``--seqid-map`` options."""
        old_prefix = new_prefix = ""
        if rename_prefix:
            if ":" not in rename_prefix:
                raise PrepError(f"--rename-prefix needs OLD:NEW form, got {rename_prefix!r}")
            old_prefix, _, new_prefix = rename_prefix.partition(":")
            if not old_prefix:
                raise PrepError("--rename-prefix: the OLD prefix must not be empty")
        table: dict[str, str] = {}
        if seqid_map is not None:
            with seqid_map.open() as fh:
                for lineno, raw in enumerate(fh, 1):
                    line = raw.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) < 2:
                        raise PrepError(f"{seqid_map}:{lineno}: expected 'old new', got {line!r}")
                    table[parts[0]] = parts[1]
        return cls(table=table, old_prefix=old_prefix, new_prefix=new_prefix)


# -- ordering ---------------------------------------------------------


def coalesce(records: Iterable[Record]) -> list[Record]:
    """Merge abutting records that share a sequence and a label."""
    out: list[list] = []
    for seqid, start, end, label in records:
        if out and out[-1][0] == seqid and out[-1][3] == label and out[-1][2] == start:
            out[-1][2] = end
        else:
            out.append([seqid, start, end, label])
    return [(c, s, e, lbl) for c, s, e, lbl in out]


# -- writers ----------------------------------------------------------


def write_bed(path: Path, records: Iterable[Record]) -> int:
    """Write 4-column BED, returning the record count."""
    n = 0
    with path.open("w") as out:
        for seqid, start, end, label in records:
            if end <= start:
                continue
            out.write(f"{seqid}\t{start}\t{end}\t{label}\n")
            n += 1
    if n == 0:
        raise PrepError(f"{path}: no records produced — check the input and any seqid renaming")
    return n


def write_hierarchy(path: Path, edges: Iterable[Edge]) -> int:
    """Write a 2-column ``child<TAB>parent`` hierarchy, returning the edge count."""
    n = 0
    with path.open("w") as out:
        for child, parent in edges:
            out.write(f"{child}\t{parent}\n")
            n += 1
    return n


def write_priority(
    path: Path,
    edges: Iterable[Edge],
    priority: int = 1,
    background: str | None = None,
    background_priority: int | None = None,
) -> int:
    """Write the same tree as 3-column ``child<TAB>priority<TAB>parent``.

    ``background`` names the gap-fill leaf, which ``build`` adds to the
    hierarchy itself and which therefore does not appear in ``edges``. It has to
    be written anyway: ``build`` requires a priority for every node, and if it
    did not, an unlisted node would take the *best* priority and the gap-fill
    would outrank every real feature.

    It is given ``background_priority`` — by default one step worse than
    ``priority``, so a k-mer a feature shares with unannotated sequence stays
    with the feature.
    """
    n = 0
    with path.open("w") as out:
        for child, parent in edges:
            out.write(f"{child}\t{priority}\t{parent}\n")
            n += 1
        if background is not None:
            worse = priority + 1 if background_priority is None else background_priority
            out.write(f"{background}\t{worse}\tcategorized\n")
            n += 1
    return n


def write_colors(path: Path, set_name: str, rows: Iterable[ColorRow]) -> int:
    """Write a colours file in the 4-column shape ``build`` both reads and emits.

    The ``legend_group`` column is written only when some feature declares one,
    so an ungrouped set stays a plain 3-column file.
    """
    rows = list(rows)
    grouped = any(group for _f, _c, group in rows)
    with path.open("w") as out:
        out.write(
            "feature_set\tfeature\tcolor\tlegend_group\n"
            if grouped
            else "feature_set\tfeature\tcolor\n"
        )
        for feature, color, group in rows:
            if grouped:
                out.write(f"{set_name}\t{feature}\t{color}\t{group}\n")
            else:
                out.write(f"{set_name}\t{feature}\t{color}\n")
    return len(rows)


# -- result and build-spec stanza -------------------------------------


@dataclass
class PrepResult:
    """What a converter produced, and how ``build`` should be told about it."""

    name: str
    bed: Path
    n_records: int
    #: ``None`` means the set already tiles every base, so gap-fill must be off.
    background: str | None = None
    hierarchy: Path | None = None
    n_edges: int = 0
    priority: Path | None = None
    #: Rows in the priority file. Differs from ``n_edges`` when a background is
    #: written, since ``build`` adds that leaf to the hierarchy itself.
    n_priority: int = 0
    colors: Path | None = None
    n_colors: int = 0
    #: Sequences present in the assembly that this set deliberately does not
    #: cover — suggested for the spec's top-level ``exclude:``.
    exclude: list[str] = field(default_factory=list)
    #: Human-readable remarks for stderr (unmapped classes, dropped rows, ...).
    notes: list[str] = field(default_factory=list)

    def summary(self) -> str:
        """One line naming every file written, for stderr."""
        parts = [f"{self.bed} ({self.n_records:,} records)"]
        if self.hierarchy is not None:
            parts.append(f"{self.hierarchy} ({self.n_edges:,} edges)")
        if self.priority is not None:
            parts.append(f"{self.priority} ({self.n_priority or self.n_edges:,} rows)")
        if self.colors is not None:
            parts.append(f"{self.colors} ({self.n_colors:,} features)")
        return "wrote " + ", ".join(parts)


def render_stanza(result: PrepResult) -> str:
    """Render the ``feature_sets:`` fragment to paste into a build spec.

    Hand-rendered rather than dumped through PyYAML so the keys stay in the
    order the docs present them and the paths appear exactly as the user typed
    them.
    """
    lines = [
        "# add to your build spec:",
        "feature_sets:",
        f"  - name: {result.name}",
        f"    bed: {result.bed}",
    ]
    if result.hierarchy is not None:
        lines.append(f"    hierarchy: {result.hierarchy}")
    if result.priority is not None:
        lines.append(f"    priority: {result.priority}")
    if result.colors is not None:
        lines.append(f"    colors: {result.colors}")
    if result.background is None:
        lines.append("    background: null   # this set already tiles every base")
    else:
        lines.append(f"    background: {result.background}")
    if result.exclude:
        lines.append("")
        lines.append("# and at the top level, so no set covers these sequences:")
        lines.append("exclude:")
        lines.extend(f"  - {seqid}" for seqid in result.exclude)
    return "\n".join(lines)
