"""PaVE genome records -> a papillomavirus `gene` set, taxonomy, and the FASTA.

The Papillomavirus Episteme (PaVE) publishes one JSON record per genome from
``/api/genome/{id}``, and that record carries three things nothing else has to
be downloaded for: the genome sequence, the ORF coordinates, and the ICTV
lineage. So this converter writes all three — the `gene` feature set it is named
for, plus ``--fasta`` and ``--taxonomy`` outputs that would otherwise each need
a bespoke script. There is no bulk FASTA endpoint that works, so re-deriving
them from the same records is also what keeps a rebuild reproducible.

Two properties of papillomavirus genomes shape the conversion:

**They are circular.** The URR spans the origin, and PaVE expresses that as two
locations (``7157..7906`` and ``1..103``). Both become BED records with the same
label. The k-mers that straddle the junction are not recoverable from a linear
FASTA and are simply absent — about 30 bp per genome.

**The reading frames overlap.** 10.2% of bases are claimed by more than one
annotated feature, so the set is built in priority mode rather than flattened.
Almost all of that overlap is the three *spliced* features, which this converter
drops: E1^E4, E8^E2 and E6* are 100% contained within E1, E2 and E6 across all
224 human reference genomes, so under any ranking that prefers the primary ORF
they would be leaves holding no sequence at all. What remains is 1.2% of bases,
overlapping in genome order, which is what the emitted priority file encodes.
"""

from __future__ import annotations

import json
from collections import defaultdict
from itertools import pairwise
from pathlib import Path
from typing import Any

from karyoscope.core.prep.common import (
    ColorRow,
    Edge,
    PrepError,
    PrepResult,
    Record,
    SeqidRewriter,
    open_text,
    write_bed,
    write_colors,
    write_hierarchy,
)

#: Spliced transcripts, dropped because they are wholly contained in the ORFs
#: they are spliced from and would otherwise be empty leaves. See the module
#: docstring.
SPLICED = frozenset({"E1^E4", "E8^E2", "E6*"})

#: Regulatory binding motifs. At 12-20 bp they are shorter than any usable k,
#: so a k-mer index can never resolve them; including them would only take
#: bases away from the ORFs that contain them.
BINDING_MOTIF = "BindingMotif"

#: PaVE writes one genus-specific E5 per genome plus a bare ``E5`` in one
#: genome. They are grouped under an interior ``E5`` node, so the bare form is
#: renamed to keep leaf and node names distinct.
E5_OTHER = "E5_OTHER"

#: Gene-set tree. Grouping the ORFs into early and late is the standard
#: papillomavirus genome map, and it gives a k-mer shared by two early genes
#: somewhere to resolve to other than the root.
GENE_HIERARCHY: list[Edge] = [
    ("E6", "early"),
    ("E7", "early"),
    ("E1", "early"),
    ("E2", "early"),
    ("E5", "early"),
    ("E10", "early"),
    ("E5_ALPHA", "E5"),
    ("E5_BETA", "E5"),
    ("E5_GAMMA", "E5"),
    ("E5_DELTA", "E5"),
    (E5_OTHER, "E5"),
    ("L2", "late"),
    ("L1", "late"),
    ("early", "categorized"),
    ("late", "categorized"),
    ("URR", "categorized"),
]

#: Priorities, in genome order, which resolves every observed overlap the way
#: the earlier ORF wants: E1 over E2, E7 over E1, E6 over E7, L2 over L1, and
#: (through ``early`` beating ``URR``) E10 over URR. Siblings are distinct, as
#: HKS requires; the background is added by :func:`_priority_rows`.
GENE_PRIORITY: dict[str, int] = {
    "early": 1,
    "late": 2,
    "URR": 3,
    "E6": 1,
    "E7": 2,
    "E1": 3,
    "E2": 4,
    "E5": 5,
    "E10": 6,
    "E5_ALPHA": 1,
    "E5_BETA": 2,
    "E5_GAMMA": 3,
    "E5_DELTA": 4,
    E5_OTHER: 5,
    "L2": 1,
    "L1": 2,
}

#: Early genes cool, late genes warm, regulatory grey — the convention every
#: published papillomavirus genome map uses. The five E5 variants share one
#: colour and one legend row, since only one occurs in any given genome.
GENE_COLORS: dict[str, tuple[str, str]] = {
    "E6": ("#1F77B4", ""),
    "E7": ("#4C9BD4", ""),
    "E1": ("#2E7D32", ""),
    "E2": ("#66A61E", ""),
    "E10": ("#17BECF", ""),
    "E5_ALPHA": ("#9467BD", "E5"),
    "E5_BETA": ("#9467BD", "E5"),
    "E5_GAMMA": ("#9467BD", "E5"),
    "E5_DELTA": ("#9467BD", "E5"),
    E5_OTHER: ("#9467BD", "E5"),
    "L2": ("#D62728", ""),
    "L1": ("#FF7F0E", ""),
    "URR": ("#808080", ""),
}

#: Every genome is a Papillomaviridae, so that level would be a single child of
#: the root carrying no information. The genus becomes the root's child instead.
FAMILY = "Papillomaviridae"

_FASTA_WRAP = 60


def sanitise(label: str) -> str:
    """Make a taxon name usable as a feature label.

    Hierarchy and priority files are whitespace-separated, so a species like
    ``Alphapapillomavirus 9`` would parse as two columns and be silently
    truncated. Collapsing runs of whitespace to a single underscore is enough:
    PaVE taxon names contain no other separator.
    """
    return "_".join(label.split())


def load_records(path: Path) -> list[dict[str, Any]]:
    """Read the PaVE genome records.

    Accepts a JSON array of ``/api/genome/{id}`` responses, a single such
    response, or an object wrapping them under ``data``, so the recipe can
    concatenate downloads the obvious way without a reshaping step.
    """
    with open_text(path) as fh:
        try:
            payload = json.load(fh)
        except json.JSONDecodeError as e:
            raise PrepError(f"{path}: not valid JSON: {e}") from e

    if isinstance(payload, dict):
        payload = payload.get("data", payload)
    records = payload if isinstance(payload, list) else [payload]

    out: list[dict[str, Any]] = []
    for i, rec in enumerate(records):
        if not isinstance(rec, dict) or "itemSequence" not in rec:
            raise PrepError(
                f"{path}: entry {i} is not a PaVE genome record (no 'itemSequence'). "
                "These come from /api/genome/{id}; the /api/genome list endpoint "
                "returns summaries, which carry no sequence or coordinates."
            )
        out.append(rec)
    if not out:
        raise PrepError(f"{path}: no genome records found")
    return out


def genome_id(record: dict[str, Any]) -> str:
    """The record's PaVE locus id, used as the sequence name."""
    for key in ("genome", "item"):
        value = (record.get(key) or {}).get("pave_id")
        if value:
            return str(value)
    raise PrepError("a genome record has no pave_id")


def lineage(record: dict[str, Any]) -> list[str]:
    """Taxon names from the record, most specific first, family dropped."""
    names = [str(t["name"]) for t in record.get("virus_lineage") or [] if t.get("name")]
    return [n for n in names if n != FAMILY]


def _gene_label(feature: dict[str, Any]) -> str | None:
    """Leaf label for a PaVE feature, or ``None`` if it is not part of the set."""
    name = str(feature.get("featureName") or "")
    if not name or feature.get("featureType") == BINDING_MOTIF or name in SPLICED:
        return None
    return E5_OTHER if name == "E5" else name


def _priority_rows(background: str | None) -> list[tuple[str, int, str]]:
    """The tree as ``(child, priority, parent)``, with the gap-fill leaf added.

    ``build`` attaches the background to the root itself, so the hierarchy file
    does not name it — but the priority file has to, or it would default to the
    best priority and outrank every real ORF.
    """
    rows = [(child, GENE_PRIORITY[child], parent) for child, parent in GENE_HIERARCHY]
    if background is not None:
        worst = max(GENE_PRIORITY[c] for c, p in GENE_HIERARCHY if p == "categorized")
        rows.append((background, worst + 1, "categorized"))
    return rows


def write_fasta(path: Path, records: list[dict[str, Any]], rename: SeqidRewriter) -> int:
    """Write the genome sequences carried in the records, returning total bases."""
    total = 0
    with path.open("w") as out:
        for record in records:
            seqid = rename(genome_id(record))
            sequence = str((record.get("itemSequence") or {}).get("sequence") or "").upper()
            if not sequence:
                raise PrepError(f"{seqid}: record carries no sequence")
            description = str((record.get("genome") or {}).get("description") or "").strip()
            out.write(f">{seqid} {description}\n" if description else f">{seqid}\n")
            for i in range(0, len(sequence), _FASTA_WRAP):
                out.write(sequence[i : i + _FASTA_WRAP] + "\n")
            total += len(sequence)
    return total


def write_taxonomy(path: Path, records: list[dict[str, Any]], rename: SeqidRewriter) -> int:
    """Write the ICTV lineage as a hierarchy for the `type` feature set.

    A ``.fai`` cannot supply this — it is exactly the curation `prep-bed fai`
    declines to invent — but PaVE states it, so it is derived rather than
    hand-written.
    """
    edges: list[Edge] = []
    seen: set[Edge] = set()
    for record in records:
        names = [sanitise(n) for n in lineage(record)]
        if not names:
            raise PrepError(f"{genome_id(record)}: record has no virus_lineage")
        # The most specific name is the genome itself; use the sequence name so
        # the leaf matches the `type` BED that `prep-bed fai` writes.
        names[0] = rename(genome_id(record))
        chain = [*names, "categorized"]
        for child, parent in pairwise(chain):
            edge = (child, parent)
            if edge not in seen:
                seen.add(edge)
                edges.append(edge)
    return write_hierarchy(path, edges)


def from_pave(
    *,
    input_path: Path,
    output: Path,
    hierarchy: Path,
    priority: Path | None = None,
    colors: Path | None = None,
    fasta: Path | None = None,
    taxonomy: Path | None = None,
    name: str = "gene",
    background: str = "intergenic",
    rename: SeqidRewriter | None = None,
) -> PrepResult:
    """Convert PaVE genome records into a papillomavirus `gene` feature set."""
    rename = rename or SeqidRewriter()
    records = load_records(input_path)

    bed: list[Record] = []
    per_label: defaultdict[str, int] = defaultdict(int)
    dropped: defaultdict[str, int] = defaultdict(int)
    covered = 0
    total_bases = 0

    for record in records:
        seqid = rename(genome_id(record))
        length = int((record.get("itemSequence") or {}).get("length") or 0)
        if length <= 0:
            raise PrepError(f"{seqid}: record has no sequence length")
        total_bases += length
        claimed: set[int] = set()

        for feature in record.get("features") or []:
            label = _gene_label(feature)
            if label is None:
                dropped[str(feature.get("featureName") or "?")] += 1
                continue
            if label not in GENE_PRIORITY:
                raise PrepError(
                    f"{seqid}: feature {label!r} has no leaf in the papillomavirus gene "
                    "tree. PaVE added a feature name this converter does not know; it "
                    "needs a leaf, a priority and a colour before the set can be built."
                )
            for location in feature.get("featureLocations") or []:
                # PaVE is 1-based inclusive; BED is 0-based half-open. A feature
                # spanning the origin of a circular genome arrives as two
                # locations and becomes two records with the same label.
                start = int(location["start"]) - 1
                end = min(int(location["end"]), length)
                if end <= start:
                    continue
                bed.append((seqid, start, end, label))
                per_label[label] += end - start
                claimed.update(range(start, end))
        covered += len(claimed)

    n_records = write_bed(output, bed)
    n_edges = write_hierarchy(hierarchy, GENE_HIERARCHY)

    n_priority = 0
    if priority is not None:
        rows = _priority_rows(background)
        with priority.open("w") as out:
            for child, value, parent in rows:
                out.write(f"{child}\t{value}\t{parent}\n")
        n_priority = len(rows)

    n_colors = 0
    if colors is not None:
        used = [leaf for leaf, _p in GENE_HIERARCHY if leaf in GENE_COLORS]
        rows_c: list[ColorRow] = [(leaf, *GENE_COLORS[leaf]) for leaf in used]
        rows_c.append((background, "#C7C7C7", ""))
        n_colors = write_colors(colors, name, rows_c)

    notes = [
        f"{len(records)} genome(s), {total_bases:,} bp; "
        f"{total_bases - covered:,} bp ({100 * (total_bases - covered) / total_bases:.2f}%) "
        f"uncovered and left to the '{background}' gap-fill"
    ]
    if dropped:
        shown = ", ".join(f"{k} x{v}" for k, v in sorted(dropped.items()))
        notes.append(f"dropped spliced transcripts and binding motifs: {shown}")
    notes.append(
        "leaf coverage: " + ", ".join(f"{k} {v:,} bp" for k, v in sorted(per_label.items()))
    )

    if fasta is not None:
        written = write_fasta(fasta, records, rename)
        notes.append(f"wrote {fasta} ({written:,} bp) — bgzip and index it for build's `sequence:`")
    if taxonomy is not None:
        n_taxa = write_taxonomy(taxonomy, records, rename)
        notes.append(f"wrote {taxonomy} ({n_taxa:,} edges) — the hierarchy for the `type` set")

    return PrepResult(
        name=name,
        bed=output,
        n_records=n_records,
        background=background,
        hierarchy=hierarchy,
        n_edges=n_edges,
        priority=priority,
        n_priority=n_priority,
        colors=colors,
        n_colors=n_colors,
        notes=notes,
    )
